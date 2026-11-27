#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = ["typed-argument-parser"]
# ///
"""Raw textDocument/hover text from a language server.

Reads a source file and a position in it. Writes the hover's `contents.value`
to stdout, before any float rendering, or the whole result with --json.

The server is launched with the command and settings nvim resolves for it, so
extraPaths, stubPath and whatever mason-lspconfig merged all apply.
"""
import argparse
import json
import os
import subprocess
import sys
import time
from collections.abc import Callable
from typing import Any

from tap import Tap


def nvim_config(server: str) -> dict[str, Any]:
    """Command and settings nvim resolves for a language server.

    Reading them from nvim rather than restating `lsp/<server>.lua` keeps one
    source of truth, and picks up whatever mason-lspconfig merged on top.

    params:
    - server: server name, as in `lsp/<server>.lua`.

    returns: the decoded `{cmd, settings}` table.
    """
    lua = (
        "local c = vim.lsp.config[%r]; "
        "io.stdout:write(vim.json.encode({ cmd = c.cmd, settings = c.settings }))"
    ) % server
    out = subprocess.run(
        ["nvim", "--headless", "-c", "lua " + lua, "-c", "qa!"],
        capture_output=True, text=True, check=True,
    )
    return json.loads(out.stdout)


def poll(query: Callable[[], Any], timeout: float) -> Any:
    """Call a query until it answers with something, or time runs out.

    A server accepts requests before it has indexed, and answers those with
    null rather than making the client wait.

    params:
    - query: called repeatedly; any falsy result counts as not ready.
    - timeout: seconds to keep trying.

    returns: the first truthy result, else None.
    """
    deadline = time.monotonic() + timeout
    while True:
        result = query()
        if result or time.monotonic() > deadline:
            return result
        time.sleep(0.2)


class Client:
    """An LSP session over a server's stdio."""

    def __init__(self, cmd: list[str]) -> None:
        self.proc = subprocess.Popen(
            cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        # Both pipes exist because PIPE was asked for; bind them so the reads
        # and writes below are not each an Optional access.
        assert self.proc.stdin is not None and self.proc.stdout is not None
        self.stdin = self.proc.stdin
        self.stdout = self.proc.stdout
        self.last_id = 0

    def send(self, method: str, params: dict[str, Any], notify: bool = False) -> int | None:
        """Write one request or notification.

        params:
        - method: LSP method name.
        - params: method parameters.
        - notify: send without an id, so there is no reply to wait on.

        returns: the request id, or None for a notification.
        """
        message: dict[str, Any] = {"jsonrpc": "2.0", "method": method, "params": params}
        if not notify:
            self.last_id += 1
            message["id"] = self.last_id
        body = json.dumps(message).encode()
        self.stdin.write(b"Content-Length: %d\r\n\r\n" % len(body) + body)
        self.stdin.flush()
        return message.get("id")

    def receive(self) -> dict[str, Any] | None:
        """Read one message, or None once the server closes its stdout."""
        length = 0
        while True:
            line = self.stdout.readline()
            if not line:
                return None
            line = line.strip()
            if not line:
                break
            name, _, value = line.decode().partition(":")
            if name.lower() == "content-length":
                length = int(value)
        return json.loads(self.stdout.read(length))

    def result(self, request_id: int | None) -> Any:
        """Read until a request's reply arrives, skipping the server's own traffic.

        params:
        - request_id: id returned by `send`.

        returns: the reply's `result`.

        raises:
        - SystemExit: the server exited, or answered with an error.
        """
        while True:
            message = self.receive()
            if message is None:
                sys.exit("server closed the connection")
            if message.get("id") == request_id:
                if "error" in message:
                    sys.exit("server error: %s" % message["error"])
                return message.get("result")

    def initialize(self, root: str, settings: dict[str, Any]) -> None:
        """Complete the handshake and push the workspace configuration.

        params:
        - root: workspace root directory.
        - settings: the server's `workspace/didChangeConfiguration` payload.
        """
        self.result(self.send("initialize", {
            "processId": os.getpid(),
            "rootUri": "file://" + os.path.abspath(root),
            "capabilities": {
                "textDocument": {"hover": {"contentFormat": ["markdown", "plaintext"]}}
            },
        }))
        self.send("initialized", {}, notify=True)
        self.send("workspace/didChangeConfiguration", {"settings": settings}, notify=True)

    def open_file(self, path: str, language: str) -> None:
        """Announce a file's contents to the server.

        params:
        - path: file to read.
        - language: LSP language id.
        """
        with open(path, encoding="utf-8") as handle:
            text = handle.read()
        self.send("textDocument/didOpen", {"textDocument": {
            "uri": "file://" + os.path.abspath(path),
            "languageId": language,
            "version": 1,
            "text": text,
        }}, notify=True)

    def hover(self, path: str, line: int, column: int) -> dict[str, Any] | None:
        """Hover response at a position.

        params:
        - path: an already-opened file.
        - line: 1-based, as an editor displays it.
        - column: 1-based.

        returns: the hover result, or None where the server has nothing to say.
        """
        return self.result(self.send("textDocument/hover", {
            "textDocument": {"uri": "file://" + os.path.abspath(path)},
            "position": {"line": line - 1, "character": column - 1},
        }))


class Args(Tap):
    file: str
    "Source file to hover in."

    line: int
    "1-based line, as displayed."

    column: int
    "1-based column, as displayed."

    server: str = "basedpyright"
    "Name of an lsp/<server>.lua."

    language: str = "python"
    "LSP language id sent with the file."

    root: str = "."
    "Workspace root."

    timeout: float = 20.0
    "Seconds to keep asking while the server indexes."

    json: bool = False
    "Write the whole hover result, not just its text."

    def configure(self) -> None:
        self.add_argument("file")
        self.add_argument("line")
        self.add_argument("column")
        self.add_argument("-s", "--server")
        self.add_argument("-l", "--language")
        self.add_argument("-r", "--root")
        self.add_argument("-t", "--timeout")
        self.add_argument("-j", "--json")


if __name__ == "__main__":
    args = Args(
        underscores_to_dashes=True,
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    ).parse_args()

    config = nvim_config(args.server)
    client = Client(config["cmd"])
    client.initialize(args.root, config["settings"])
    client.open_file(args.file, args.language)
    result = poll(lambda: client.hover(args.file, args.line, args.column), args.timeout)
    client.proc.kill()

    if result is None:
        sys.exit("no hover at %s:%d:%d" % (args.file, args.line, args.column))
    print(json.dumps(result, indent=2) if args.json else result["contents"]["value"])
