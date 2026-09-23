-- https://github.com/mason-org/mason-registry/blob/main/packages/clice/package.yaml
-- with a newer version. mason-registry pins v0.1.2026082603, which predates
-- `rules.default_command` (clice-io/clice#664).
-- A registry entry rather than `clice@<version>` in ensure_installed: Mason marks a
-- package outdated when its installed version differs from its registry's, so an
-- update would reinstall the older one.
-- Delete once mason-registry's clice is at v0.1.2026092207 or later.

--- One `source.asset` entry: the release archive for a mason target.
--- @param target string mason target, e.g. "darwin_arm64"
--- @param triple string target triple in the archive name
--- @param ext string archive extension
--- @param bin string path of the clice binary inside the archive
--- @return table
local function asset(target, triple, ext, bin)
    return {
        target = target,
        file = ('clice-{{ version | strip_prefix "v" }}.%s.%s'):format(triple, ext),
        bin = bin,
    }
end

return {
    name = "clice",
    description = "A next-generation C++ language server for modern C++, focused on high performance and deep code intelligence.",
    homepage = "https://github.com/clice-io/clice",
    licenses = { "Apache-2.0" },
    languages = { "C", "C++" },
    categories = { "LSP" },
    source = {
        id = "pkg:github/clice-io/clice@v0.1.2026092207",
        asset = {
            asset("darwin_arm64", "aarch64-apple-darwin", "tar.gz", "clice/bin/clice"),
            asset("darwin_x64", "x86_64-apple-darwin", "tar.gz", "clice/bin/clice"),
            asset("linux_arm64_gnu", "aarch64-unknown-linux-gnu", "tar.gz", "clice/bin/clice"),
            asset("linux_x64_gnu", "x86_64-unknown-linux-gnu", "tar.gz", "clice/bin/clice"),
            asset("win_arm64", "aarch64-pc-windows-msvc", "zip", "clice/bin/clice.exe"),
            asset("win_x64", "x86_64-pc-windows-msvc", "zip", "clice/bin/clice.exe"),
        },
    },
    bin = { clice = "{{source.asset.bin}}" },
    neovim = { lspconfig = "clice" },
}
