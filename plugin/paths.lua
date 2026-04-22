-- append $PATH to vim option &path (&path uses comma delim in text form and $PATH uses colon)
vim.opt_local.path:append(vim.split(vim.env.PATH, ":"))

-- edit-in-kitty on remotes doesn't copy the env variables that are normally
-- present locally so we have to set the necessary ones. These are for julia
-- and python lsp to find the right binaries. These were found by copying all
-- of `env` and checking which fixed the problem by adding them to
-- `edit-in-kitty --env ...` which can be given multiple times like the kitty @
-- launch command.
-- Hardcoding miniforge base env path. It would be cooler to use activated env
-- somehow.
-- ~/.local/bin is location of jedi-language-server
vim.env.PATH = vim.env.HOME
	.. "/.local/bin/:/opt/homebrew/Caskroom/miniforge/base/bin:/opt/homebrew/bin:"
	.. vim.env.PATH

-- PyMOL stubs for LSP are in basedpyright extraPaths (lsp/basedpyright.lua),
-- not PYTHONPATH — global PYTHONPATH breaks the real pymol binary.

-- Git project roots. vim.fs.root walks the tree in pure Lua — no subprocess.
-- `source` is a bufnr or path. git_root matches .git as file or directory, so
-- inside a submodule it returns the submodule root (same as `git rev-parse
-- --show-toplevel`).
local function git_root(source)
	return vim.fs.root(source, ".git")
end

-- Like git_root but filters to .git *directories*, which means it skips
-- submodule gitlinks and returns the outermost enclosing repo.
local function git_root_top(source)
	local path = type(source) == "number" and vim.api.nvim_buf_get_name(source) or source
	path = path ~= "" and vim.fs.dirname(path) or vim.uv.cwd()
	local matches = vim.fs.find(function(name, p)
		if name ~= ".git" then
			return false
		end
		local stat = vim.uv.fs_stat(p .. "/" .. name)
		return stat and stat.type == "directory"
	end, { upward = true, path = path, limit = math.huge })
	return matches[#matches] and vim.fs.dirname(matches[#matches]) or nil
end

-- $ROOT / $ROOTTOP from nvim's startup cwd. Exposed as env vars so shell-style
-- path completion in other plugins (zsh, julia) picks them up.
local root = git_root(vim.uv.cwd())
if root then
	vim.env.ROOT = root
	vim.opt.path:prepend(root)
	vim.opt.path:append({ root .. "/src", root .. "/src/*" })
end
local roottop = git_root_top(vim.uv.cwd())
if roottop and roottop ~= root then
	vim.env.ROOTTOP = roottop
	vim.opt.path:append(roottop)
end

-- Look up a shell/Makefile-style `NAME=value` assignment in the buffer.
-- Matches line-anchored assignments with no spaces around `=` (zsh/bash/make
-- convention), strips surrounding single or double quotes. Returns nil if not
-- found. Scans from the bottom so later assignments win.
local function buffer_var(name, bufnr)
	local lines = vim.api.nvim_buf_get_lines(bufnr, 0, -1, false)
	for i = #lines, 1, -1 do
		local val = lines[i]:match("^%s*" .. name .. "=(.+)")
		if val then
			val = val:match("^%s*(.-)%s*$")
			return val:match("^['\"](.-)['\"]$") or val
		end
	end
end

-- gf with $(git root) / $VAR / buffer-local var expansion. Stock gf handles
-- $VAR/path via isfname + env, but breaks on $(...) (chars not in isfname, no
-- command substitution) and on vars that exist only as assignments in the
-- current file (e.g. `outdir=fpath` in a zsh script). We resolve:
--   * $(git root) and $(git rev-parse --show-toplevel) only — arbitrary shell
--     from a file would be a footgun.
--   * $VAR from env (vim.fn.expand), falling back to `VAR=value` in the buffer.
vim.keymap.set("n", "gf", function()
	local line = vim.api.nvim_get_current_line()
	local col = vim.api.nvim_win_get_cursor(0)[2] + 1
	local bufnr = vim.api.nvim_get_current_buf()
	local match
	for _, pat in ipairs({
		"%$%b()[%w/%._%-+~=@#%%]*",
		"%$[%w_]+[%w/%._%-+~=@#%%]*",
	}) do
		local start = 1
		while true do
			local s, e = line:find(pat, start)
			if not s then
				break
			end
			if col >= s and col <= e then
				match = line:sub(s, e)
				break
			end
			start = e + 1
		end
		if match then
			break
		end
	end
	if match then
		local token = match:gsub("%$%b()", function(m)
			local cmd = m:sub(3, -2):match("^%s*(.-)%s*$")
			if cmd == "git root" or cmd == "git rev-parse --show-toplevel" then
				return git_root(bufnr) or vim.env.ROOT or m
			end
			return m
		end)
		-- Env var first; only fall back to buffer assignment if unset.
		token = token:gsub("%$([%w_]+)", function(name)
			if vim.env[name] then
				return "$" .. name
			end
			return buffer_var(name, bufnr) or ("$" .. name)
		end)
		local path = vim.fn.expand(token)
		if vim.uv.fs_stat(path) then
			vim.cmd.edit(vim.fn.fnameescape(path))
			return
		end
	end
	vim.cmd("normal! gf")
end, { desc = "gf with $(git root), $VAR, buffer var expansion" })
