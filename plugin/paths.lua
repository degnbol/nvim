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

-- gf for $(git root)/path. Stock gf handles $VAR/path via isfname + env
-- expansion, but $(...) has chars not in isfname and neovim won't run command
-- substitution. Only `git root` and `git rev-parse --show-toplevel` are
-- resolved — arbitrary shell from a file would be a footgun. Anchored at the
-- buffer's own repo (not nvim cwd) so paths in a file mean that file's repo.
vim.keymap.set("n", "gf", function()
	local line = vim.api.nvim_get_current_line()
	local col = vim.api.nvim_win_get_cursor(0)[2] + 1
	local bufnr = vim.api.nvim_get_current_buf()
	local start = 1
	while true do
		local s, e = line:find("%$%b()[%w/%._%-+~=@#%%]*", start)
		if not s then
			break
		end
		if col >= s and col <= e then
			local token = line:sub(s, e):gsub("%$%b()", function(m)
				local cmd = m:sub(3, -2):match("^%s*(.-)%s*$")
				if cmd == "git root" or cmd == "git rev-parse --show-toplevel" then
					return git_root(bufnr) or vim.env.ROOT or m
				end
				return m
			end)
			local path = vim.fn.expand(token)
			if vim.uv.fs_stat(path) then
				vim.cmd.edit(vim.fn.fnameescape(path))
				return
			end
			break
		end
		start = e + 1
	end
	vim.cmd("normal! gf")
end, { desc = "gf with $(git root) expansion" })
