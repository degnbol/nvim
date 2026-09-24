-- append $PATH to vim option &path (&path uses comma delim in text form and $PATH uses colon)
vim.opt.path:append(vim.split(vim.env.PATH, ":"))

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

-- PyMOL stubs for LSP are in basedpyright extraPaths (after/lsp/basedpyright.lua),
-- not PYTHONPATH — global PYTHONPATH breaks the real pymol binary.

local paths = require("utils.paths")
local util = require("utils/init")
local git_root = paths.git_root

-- Like git_root but filters to .git *directories*, which means it skips
-- submodule gitlinks and returns the outermost enclosing repo.
local function git_root_top(source)
	local path = type(source) == "number" and vim.api.nvim_buf_get_name(source) or source
	path = path ~= "" and vim.fs.dirname(path) or vim.uv.cwd()
	local matches = vim.fs.find(function(name, p)
		if name ~= ".git" then
			return false
		end
		return paths.is_dir(vim.fs.joinpath(p, name))
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

-- gf with $(git root) / $VAR / buffer-local var expansion and a :lnum[:col]
-- suffix to land on (resolver in utils.paths). `normal! gf` stays as the final
-- net for what the resolver leaves: includeexpr and the count'th match in
-- 'path' (:h gf).
vim.keymap.set("n", "gf", function()
	local path, lnum, col = paths.resolve_location_under_cursor(0)
	if not path then
		vim.cmd("normal! " .. vim.v.count1 .. "gf")
		return
	end
	util.jumplist_add()
	util.edit(path)
	if lnum then
		vim.fn.cursor(lnum, col or 1) -- clamps; nvim_win_set_cursor would throw
		vim.cmd("normal! zv") -- the target may sit inside a closed fold
	end
end, { desc = "gf with $(git root), $VAR, buffer var, :lnum" })
