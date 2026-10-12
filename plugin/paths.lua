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

-- 'includeexpr' is local to the buffer and an ftplugin setting it replaces the
-- global value, so FileType moves a filetype's own expression to
-- b:ft_includeexpr, which paths.includeexpr tries first. The ftplugin FileType
-- autocmd is defined before plugin/ is sourced (:h startup), so this one runs
-- after it.
vim.o.includeexpr = paths.INCLUDEEXPR
vim.api.nvim_create_autocmd("FileType", {
	group = vim.api.nvim_create_augroup("paths_includeexpr", {}),
	callback = function(args)
		local bo, b = vim.bo[args.buf], vim.b[args.buf]
		if bo.includeexpr ~= "" and bo.includeexpr ~= paths.INCLUDEEXPR then
			b.ft_includeexpr = bo.includeexpr
			b.ft_includeexpr_ft = bo.filetype
		-- Our value with the filetype unchanged: FileType fired again (lz.n
		-- does) and the ftplugin's b:did_ftplugin guard skipped it, as there was
		-- no b:undo_ftplugin to clear the guard. The stored value still holds.
		elseif b.ft_includeexpr_ft ~= bo.filetype then
			b.ft_includeexpr = nil
			b.ft_includeexpr_ft = nil
		end
		bo.includeexpr = paths.INCLUDEEXPR
	end,
})

-- The gf family with $(git root) / $VAR / buffer-local var expansion and a
-- :lnum[:col] suffix to land on. The cursor location is tried before native gf
-- because native gf often finds something else in a $ token: on `git` in
-- `$(git root)/x` it opens /usr/bin/git, as 'path' holds $PATH (top of file).
-- The native key stays as the net for the rest: 'path', the count,
-- includeexpr and directories.

--- A mapping callback: open the location under the cursor with `open` and land
--- on its :lnum[:col], else run the native `key` with the count.
---@param key string the native key
---@param open fun(path: string) shows `path` in a window
---@return fun()
local function goto_location(key, open)
	return function()
		local path, lnum, col = paths.cursor_location()
		-- The buffer's directory first: native gf searches 'path', where $ROOT
		-- comes before ".". A count picks among 'path' matches, so it stays
		-- native.
		if not path and vim.v.count == 0 then
			local cfile = vim.fn.expand("<cfile>") --[[@as string]]
			path = paths.file_in_buf_dir(cfile, 0)
			if not path then
				local name = paths.transform_name(cfile, 0)
				path = name and paths.file_in_buf_dir(name, 0)
			end
		end
		if not path then
			local ok, err = pcall(vim.cmd.normal, { vim.v.count1 .. vim.keycode(key), bang = true })
			-- Shown as the native key shows it, not as a Lua traceback.
			if not ok then
				vim.api.nvim_echo({ { tostring(err):match("E%d+:.*") or tostring(err) } }, true, { err = true })
			end
			return
		end
		util.jumplist_add()
		open(path)
		if lnum then
			vim.fn.cursor(lnum, col or 1) -- clamps; nvim_win_set_cursor would throw
			vim.cmd.normal { "zv", bang = true } -- the target may sit inside a closed fold
		end
	end
end

--- An opener that shows a path in a new window made by the Ex command
--- `split_cmd`, and closes that window again if showing the path fails.
---@param split_cmd string e.g. "split" or "tab split"
---@return fun(path: string)
local function in_new_window(split_cmd)
	return function(path)
		vim.cmd(split_cmd)
		local ok, err = pcall(util.edit, path)
		if not ok then
			vim.cmd.close()
			error(err, 0)
		end
	end
end

local DESC = " with $(git root), $VAR, buffer var, :lnum"
local split_edit = in_new_window("split")
vim.keymap.set("n", "gf", goto_location("gf", util.edit), { desc = "gf" .. DESC })
vim.keymap.set("n", "<C-w>f", goto_location("<C-w>f", split_edit), { desc = "Split gf" .. DESC })
vim.keymap.set("n", "<C-w><C-f>", goto_location("<C-w><C-f>", split_edit), { desc = "Split gf" .. DESC })
vim.keymap.set("n", "<C-w>gf", goto_location("<C-w>gf", in_new_window("tab split")), { desc = "Tab gf" .. DESC })
