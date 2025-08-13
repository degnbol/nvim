-- append $PATH to vim option &path (&path uses comma delim in text form and $PATH uses colon)
vim.opt_local.path:append(vim.split(vim.env.PATH, ':'))

-- prepend git ROOT
-- Note the '$'. It makes ROOT and env var which means it is available for e.g. 
-- path completion so that $ROOT/ will have completion, e.g. in julia or zsh.
-- set env variable ROOT to git root:
vim.env.ROOT = vim.system({'git', 'root'}, {text=true}):wait().stdout:gsub("\n$", "")
-- Note that for a nested git repo, we will focus on the immediate root (found with `git root`), 
-- but also add a reference to the top level root:
vim.env.ROOTTOP = vim.fn.finddir('.git/..', vim.fn.expand('%:p:h') .. ';')

if vim.env.ROOT ~= nil then
    vim.opt_local.path:prepend(vim.env.ROOT)
    vim.opt_local.path:append {
        vim.env.ROOT .. "/src",
        vim.env.ROOT .. "/src/*",
    }
end
if vim.env.ROOTTOP ~= nil then
    -- will automatically not append if ROOTTOP == ROOT (no duplicate entries)
    vim.opt_local.path:append(vim.env.ROOTTOP)
end

-- edit-in-kitty on remotes doesn't copy the env variables that are normally 
-- present locally so we have to set the necessary ones. These are for julia 
-- and python lsp to find the right binaries. These were found by copying all 
-- of `env` and checking which fixed the problem by adding them to 
-- `edit-in-kitty --env ...` which can be given multiple times like the kitty @ 
-- launch command.
-- Hardcoding miniforge base env path. It would be cooler to use activated env 
-- somehow.
-- ~/.local/bin is location of jedi-language-server
vim.env.PATH = vim.env.HOME .. '/.local/bin/:/opt/homebrew/Caskroom/miniforge/base/bin:/opt/homebrew/bin:' .. vim.env.PATH
vim.env.CONDA_PREFIX = vim.env.HOME .. '/miniconda3'

-- add lsp/completion for pymol.
-- In .py script, write `from pymol import cmd`
local rtp = vim.opt.runtimepath:get()[1]
if vim.env.PYTHONPATH == nil then
    vim.env.PYTHONPATH = rtp .. "/lsp/pymol-open-source/modules"
else
    vim.env.PYTHONPATH = vim.env.PYTHONPATH .. ':' .. rtp .. "/lsp/pymol-open-source/modules"
end

