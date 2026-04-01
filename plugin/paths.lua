-- append $PATH to vim option &path (&path uses comma delim in text form and $PATH uses colon)
vim.opt_local.path:append(vim.split(vim.env.PATH, ':'))

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

-- PyMOL stubs for LSP are in basedpyright extraPaths (lsp/basedpyright.lua),
-- not PYTHONPATH — global PYTHONPATH breaks the real pymol binary.

-- Defer git root lookup — avoid blocking startup with subprocess
vim.schedule(function()
    -- prepend git ROOT
    -- Note the '$'. It makes ROOT and env var which means it is available for e.g.
    -- path completion so that $ROOT/ will have completion, e.g. in julia or zsh.
    vim.system({'git', 'root'}, {text=true}, function(obj)
        vim.schedule(function()
            local root = obj.stdout and obj.stdout:gsub("\n$", "") or ""
            if root == "" then return end
            vim.env.ROOT = root
            vim.opt.path:prepend(root)
            vim.opt.path:append {
                root .. "/src",
                root .. "/src/*",
            }
            -- Note that for a nested git repo, we will focus on the immediate root (found with `git root`),
            -- but also add a reference to the top level root:
            local roottop = vim.fn.finddir('.git/..', vim.fn.expand('%:p:h') .. ';')
            if roottop ~= "" then
                vim.env.ROOTTOP = roottop
                -- will automatically not append if ROOTTOP == ROOT (no duplicate entries)
                vim.opt.path:append(roottop)
            end
        end)
    end)
end)
