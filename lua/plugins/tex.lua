local util = require "utils/init"
local latexmk = require "tex.latexmk"
local g = vim.g

return {
    -- more aggressive conceal
    {
        "KeitaNakamura/tex-conceal.vim",
        ft = "tex",
        config = function()
            -- https://github.com/gillescastel/latex-snippets
            g.tex_conceal = 'abdmg'
        end,
    },
    {
        "lervag/vimtex",
        -- :checkhealth suggests pstree for inverse search.
        -- 2>/dev/null in case we don't have brew
        build = "brew install pstree 2> /dev/null",
        init = function()
            -- instead of <localleader>l
            -- We want <leader>l for LSP and double leader (i.e. localleader) for filetype specific mappings.
            g.vimtex_mappings_prefix = "<localleader>"

            -- g.vimtex_view_method = 'sioyek' -- not working
            g.vimtex_view_method = 'skim'
            -- inspiration from https://dr563105.github.io/blog/skim-vimtex-setup/
            -- forward search after every successful compilation.
            -- Changed back to default, can get annoying with jumping around with slow compilation.
            -- Use manual VimtexView (<leader>cv) instead to do forward search.
            g.vimtex_view_skim_sync = false
            -- change focus to skim after command `:VimtexView` is given
            g.vimtex_view_skim_activate = true
            g.vimtex_view_skim_reading_bar = true

            local ignore = {
                "Underfull",
                "Overfull",
                "biblatex: Duplicate entry key",          -- same paper in two zotero folders
                "Text page [0-9]* contains only floats.", -- on purpose in some cases for large full page figure
                -- It's ok to have glossaries package in preamble without using it (e.g. yet).
                -- \\ in literal tested to work
                [[No \\printglossary or \\printglossaries found]],
                [[Warning: Command \\underbar  has changed.]],
                [[Warning: Command \\overline  has changed.]],
                [[Warning: Empty bibliography]],
            }
            g.vimtex_log_ignore = ignore
            g.vimtex_quickfix_ignore_filters = ignore
            g.vimtex_quickfix_open_on_warning = false

            -- :h vimtex-af-enhanced-matchparen
            -- wrong matching for some {} with this on:
            -- vim.g.matchup_override_vimtex = true

            -- defaults at :h vimtex_compiler_latexmk
            g.vimtex_compiler_latexmk = {
                aux_dir = "aux",
                -- out_dir = "out", -- compile selected only works if using default
                -- callback = true,
                options = {
                    -- defaults:
                    "-verbose",
                    "-file-line-error",
                    "-synctex=1",
                    "-interaction=nonstopmode",
                    -- not defaults:
                    "--shell-escape", -- to let package minted run pygmentize
                },
            }
            -- set default latex engine to the modern lualatex over pdflatex
            g.vimtex_compiler_latexmk_engines = {
                _ = "-lualatex",
            }

            -- formatter when calling gq, but note that autoformatting calls the
            -- builtin vim formatter which you can call with gw.
            -- In order to not autoformat in math mode I made an autocmd to
            -- disable/enable auto formatoption based on vimtex's detection in_mathzone.
            -- See ftplugin/tex.vim
            g.vimtex_format_enabled = true
            -- trying to avoid indent after \item\n. Why? when the "\item" is
            -- concealed the text doesn't align nicely which must have been the
            -- original intent. Also, my snippet expanding a '-' for a new item becomes more complicated.
            -- vimtex solution doesn't work: https://github.com/lervag/vimtex/issues/2599
            -- use builtin indent function GetTexIndent since the solution for that does work:
            -- https://github.com/vim/vim/blob/eb3dc87f01391bb075d97aef3d00f91b4e08a25c/runtime/indent/tex.vim#L70-L122
            -- Update: I don't use conceal anymore, and we want \if \else etc. to have indent.
            g.vimtex_indent_enabled = true
            -- builtin vim setting for tex files, not a vimtex setting
            g.tex_indent_items = 0

            g.vimtex_toc_config = {
                indent_levels = 1, -- has no effect but would be nice
                mode = 4,          -- the only mode that works
                -- layers = {'content', 'todo', 'include'}, -- don't list labels
            }

            -- if in a subfile, by default we compile only that.
            -- <localleader>m to toggle compiling main instead.
            g.vimtex_subfile_start_local = true -- doesn't seem to work.
            -- Also toggling is inconsistent so just use it and look at the message.

            -- Since we use function as bold func def and function.call as unbold, we relink:
            vim.defer_fn(function()
                vim.api.nvim_set_hl(0, "texCmd", { link = "@function.call", force = true })
            end, 1000)

            local grp = vim.api.nvim_create_augroup("vimtex", { clear = true })

            vim.api.nvim_create_autocmd("User", {
                pattern = "*VimtexEventCompileStarted*",
                group = grp,
                callback = function()
                    local VimtexCompiling = true
                end
            })

            vim.api.nvim_create_autocmd("User", {
                pattern = "VimtexEventCompileFailed",
                group = grp,
                callback = function()
                    -- Make function arg if needed
                    local opts = {
                        main = "main",
                        aux  = "aux",
                    }
                    local auxs = vim.fs.find(opts.aux, { upward = true, limit = 5 })
                    if #auxs == 0 then return end
                    local search_pattern = "ERROR - " .. opts.aux .. "/" .. opts.main .. ".bcf is malformed"
                    vim.system({ "grep", search_pattern, "main.blg" }, { cwd = auxs[1] }, function(obj)
                        if obj.code == 0 then -- search pattern found, i.e. main.bcf is malformed
                            print(opts.main .. ".bcf malformed. Cleaning...")
                            vim.schedule(function()
                                vim.fn["vimtex#compiler#clean"](0)
                                vim.defer_fn(function()
                                    latexmk.is_running(function(is_running)
                                        -- Since this autocmd is triggered by failed compile it means we just tried to compile.
                                        -- Then we should either be in continuous mode with a latexmk process running (is_running == true),
                                        -- or it was a single shot compile that failed, hence we redo single shot compile here.
                                        if not is_running then
                                            vim.fn["vimtex#compiler#compile_ss"]()
                                        end
                                    end)
                                end, 500)
                            end)
                        end
                    end)
                end
            })

            -- Set tex main file by looking for a file upwards named "main.tex"
            vim.api.nvim_create_autocmd("BufReadPre", {
                pattern = "*.tex",
                group = vim.api.nvim_create_augroup("TexMain", { clear = true }),
                callback = function()
                    local mainfile = vim.fs.find("main.tex", { type = "file", upward = true })[1]
                    if mainfile ~= nil then vim.b.vimtex_main = mainfile end
                end
            })

            -- TODO: Do we still need this as a bdelete alternative mapping?
            ---Delete buffer. Repeat for unnamed empty buffers.
            ---@param opts table with bool key force (passed to vim.api.nvim_buf_delete)
            ---@param lastbufnr integer? for recursion
            local function bufdel(opts, lastbufnr)
                opts = opts or {}
                local bufnr = vim.api.nvim_get_current_buf()
                -- make sure to not retry if a previous call failed
                if bufnr == lastbufnr then return end
                vim.api.nvim_buf_delete(0, opts)
                -- repeat if next buffer is empty (stop annoying behaviour of vimtex)
                if not util.is_named() and util.is_empty() then
                    vim.schedule(function() bufdel(opts, bufnr) end)
                end
            end
        end,
    },
}
