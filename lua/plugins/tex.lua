local util = require "utils/init"
local hi = require "utils/highlights"
local latexmk = require "tex.latexmk"
local g = vim.g

return {
    -- more aggressive conceal
    {
        "tex-conceal.vim",
        ft = "tex",
        after = function()
            -- https://github.com/gillescastel/latex-snippets
            g.tex_conceal = 'abdmg'
        end,
    },
    {
        "vimtex",
        -- :checkhealth suggests pstree for inverse search.
        -- 2>/dev/null in case we don't have brew
        before = function()
            -- instead of <localleader>l
            -- We want <leader>l for LSP and double leader (i.e. localleader) for filetype specific mappings.
            g.vimtex_mappings_prefix = "<localleader>"

            if util.is_mac() then
                g.vimtex_view_method = 'skim'
                -- inspiration from https://dr563105.github.io/blog/skim-vimtex-setup/
                -- forward search after every successful compilation.
                -- Changed back to default, can get annoying with jumping around with slow compilation.
                -- Use manual VimtexView (<leader>cv) instead to do forward search.
                g.vimtex_view_skim_sync = false
                -- change focus to skim after command `:VimtexView` is given
                g.vimtex_view_skim_activate = true
                g.vimtex_view_skim_reading_bar = true
            else
                -- :h vimtex-view-zathura
                -- There's also 'zathura', which requires xdotool to check for other zathura viewers running with
                -- `xdotool search --class Zathura`
                -- This only works on x11, and we are currently on wayland (hyprland).
                -- https://github.com/lervag/vimtex/issues/2046
                -- g.vimtex_view_method = 'zathura'
                g.vimtex_view_method = 'zathura_simple'
            end

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

            require("tex.textcolor").setup()

            -- fix missing or inconsistent hl links on every ColorScheme + now
            hi.onColorScheme(function()
                -- Pick a reduced colour for removing emphasis on things like \cite{...} where the body's color and underline gives it emphasis by itself.
                -- We want to differentiate from comment and nontext, and nontext is bold so the fg with italic should be enough differentiation, plus we would write comments more that using nontext.
                local gray = hi.fg("NonText")
                -- Since we use function as bold func def and function.call as unbold, we relink:
                hi.link("texCmd", "@function.call")
                hi.link("texCmdEnv", "@keyword.function") -- italic instead of bold for begin end
                hi.link("texCmdRef", "@function.builtin") -- italic
                -- italic \section{...}, bold etc. Gray a bit since the "..." shows aesthetic
                hi.set("texCmdPart", { fg = gray, italic = true })
                hi.set("texCmdStyleBold", { fg = gray, italic = true })
                hi.set("texCmdStyleItal", { fg = gray, italic = true })
                hi.set("texTypeStyle", { fg = gray, italic = true })       -- e.g. \underline
                hi.set("texItalStyle", { italic = true })                  -- contents of \emph{...}
                hi.set("texCmdRefConcealed", { fg = gray, italic = true }) -- italic \cite
                hi.set("texCmdRef", { fg = gray, italic = true })
                hi.set("texCmdCRef", { fg = gray, italic = true })
                hi.set("texCmdAcro", { fg = gray })           -- custom cmd defined in after/syntax/tex.vim
                hi.link("texCmdPackage", "@function.builtin") -- italic \package
                hi.link("texCmdInput", "@function.builtin")   -- italic \inputgraphics
                hi.link("texCmdTitle", "@function.builtin")   -- italic \title
                hi.link("texCmdAuthor", "@function.builtin")  -- italic \author
                hi.link("texCmdLet", "@function.builtin")     -- italic \let
                hi.link("texStatement", "@function.builtin")  -- only seen for \mathrm so far
                hi.mod("texMatcher", { underline = true })    -- matched parenthesis, \underline body, etc.
                hi.link("texEnvArgName", "@method")           -- bold and shine instead of nothing
                hi.link("texCmdBeamer", "@function")
                hi.link("texOpt", "@parameter")
                hi.link("texBeamerOpt", "@parameter")
                hi.link("texOptEqual", "@operator")
                hi.link("texArg", "@parameter")
                hi.link("texFileArg", "@string")
                hi.link("texFilesArg", "@string")
                hi.link("texFileOpt", "@parameter")
                hi.link("TexBeamerDelim", "Delimiter")
                hi.link("superscript", "Type")                                               -- like \huge, \normalsize etc
                hi.link("subscript", "Type")                                                 -- like \huge, \normalsize etc
                hi.set("texRefConcealedArg", { fg = hi.fg("TexFileArg"), underline = true }) -- body of \cite{...}
                hi.link("texTitleArg", "Title")
                hi.link("texPartArgTitle", "Title")
                hi.link("texRefArg", "@tag")                                    -- body of \label
                hi.link("texSpecialChar", "@comment")                           -- unbreakable space ~, and \&
                hi.link("texMathZone", "@number")                               -- Most of tex math zone that isn't captured by anything else (such as math functions) is numbers and we don't use numbers much elsewhere.
                hi.set("texMathCmdText", { fg = gray, italic = true })          -- italic \text in math mode
                hi.set("texMathSymbol", { fg = hi.fg("@type"), italic = true }) -- type is similar colour to number
                hi.set("texMathSymbol", { fg = hi.fg("@type"), italic = true }) --
                hi.link("texSICmd", "@number")                                  -- not bold SI. Color like math mode
                hi.set("texLigature", { bold = true })                          -- bold instead of strong color to only give subtle focus to ``'', --, and the ' in don't
                hi.link("texCmdLigature", "@function.call")
                hi.mod("texCmdLigature", { italic = true })
                hi.link("texTabularChar", "Operator") -- & and \\ in tables. Could also use Delimiter but this makes them bold.
                hi.mod("texCmdClass", { italic = true, bold = true })
                hi.link("texOptSep", "Delimiter")
                hi.mod("texCmdDef", { bold = true, italic = true })    -- an actual function definition. \def. TeX primitive.
                hi.mod("texCmdNewcmd", { bold = true, italic = true }) -- an actual function definition. \newcommand. LaTeX wrapper on def.
                hi.link("texNewcmdArgName", "@parameter")
            end)

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
