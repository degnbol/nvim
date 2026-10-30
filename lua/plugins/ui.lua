return {
    -- :Capture hi to call :hi where you can search etc.
    { "capture.vim", cmd = "Capture" },
    {
        "kittyREPL.nvim",
        load = function() end,
        after = function()
            require("kittyREPL").setup {
                keymap = {
                    focus = "<C-CR>",
                    setlast = "<leader>rr",
                    new = "<leader>rs",
                    run = "<CR>",
                    paste = "<S-CR>",
                    help = "<leader>K",
                    runLine = "<CR><CR>",
                    runLineFor = "<leader>ro",
                    runLineForI = "<leader>ri",
                    pasteLine = "<S-CR><S-CR>",
                    runVisual = "<CR>",
                    pasteVisual = "<S-CR>",
                    q = "<leader>rq",
                    cr = "<leader>r<S-CR>",
                    ctrld = "<leader>rd",
                    ctrlc = "<leader>rc",
                    interrupt = "<leader>rk",
                    scrollStart = "[r",
                    scrollUp = "[r",
                    scrollDown = "]r",
                    scrollOutputAbove = "[R",
                    scrollOutputBelow = "]R",
                    progress = "<leader>rp",
                    editPaste = "<leader>re",
                },
                exclude = { tex = true, tsv = true, markdown = true },
                progress = true,
                editpaste = true,
                closepager = true,
            }
        end,
    },
    -- File explorer as a buffer with manipulation abilities.
    -- Differes from mini.files by only having a single view taking up the whole screen and has different default keymaps and preview behaviour.
    {
        "oil.nvim",
        keys = {
            { "<leader>-", function() require("oil").open() end, desc = "Oil" },
        },
        after = function()
            require("oil").setup {
                -- no prompt on rename, including folder change.
                -- The prompt is a nice list of the performed changes.
                -- skip_confirm_for_simple_edits = true,
                keymaps = {
                    -- I regularly /-search for filename and don't need to hl that inside the buffer after opening
                    ["<CR>"] = function()
                        vim.cmd.nohlsearch()
                        require "oil".select()
                    end
                }
            }
        end,
    },
    {
        "fyler.nvim",
	lazy = true,
	cmd = "Fyler",
        after = function()
            require("fyler").setup {
                icon_provider = "nvim-web-devicons",
            }
        end,
    },
    -- https://github.com/jake-stewart/auto-cmdheight.nvim
    -- Auto-sizes the cmdline height to fit content, avoiding many "Hit Enter" prompts.
    -- Redudant when using the new ui2, and in fact interferes and breaks the UI when both are active.
    -- With both active, :Inspect and other multi-line messages show blank.
    {
        "auto-cmdheight.nvim",
        enabled = not package.loaded["vim._core.ui2"],
        lazy = false,
        after = function()
            require("auto-cmdheight").setup {
                -- max cmdheight before displaying hit enter prompt.
                max_lines = 5,

                -- number of seconds until the cmdheight can restore.
                duration = 2,

                -- whether key press is required to restore cmdheight.
                remove_on_key = true,

                -- always clear the cmdline after duration and key press.
                -- by default it will only happen when cmdheight changed.
                clear_always = false,
            }
        end,
    },
    {
        'incline.nvim',
        before = function()
            -- This plugin only makes sense with global statusline,
            -- or no statusline at all.
            -- The latter is possible with a trick.
            -- By setting cmdheight=0 and (optionally) clearing the
            -- 'statusline' setting we can use a single line for cmdline and
            -- (secret) global statusline.
            -- The problem with this approach is how buggy cmdheight=0 still is.
            -- A different solution is laststatus=0 where the statusline is redesigned as a simple border.
            -- The first agentic session switches to laststatus=3 for the budget (see agents.lua).
            vim.o.laststatus = 0
            -- Draw the border with the statusline FILL char, not fixed-width
            -- content. The fill auto-sizes to each window, so it never truncates
            -- (no stray '<' in narrow splits) and needs no WinResized recompute.
            vim.opt.fillchars:append("stl:─,stlnc:─")
            vim.o.statusline = "%="
        end,
        after = function()
            local helpers = require 'incline.helpers'
            local devicons = require 'nvim-web-devicons'
            local hi = require 'utils/highlights'
            -- Refresh incline when agentic.nvim updates header state
            vim.api.nvim_create_autocmd("User", {
                pattern = "AgenticHeadersChanged",
                callback = function() require('incline').refresh() end,
            })
            require('incline').setup {
                ignore = {
                    -- Also display for help etc.
                    buftypes = { 'quickfix' },
                    wintypes = {},
                    unlisted_buffers = false,
                    filetypes = { 'AgenticFiles', 'AgenticTodos', 'AgenticActivity', 'DiffviewFiles' },
                },
                window = {
                    padding = 0,
                    margin = { horizontal = 0, vertical = 0, },
                    placement = {
                        -- Keep at top since it fits well with TS-context that also provides breadcrumps.
                        -- vertical = "bottom",
                    },
                    overlap = {
                        borders = false,
                        -- statusline = true,
                        -- tabline = true,
                        -- winbar = true,
                    },
                    -- Above fidget if margin.vertical =1 and placement.vertical = "bottom" but not if margin.vertical = 0
                    -- zindex = 100,
                },
                render = function(props)
                    local ft = vim.bo[props.buf].filetype
                    -- The main input has its own place in the widget; a
                    -- subagent's message input does not.
                    if ft == 'AgenticInput' and vim.b[props.buf].agentic_window ~= 'message' then
                        return nil
                    end
                    if ft == 'AgenticChat' or ft == 'AgenticInput' then
                        -- The header parts are owned by agentic.nvim;
                        -- surface what it publishes rather than rebuild it here.
                        -- Permission and idle states wait for the user, so they get the warn colour.
                        local deco = require('agentic.ui.window_decoration')
                        local part = deco.get_header(props.buf)
                        local state = deco.header_state(part)
                        local color = props.focused and hi.fg("StatusLine") or "gray"
                        local state_color = (state == 'permission' or state == 'idle') and hi.fg("DiagnosticWarn") or color
                        return {
                            part.icon and { part.icon .. ' ', guifg = color } or '',
                            { part.title, guifg = color },
                            state and { ' ' .. require('agentic.glyphs').HEADER_STATE[state], guifg = state_color } or '',
                            part.context and { ' · ' .. part.context, guifg = color } or '',
                        }
                    end
                    -- In diffview tabs: show rev info instead of filename
                    local dv_ok, dv_lib = pcall(require, 'diffview.lib')
                    if dv_ok then
                        local view = dv_lib.get_current_view()
                        if view and view.cur_layout then
                            for _, key in ipairs({ 'a', 'b', 'c', 'd' }) do
                                local w = view.cur_layout[key]
                                if w and w.id == props.win and w.file and w.file.winbar then
                                    local label = w.file.winbar:gsub("^ ", ""):gsub("^WORKING TREE %- ", "")
                                    local filename = vim.fn.fnamemodify(vim.api.nvim_buf_get_name(props.buf), ':t')
                                    local ft_icon, ft_color = devicons.get_icon_color(filename, nil, { default = true })
                                    local color = props.focused and hi.fg("StatusLine") or "gray"
                                    return {
                                        ft_icon and { ft_icon, guifg = ft_color } or '',
                                        ' ',
                                        { label, guifg = color },
                                    }
                                end
                            end
                            -- Diffview window but no winbar (e.g. file panel) — hide
                            return ''
                        end
                    end
                    local filename = vim.fn.fnamemodify(vim.api.nvim_buf_get_name(props.buf), ':t')
                    if filename == '' then
                        filename = '[No Name]'
                    end
                    local ft_icon, ft_color = devicons.get_icon_color(filename, nil, { default = true })
                    local modified = vim.bo[props.buf].modified
                    return {
                        ft_icon and { ft_icon, guifg = ft_color } or '',
                        ' ',
                        { filename, gui = modified and 'bold' or nil, guifg = props.focused and hi.fg("StatusLine") or "gray" },
                    }
                end,
            }
        end,
        event = 'DeferredUIEnter',
    },
}
