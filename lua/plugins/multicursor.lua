local hi = require "utils/highlights"
local util = require "utils/init"

return {
    {
        "mg979/vim-visual-multi",
        init = function ()
            vim.g.VM_mouse_mappings = true
            vim.g.VM_leader = "<leader>m"
            vim.g.VM_theme = "neon"
        end,
        config = function ()
            -- use alt instead of ctrl since ctrl arrows moves mac os windows.
            -- Other keybindings (e.g. [] and y) cannot be set by multicursor since they are set elsewhere, hence startup warnings.
            -- Ignore for now since we can live without them.
            local set_keymap_desc = function(...) pcall(require"mini.clue".set_keymap_desc, ...) end
            set_keymap_desc('n', "<C-Down>",  "Multi cursor down")
            set_keymap_desc('n', "<C-Up>",  "Multi cursor up")
        end,
    },
    {
        "smoka7/multicursors.nvim",
        enabled = false, -- doesn't insert text correctly, e.g. try ""<S-space>hello
        event = "VeryLazy",
        dependencies = {
            'nvim-treesitter/nvim-treesitter',
            'smoka7/hydra.nvim',
        },
        opts = {
            generate_hints = {
                normal = true,
                insert = true,
                extend = true,
            },
        },
        cmd = { 'MCstart', 'MCvisual', 'MCclear', 'MCpattern', 'MCvisualPattern', 'MCunderCursor' },
        keys = {
            {
                mode = { 'v', 'n' },
                '<Leader>m',
                '<Cmd>MCstart<CR>',
                desc = 'Create a selection for selected text or word under the cursor',
            },
            {
                mode = { 'v', 'n' },
                '<Leader>M',
                '<Cmd>MCunderCursor<CR>',
                desc = 'Select the char under the cursor and start listening for the actions.',
            },
        },
    },
    {
        "jake-stewart/multicursor.nvim",
        enabled = false, -- not copying well enough what is inserted at main cursor
        -- copied default config, then modified
        config = function()
            local mc = require("multicursor-nvim")

            mc.setup()

            local set = vim.keymap.set

            -- Add or skip cursor above/below the main cursor.
            set({"n", "v"}, "<S-up>", util.fncount(mc.lineAddCursor, -1))
            set({"n", "v"}, "<S-down>", util.fncount(mc.lineAddCursor,1))
            set({"n", "v"}, "<S-A-up>", util.fncount(mc.lineSkipCursor,-1))
            set({"n", "v"}, "<S-A-down>", util.fncount(mc.lineSkipCursor,1))

            -- Add or skip adding a new cursor bymatching word/selection
            set({"n", "v"}, "<C-n>", util.fncount(mc.matchAddCursor,1))
            set({"n", "v"}, "<C-S-n>", util.fncount(mc.matchSkipCursor,1))
            set({"n", "v"}, "<C-p>", util.fncount(mc.matchAddCursor,-1))
            set({"n", "v"}, "<C-S-p>", util.fncount(mc.matchSkipCursor,-1))

            -- You can also add cursors with any motion you prefer:
            -- set("n", "<right>", function()
            --     mc.addCursor("w")
            -- end)
            -- set("n", "<leader><right>", function()
            --     mc.skipCursor("w")
            -- end)

            -- Rotate the main cursor.
            set({"n", "v"}, "<S-left>",  util.fncount(mc.prevCursor), {desc="Multicursor rotate previous"})
            set({"n", "v"}, "<S-right>", util.fncount(mc.nextCursor), {desc="Multicursor rotate next"})

            -- Delete the main cursor.
            set({"n", "v"}, "<leader>mx", mc.deleteCursor, {desc="Delete main cursor"})

            -- Add and remove cursors with control + left click.
            set("n", "<c-leftmouse>", mc.handleMouse)

            set({"n", "v"}, "<leader>mm", function()
                if mc.cursorsEnabled() then
                    -- Stop other cursors from moving.
                    -- This allows you to reposition the main cursor.
                    mc.disableCursors()
                else
                    mc.addCursor()
                end
            end, {desc="Lock multi/mark new if locked"})

            -- clone every cursor and disable the originals
            set({"n", "v"}, "<leader>mc", mc.duplicateCursors, {desc="Clone+lock cursors"})

            set("n", "<esc>", function()
                if not mc.cursorsEnabled() then
                    mc.enableCursors()
                elseif mc.hasCursors() then
                    mc.clearCursors()
                else
                    -- Default <esc> handler placed here since we override our custom esc remap that stops /-search hl
                    vim.cmd.nohlsearch()
                end
            end)

            -- Align cursor columns.
            set("v", "<leader>ma", mc.alignCursors, {desc="Align"})

            -- Split visual selections by regex.
            set("v", "<leader>ms", mc.splitCursors, {desc="Split by regex"})

            -- Append/insert for each line of visual selections.
            set("v", "I", mc.insertVisual, {desc="Insert (multi)"})
            set("v", "A", mc.appendVisual, {desc="Append (multi)"})

            -- match new cursors within visual selections by regex.
            set("v", "<leader>mr", mc.matchCursors, {desc="match new cursors within visual selections by regex"})

            -- Rotate visual selection contents.
            set("v", "]m", function() mc.transposeCursors( util.count()) end, {desc="Rotate contents forward"})
            set("v", "[m", function() mc.transposeCursors(-util.count()) end, {desc="Rotate contents backward"})

            local grp = vim.api.nvim_create_augroup("colorscheme", {clear=true})
            vim.api.nvim_create_autocmd("ColorScheme", {
                pattern = "*",
                group = grp,
                callback = function ()
                    -- A dim version of the cursor
                    hi.set("MultiCursorCursor", {fg=hi.getfg("Cursor"), bg="gray"})
                    hi.link("MultiCursorVisual", "Visual")
                    hi.link("MultiCursorSign", "SignColumn")
                    -- disable aka locked cursor. Changed default link to visual (which feels misleading) to an even dimmer cursor.
                    hi.set("MultiCursorDisabledCursor", {fg=hi.getfg("Cursor"), bg=hi.getfg("NonText")})
                    hi.link("MultiCursorDisabledVisual", "Visual")
                    hi.link("MultiCursorDisabledSign", "SignColumn")
                end
            })
        end
    },
}
