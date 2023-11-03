#!/usr/bin/env lua

-- Auto disable hlsearch after use. The following works but is overkill to check every key press, just add :noh to nmap esc.
-- vim.on_key(function(char)
--   if vim.fn.mode() == "n" then
--     if vim.tbl_contains({ "<Esc>", "<Left>", "<Right>", "<Up>", "<Down>" }, vim.fn.keytrans(char)) then
--         vim.cmd.nohlsearch()
--     end
--   end
-- end, vim.api.nvim_create_namespace "auto_hlsearch")
vim.keymap.set('n', '<Esc>', ":noh<CR><Esc>", { silent=true, remap=false, desc="hello" })
-- disable when entering visual
local grp = vim.api.nvim_create_augroup("auto_hlsearch", {clear=true})
vim.api.nvim_create_autocmd("ModeChanged", {
    -- Mode changed is indicated as <before>:<after> with the codes:
    -- i=insert, c=cmdline, n=normal, v=visual, V=line visual, \x16=block visual (guessing it means <C-V>), no=normal operator pending.
    -- we need to avoid responding to changes between normal and cmdline mode since that change is triggered frequently by plugins etc.
    -- these patterns should capture the same:
    -- pattern = "*:*[oivV\x16]*",
    pattern = "*:*[^nc]",
    group = grp,
    callback = function ()
        vim.schedule(vim.cmd.nohlsearch)
    end
})


return {
    -- let's search result box show number of matches when there's >99 matches
    {
        "google/vim-searchindex",
        -- keys={"/", "g/"} -- keys doesn't seem to work
    },
    -- show a scrollbar, mostly in order to show search results far away in file
    -- requires hlslens to show search results in scrollbar
    {"petertriho/nvim-scrollbar", cmd={"ScrollbarToggle", "ScrollbarShow"},
    dependencies={"kevinhwang91/nvim-hlslens"},
    init = function ()
        vim.keymap.set("n", "yoS", "<Cmd>ScrollbarToggle<CR>", { desc="Scrollbar" })
    end,
    opts = {
        show = false, -- enable with :ScrollbarToggle etc.
        marks = { Search = { color = "yellow" } },
        handlers = {
            diagnostics = false,
            search = true
        }
    }},
    -- show counter for how many n or N's a search result is away from cursor
    -- also lots of other search highlight customizations possible
    {
        "kevinhwang91/nvim-hlslens",
        cmd = {"HlSearchLensEnable", "HlSearchLensToggle"},
        init = function ()
            vim.keymap.set("n", "yoH", function () require'hlslens'.toggle() end, { desc="HlSearchLens", silent=true })
            -- integrate with haya14busa/vim-asterisk
            vim.api.nvim_set_keymap('n', '*', [[<Plug>(asterisk-z*)<Cmd>lua require('hlslens').start()<CR>]],   {desc="Search word under cursor"})
            vim.api.nvim_set_keymap('n', '#', [[<Plug>(asterisk-z#)<Cmd>lua require('hlslens').start()<CR>]],   {desc="Search word under cursor"})
            vim.api.nvim_set_keymap('n', 'g*', [[<Plug>(asterisk-gz*)<Cmd>lua require('hlslens').start()<CR>]], {desc="Search word under cursor"})
            vim.api.nvim_set_keymap('n', 'g#', [[<Plug>(asterisk-gz#)<Cmd>lua require('hlslens').start()<CR>]], {desc="Search word under cursor"})
            vim.api.nvim_set_keymap('x', '*', [[<Plug>(asterisk-z*)<Cmd>lua require('hlslens').start()<CR>]],   {desc="Search word under cursor"})
            vim.api.nvim_set_keymap('x', '#', [[<Plug>(asterisk-z#)<Cmd>lua require('hlslens').start()<CR>]],   {desc="Search word under cursor"})
            vim.api.nvim_set_keymap('x', 'g*', [[<Plug>(asterisk-gz*)<Cmd>lua require('hlslens').start()<CR>]], {desc="Search word under cursor"})
            vim.api.nvim_set_keymap('x', 'g#', [[<Plug>(asterisk-gz#)<Cmd>lua require('hlslens').start()<CR>]], {desc="Search word under cursor"})
        end,
        opts = {
            auto_enable = false,
            -- :nohlsearch when moving out of search term.
            -- cons: disables when scrolling far since cursor moves. Doesn't disable right away on insert, only after a change.
            -- better solution made manually in ModeChanged and specific keypresses, such as Esc.
            -- calm_down = true,
        },
    },
}
