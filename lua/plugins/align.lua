
-- mini.align is lua version of junegunn/vim-easy-align with config
-- require'mini.align'.setup {}
return {
    -- the only plugin I could get to work with auto aligning tables in latex.
    {"junegunn/vim-easy-align", config=function ()
        -- ga is normally for seeing info about char under cursor, but I rarely use it.
        -- Just use gA for that instead, set in init.lua for tpope/vim-characterize.
        vim.keymap.set({"n", "x"}, "ga", "<Plug>(EasyAlign)", { desc="Align" })
    end},
    -- "godlygeek/tabular", -- one of the first ones. No operator motion 
    -- use which is needed to align entire table env in latex.
    -- "godlygeek/tabular",
    {"Vonr/align.nvim", enabled=false, config=function ()
        -- copy of default config
-- align_to_char(length, reverse, preview, marks)
-- align_to_string(is_pattern, reverse, preview, marks)
-- align(str, reverse, marks)
-- operator(fn, opts)
-- Where:
--      length: integer
--      reverse: boolean
--      preview: boolean
--      marks: table (e.g. {1, 0, 23, 15})
--      str: string (can be plaintext or Lua pattern if is_pattern is true)
local NS = { noremap = true, silent = true }
vim.keymap.set('x', 'aa', function() require'align'.align_to_char(1, true)             end, NS) -- Aligns to 1 character, looking left
vim.keymap.set('x', 'as', function() require'align'.align_to_char(2, true, true)       end, NS) -- Aligns to 2 characters, looking left and with previews
vim.keymap.set('x', 'aw', function() require'align'.align_to_string(false, true, true) end, NS) -- Aligns to a string, looking left and with previews
vim.keymap.set('x', 'ar', function() require'align'.align_to_string(true, true, true)  end, NS) -- Aligns to a Lua pattern, looking left and with previews
-- Example gawip to align a paragraph to a string, looking left and with previews
vim.keymap.set( 'n', 'gaw', function()
        local a = require'align'
        a.operator( a.align_to_string, { is_pattern = false, reverse = true, preview = true })
    end, NS)
-- Example gaaip to aling a paragraph to 1 character, looking left
vim.keymap.set( 'n', 'gaa', function()
        local a = require'align'
        a.operator( a.align_to_char, { length = 1, reverse = true })
    end, NS)
end}

}
