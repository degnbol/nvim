return {
    {
        'LuaSnip',
        -- the make command is optional: https://github.com/L3MON4D3/LuaSnip
        event = "InsertEnter",
        after = function()
            local luasnip = require "luasnip"
            luasnip.config.set_config {
                history = true,
                updateevents = "TextChanged,TextChangedI",
                enable_autosnippets = true,
                store_selection_keys = "<Tab>",
                delete_check_events = "TextChanged",
                ext_base_prio = 300,
                ext_prio_increase = 1,
            }
            require("luasnip.loaders.from_lua").lazy_load { paths = vim.fn.stdpath("config") .. "/luasnippets/" }
        end,
    },
}
