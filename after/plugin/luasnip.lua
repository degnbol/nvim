local rtp = vim.opt.runtimepath:get()[1]
require("luasnip.loaders.from_lua").lazy_load { paths = rtp .. "/luasnippets/" }
