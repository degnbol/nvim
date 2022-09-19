#!/usr/bin/env lua
local cmpd = require "cmp_dictionary"
local rtp = vim.opt.runtimepath:get()[1]

cmpd.setup {
    dic = {
        -- dicts generated with ./spell.sh
        ["*"] = {
            rtp .. "/spell/custom.dic",
            rtp .. "/spell/en.dic",
        },
        spelllang = { 
            da = rtp .. "/spell/da.dic",
        }
    },
    first_case_insensitive = true,
}

-- We want spelllang=en,da so we can underline bad spelling in both, 
-- but toggle completion from danish only when iminsert=1.
function CmpDictUpdate()
    if vim.bo.iminsert == 1 then
        cmpd.update()
    else
        local spelllang = vim.bo.spelllang
        vim.bo.spelllang = "en"
        cmpd.update()
        vim.bo.spelllang = spelllang
    end
end

-- not sure why defer works but simply calling the function doesn't.
vim.defer_fn(CmpDictUpdate, 0)

