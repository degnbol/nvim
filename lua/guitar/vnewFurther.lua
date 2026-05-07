local utils = require "utils"

local function getTopline()
    return vim.fn.line("w0")
end

function vnewFurther()
    for i = 1, vim.v.count1 do
        win = vim.api.nvim_get_current_win()
        vim.opt_local.scrollbind = true
        topline = getTopline()
        height = vim.api.nvim_win_get_height(win)
        vim.cmd.vnew "%"
        utils.set_view { topline = topline + height }
        vim.opt_local.scrollbind = true
    end
end

function forwardScreen()
    utils.set_view {
        topline = getTopline() + vim.api.nvim_win_get_height(win),
    }
end

function backwardScreen()
    utils.set_view {
        topline = getTopline() - vim.api.nvim_win_get_height(win),
    }
end
