local function getTopline()
    return vim.fn.line("w0")
end

local function setTopline(topline)
    vim.fn.winrestview {topline=topline, lnum=topline+vim.opt_local.scrolloff:get()}
end

function vnewFurther()
    for i = 1, math.max(vim.v.count,1) do
        win = vim.api.nvim_get_current_win()
        vim.opt_local.scrollbind = true
        topline = getTopline()
        height = vim.api.nvim_win_get_height(win)
        vim.cmd.vnew "%"
        setTopline(topline+height)
        vim.opt_local.scrollbind = true
    end
end

function forwardScreen()
    setTopline(getTopline() + vim.api.nvim_win_get_height(win))
end
function backwardScreen()
    setTopline(getTopline() - vim.api.nvim_win_get_height(win))
end


