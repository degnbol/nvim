local M = {}

function M.set(name, val)
    vim.api.nvim_set_hl(0, name, val)
end

function M.setfg(name, color)
    M.set(name, { fg = color })
end

function M.setbg(name, color)
    M.set(name, { bg = color })
end

function M.fgbg(name, fgcol, bgcol)
    M.set(name, { fg = fgcol, bg = bgcol })
end

function M.link(name, linkto)
    M.set(name, { link = linkto })
end

function M.def(name, linkto)
    M.set(name, { link = linkto, default = true })
end

function M.rev(name)
    M.set(name, { reverse = true })
end

function M.get(name)
    return vim.api.nvim_get_hl(0, { name = name, link = false })
end

function M.fg(name)
    return vim.api.nvim_get_hl(0, { name = name, link = false })['fg']
end

function M.bg(name)
    return vim.api.nvim_get_hl(0, { name = name, link = false })['bg']
end

-- update subset of settings for a highlight group instead of replacing them all
function M.mod(name, val)
    M.set(name, vim.tbl_extend("force", M.get(name), val))
end

function M.clear(name)
    vim.api.nvim_set_hl(0, name, {})
end

function M.hide(name)
    vim.api.nvim_set_hl(0, name, { fg = "bg" })
end

---Intended for adding highlights relevant to specific filetype or plugin after
---a colorscheme is set.
---@param callback function
---@return integer
function M.afterColorscheme(callback)
    return vim.api.nvim_create_autocmd("Colorscheme", {
        pattern = "*",
        group = vim.api.nvim_create_augroup("afterColorscheme", { clear = false }),
        callback = callback
    })
end

return M
