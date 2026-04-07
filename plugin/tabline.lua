function _G.CustomTabLine()
    local s = {}
    local cur_tab = vim.fn.tabpagenr()
    for i = 1, vim.fn.tabpagenr('$') do
        local is_cur = cur_tab == i
        -- tab page number for mouse clicks
        s[#s + 1] = '%' .. i .. 'T'
        s[#s + 1] = is_cur and '%1*' or '%2*'
        s[#s + 1] = '%#TabNum# '
        s[#s + 1] = is_cur and '%#TabNumSel#' or '%#TabNum#'
        s[#s + 1] = i
        local n_wins = vim.fn.tabpagewinnr(i, '$')
        if n_wins > 1 then
            s[#s + 1] = ':'
            s[#s + 1] = is_cur and '%#TabWinNumSel#' or '%#TabWinNum#'
            s[#s + 1] = n_wins
        end
        s[#s + 1] = ' %*'
        s[#s + 1] = is_cur and '%#TabLineSel#' or '%#TabLine#'
        local bufnr = vim.fn.tabpagebuflist(i)[vim.fn.tabpagewinnr(i)]
        local file = vim.fn.bufname(bufnr)
        if vim.fn.getbufvar(bufnr, 'buftype') == 'nofile' then
            file = file:match('.*/(.+)') or file
        else
            file = vim.fn.fnamemodify(file, ':p:t')
        end
        if file == '' then file = '[No Name]' end
        s[#s + 1] = file
        if is_cur then s[#s + 1] = '%m' end
    end
    s[#s + 1] = '%T%#TabLineFill#%='
    return table.concat(s)
end

vim.o.tabline = '%!v:lua.CustomTabLine()'
