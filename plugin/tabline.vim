if exists("+showtabline")
    function! CustomTabLine()
        let s = ''
        let cur_tab = tabpagenr()
        for i in range(tabpagenr('$'))
            let i=i+1
            let is_cur = cur_tab == i

            " Set the tab page number (for mouse clicks).
            let s ..= '%' .. i .. 'T'
            let s ..= (is_cur ? '%1*' : '%2*')
            let s ..= '%#TabNum# '

            let s ..= (is_cur ? '%#TabNumSel#' : '%#TabNum#')
            let s ..= i
            let n_wins = tabpagewinnr(i,'$')
            if n_wins > 1
                let s ..= ':'
                let s ..= (is_cur ? '%#TabWinNumSel#' : '%#TabWinNum#')
                let s ..= (tabpagewinnr(i,'$') > 1 ? n_wins : '')
            end

            let s ..= ' %*'
            let s ..= (is_cur ? '%#TabLineSel#' : '%#TabLine#')
            let bufnr = tabpagebuflist(i)[tabpagewinnr(i) - 1]
            let file = bufname(bufnr)
            let buftype = getbufvar(bufnr, 'buftype')
            if buftype == 'nofile'
                if file =~ '\/.'
                    let file = substitute(file, '.*\/\ze.', '', '')
                endif
            else
                let file = fnamemodify(file, ':p:t')
            endif
            if file == ''
                let file = '[No Name]'
            endif
            let s ..= file
            let s ..= (is_cur ? '%m' : '')
        endfor
        let s ..= '%T%#TabLineFill#%='
        return s
    endfunction
    set stal=2
    set tabline=%!CustomTabLine()
endif
