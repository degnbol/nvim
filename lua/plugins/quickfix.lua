-- better quickfix window. zf to open fzf inside quickfix.
return {'kevinhwang91/nvim-bqf', dependencies={'junegunn/fzf'}, config=function ()
    local fn = vim.fn
    local max = math.max
    local min = math.min

    function _G.qftf(info)
        local items
        local ret = {}
        if info.quickfix == 1 then
            items = fn.getqflist({id = info.id, items = 0}).items
        else
            items = fn.getloclist(info.winid, {id = info.id, items = 0}).items
        end
        -- find longest fname that will be printed
        local max_fname_length = 0
        for i = info.start_idx, info.end_idx do
            local e = items[i]
            if e.valid == 1 then
                if e.bufnr > 0 then
                    fname = fn.bufname(e.bufnr)
                    max_fname_length = max(max_fname_length, fname:len())
                end
            end
        end
        -- 25 is max before truncation
        local limit = min(max_fname_length, 25)
        local fnameFmt1 = '%-'  ..  limit      .. 's'
        local fnameFmt2 = '…%.' .. (limit - 1) .. 's'
        local validFmt = '%s│%5d:%-3d│%s %s'
        for i = info.start_idx, info.end_idx do
            local e = items[i]
            local fname = fnameFmt1:format('') -- default
            local str
            if e.valid == 1 then
                if e.bufnr > 0 then
                    fname = fn.bufname(e.bufnr)
                    if fname == '' then
                        fname = '[No Name]'
                    else
                        -- show home dir as just a tilde for clarity.
                        fname = fname:gsub('^' .. vim.env.HOME, '~')
                    end
                    -- char in fname may occur more than 1 width, ignore this issue in order to keep performance
                    if #fname <= limit then
                        fname = fnameFmt1:format(fname)
                    else
                        fname = fnameFmt2:format(fname:sub(1 - limit))
                    end
                end
                local lnum = e.lnum > 99999 and -1 or e.lnum
                local col = e.col > 999 and -1 or e.col
                local qtype = e.type == '' and ' ' or e.type:sub(1, 1):upper()
                str = validFmt:format(fname, lnum, col, qtype, e.text)
            else
                str = e.text
            end
            table.insert(ret, str)
        end
        return ret
    end

    vim.o.qftf = '{info -> v:lua._G.qftf(info)}'

    -- Adapt fzf's delimiter in nvim-bqf
    require('bqf').setup({
        filter = {
            fzf = {
                extra_opts = {'--bind', 'ctrl-o:toggle-all', '--delimiter', '│'}
            }
        }
    })
end}
