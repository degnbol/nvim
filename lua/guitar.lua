
--- Convert e.g. x00231 on current line to pretty version in guitar notes:
--- ╳││││●
--- │││●││
--- ││││●│
vim.api.nvim_buf_create_user_command(0, "Guitar", function()
    r, c = unpack(vim.api.nvim_win_get_cursor(0))
    line = vim.api.nvim_get_current_line()
    short = line:match(("[%dxX-]"):rep(6))
    if short == nil then return end
    short = short:gsub('-', '0')
    length = 3
    for d in short:gmatch("%d") do
        length = math.max(length, d)
    end
    long = {}
    for i = 1, length do
        long[i] = {"│","│","│","│","│","│"}
    end
    for i = 1, 6 do
        d = tonumber(short:sub(i,i))
        if d == nil then
            long[1][i] = "╳"
        elseif d ~= 0 then
            long[d][i] = "●"
        end
    end
    for i = 1, length do
        long[i] = table.concat(long[i])
    end
    vim.api.nvim_buf_set_lines(0, r-1, r, true, long)
end, {})
