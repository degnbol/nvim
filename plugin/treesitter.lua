-- use the zsh parser for all shell filetypes. tree-sitter-bash is unmaintained
-- and has cascading error bugs (parameter expansion, extglob, heredocs).
-- zsh is a superset of bash, so the zsh parser handles bash files fine.
vim.treesitter.language.register("zsh", "sh")
vim.treesitter.language.register("zsh", "bash")
vim.treesitter.language.register("zsh", "sh.zsh")
-- ```scm fenced code blocks in markdown should use the query parser
vim.treesitter.language.register("query", "scm")

-- Custom injection directive: `(#trim! @cap PREFIX_BYTES SUFFIX_BYTES)`
-- Skips the first PREFIX_BYTES (newline-aware: increments row, resets col)
-- and the last SUFFIX_BYTES (assumed not to cross a newline) of @cap, then
-- sets metadata.range to a (row, col, byte) triple consistent across all
-- three coordinates. `#offset!` does naive (row+drow, col+dcol) arithmetic
-- and breaks when col_delta would push col past the end of a short line —
-- e.g. injecting lua into `nvim -c 'lua\nCODE\n'` where the first line of
-- the raw_string ends right after `'lua`.
vim.treesitter.query.add_directive("trim!", function(match, _, source, pred, metadata)
    local capture_id = pred[2]
    local prefix_bytes = tonumber(pred[3]) or 0
    local suffix_bytes = tonumber(pred[4]) or 0
    local nodes = match[capture_id]
    if not nodes or #nodes == 0 then return end
    local node = nodes[1]
    local sr, sc, sb, er, ec, eb = node:range(true)
    local text = vim.treesitter.get_node_text(node, source)
    local new_sr, new_sc, new_sb = sr, sc, sb
    for i = 1, prefix_bytes do
        if text:byte(i) == 0x0a then
            new_sr = new_sr + 1
            new_sc = 0
        else
            new_sc = new_sc + 1
        end
        new_sb = new_sb + 1
    end
    if not metadata[capture_id] then metadata[capture_id] = {} end
    metadata[capture_id].range = {
        new_sr, new_sc, new_sb,
        er, ec - suffix_bytes, eb - suffix_bytes,
    }
end, { force = true })

-- start treesitter for each new filetype
vim.api.nvim_create_autocmd("FileType", {
    pattern = "*",
    group = vim.api.nvim_create_augroup("start_treesitter", { clear = false }),
    callback = function(args)
        local disabled = {
            -- messes with vimtex in lots of ways, e.g. conceal, detection of mathzone, cycling with ts$,
            "latex", "plaintex", "tex",
        }
        local additional_vim_regex_highlighting = {
            "vimdoc",   -- treesitter version doesn't contain useful colors from :h group-name
            "sh", "bash", "sh.zsh",
            "markdown", -- my custom comment syntax matches in after/syntax/markdown.vim
            -- Semicolon isn't currently highlighted in all cases by TS so we want to incl vim regex hl for jl.
            -- However, jl can get slowed down a lot in certain files from the syntax hl. The solution:
            -- We enable it, but avoid any default syntax hl and only set custom syntax hl in syntax/julia.vim.
            "julia",
            "sql",  -- custom postgres highlight in syntax/sql.vim
            "wgsl", -- custom in syntax/wgsl.vim
        }
        if vim.b[args.buf].largefile then return end
        if not vim.list_contains(disabled, vim.bo.filetype) then
            -- WORKAROUND: nvim 0.12 bundled markdown parser crashes during initial load.
            -- Delay treesitter.start for markdown to after buffer is fully set up.
            if vim.bo.filetype == "markdown" then
                vim.schedule(function()
                    pcall(vim.treesitter.start, args.buf)
                end)
            else
                pcall(vim.treesitter.start, args.buf)
            end
            if vim.list_contains(additional_vim_regex_highlighting, vim.bo.filetype) then
                vim.bo[args.buf].syntax = 'on'
            end
        end
    end
})
