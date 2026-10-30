
local M = {}

---Is this terminal kitty?
function M.term()
    return vim.env.KITTY_PID ~= nil
end

---Cell layout. All numbers except `aspect` are fractions of the cell height.
---@class kitty.CellFractions
---@field postscript_name string medium face
---@field aspect number cell height / cell width
---@field below number baseline to cell bottom
---@field x number ink height of "x"
---@field cap number ink height of "H"
---@field desc number ink height of "p" minus that of "x"

---`text_composition_strategy` as kitty parses and clamps it.
---@class kitty.TextComposition
---@field gamma number at least 0.01
---@field contrast number percent, 0 to 100

local launched = {} ---@type table<string, any> decoded output by script name

---Decoded JSON stdout of a config `scripts/` file run with `kitty +launch`.
---Runs synchronously, and errors if `kitty` is not on PATH.
---A success is cached for the session, a failure is not.
---@param name string file name in the config's `scripts/`
---@return any|nil decoded
---@return string|nil err kitty's stderr on failure
local function launch_json(name)
    if launched[name] == nil then
        local script = vim.fn.stdpath("config") .. "/scripts/" .. name
        local res = vim.system({ "kitty", "+launch", script }, { text = true }):wait()
        if res.code ~= 0 then return nil, res.stderr end
        launched[name] = vim.json.decode(res.stdout)
    end
    return launched[name]
end

---Kitty's cell layout for the configured font. Independent of zoom.
---Runs kitty synchronously, and errors if `kitty` is not on PATH.
---A success is cached for the session, a failure is not.
---@return kitty.CellFractions|nil fractions
---@return string|nil err kitty's stderr on failure
function M.cell_fractions()
    return launch_json("kitty_cell_fractions.py")
end

---Kitty's configured `text_composition_strategy`. Fails on `legacy` or a
---malformed value.
---Runs kitty synchronously, and errors if `kitty` is not on PATH.
---A success is cached for the session, a failure is not.
---@return kitty.TextComposition|nil composition
---@return string|nil err kitty's stderr on failure
function M.text_composition()
    return launch_json("kitty_text_composition.py")
end

---ImageMagick `-fx` expression that curves coverage alpha as kitty does for a
---glyph before blending (`foreground_contrast_new` in kitty's `cell.frag.glsl`):
---a′ = clamp((a + (a^(1/γ) − a)·w)·(1 + c/100), 0, 1), w = (1 − L_fg + L_bg) / 2,
---with a the pixel's alpha (`u` in `-fx`), a′ the curved alpha, γ and c
---`composition`'s gamma and contrast, and L_fg, L_bg the two luminances.
---@param composition kitty.TextComposition
---@param fg_luminance number relative luminance of the text colour, 0 to 1
---@param bg_luminance number relative luminance of the background colour, 0 to 1
---@return string fx
function M.text_alpha_fx(composition, fg_luminance, bg_luminance)
    local w = (1 - fg_luminance + bg_luminance) * 0.5
    return ("min(max((u+(pow(u,%.6f)-u)*%.6f)*%.6f,0),1)"):format(
        1 / composition.gamma, w, 1 + composition.contrast / 100)
end

---Enable/diable ligatures for the current kitty window.
---@param enable boolean
---@return boolean success
function M.ligatures(enable)
    if enable then
        return os.execute("kitty @ disable-ligatures never") or false
    else
        return os.execute("kitty @ disable-ligatures always") or false
    end
end

---Enable/disable ligatures for the current kitty window when a buffer with given filename pattern is in focus.
---E.g. set it for a filetype.
---@param enable boolean
---@param pattern string
---@return boolean success
function M.ligatures_pattern(enable, pattern)
    if not M.term() then return false end
    local grp = vim.api.nvim_create_augroup("kitty_ligatures", {clear=true})
    vim.api.nvim_create_autocmd("BufEnter", {
        pattern = pattern,
        group = grp,
        callback = function () M.ligatures(enable) end
    })
    vim.api.nvim_create_autocmd({"BufLeave", "BufWinLeave", "BufDelete"}, {
        pattern = pattern,
        group = grp,
        callback = function () M.ligatures(not enable) end
    })
    return true
end

return M
