#!/usr/bin/env lua
-- disabled since we don't really need the info visible all the time
return {"glepnir/galaxyline.nvim", enabled=false, config=function()
    local gl = require("galaxyline")
    local gls = gl.section
    local vim = vim
    local colors = require "themes/onedark"
    local buffer = require 'galaxyline.provider_buffer'
    local fileinfo = require 'galaxyline.provider_fileinfo'

    -- empty string is for terminal
    gl.short_line_list = {"NvimTree", ""}

    local checkwidth = function()
        -- also make sure it is a file
        local squeeze_width = vim.fn.winwidth(0) / 2
        if squeeze_width > 30 then
            return true
        end
        return false
    end

    gls.left[1] = {
        FirstElement = {
            provider = function() return '▋' end,
            highlight = { colors.darker_black, colors.darker_black }
        },
    }

    gls.left[2] = {
        current_dir = {
            provider = function()
                local dir_name = vim.fn.fnamemodify(vim.fn.getcwd(), ":t")
                return "  " .. dir_name .. " "
            end,
            highlight = {colors.grey_fg2, colors.lightbg2},
            separator = " ",
            separator_highlight = {colors.lightbg2, colors.darker_black}
        }
    }

    gls.left[3] = {
        DiffAdd = {
            provider = "DiffAdd",
            condition = checkwidth,
            icon = "  ",
            highlight = {colors.white, colors.darker_black}
        }
    }

    gls.left[4] = {
        DiffModified = {
            provider = "DiffModified",
            condition = checkwidth,
            icon = "   ",
            highlight = {colors.grey_fg2, colors.darker_black}
        }
    }

    gls.left[5] = {
        DiffRemove = {
            provider = "DiffRemove",
            condition = checkwidth,
            icon = "  ",
            highlight = {colors.grey_fg2, colors.darker_black}
        }
    }

    gls.left[6] = {
        DiagnosticError = {
            provider = "DiagnosticError",
            icon = "  ",
            highlight = {colors.red, colors.darker_black}
        }
    }

    gls.left[7] = {
        DiagnosticWarn = {
            provider = "DiagnosticWarn",
            icon = "  ",
            highlight = {colors.yellow, colors.darker_black}
        }
    }

    gls.right[1] = {
        lsp_status = {
            provider = function()
                local clients = vim.lsp.get_active_clients()
                if next(clients) ~= nil then
                    return " " .. "  " .. " LSP "
                else
                    return ""
                end
            end,
            highlight = {colors.grey_fg2, colors.darker_black}
        }
    }

    gls.right[2] = {
        ViMode = {
            provider = function()
                local alias = {
                    n = "Normal",
                    i = "Insert",
                    c = "Command",
                    V = "Visual",
                    [""] = "Visual",
                    v = "Visual",
                    R = "Replace"
                }
                local current_Mode = alias[vim.fn.mode()]

                if current_Mode == nil then
                    return "  Terminal "
                else
                    return "  " .. current_Mode .. " "
                end
            end,
            highlight = {colors.red, colors.darker_black}
        }
    }

    gls.right[3] = {
        line_percentage = {
            provider = function()
                local current_line = vim.fn.line(".")
                local total_line = vim.fn.line("$")

                if current_line == 1 then
                    return "  Top "
                elseif current_line == vim.fn.line("$") then
                    return "  Bot "
                end
                local percent, _ = math.modf((current_line / total_line) * 100)
                return string.format("  %2d", percent) .. "% "
            end,
            highlight = {colors.green, colors.darker_black}
        }
    }
end}
