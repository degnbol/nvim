local map = require "utils/keymap"

-- Compile/build keymap for JSON files
map.n("<leader>cc", function()
    local file = vim.api.nvim_buf_get_name(0)
    if file:match("complex_modifications/.*%.json$") then
        vim.system({ vim.fn.expand("~/.config/karabiner/karabiner.json.sh") }, {}, function(obj)
            vim.schedule(function()
                if obj.code == 0 then
                    vim.notify("Regenerated karabiner.json")
                else
                    vim.notify("Error: " .. (obj.stderr or "unknown"), vim.log.levels.ERROR)
                end
            end)
        end)
    else
        vim.notify("No compile action for this JSON file", vim.log.levels.WARN)
    end
end, "\"Compile\" karabiner.json", { buffer = true})
