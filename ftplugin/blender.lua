-- Add blender stubs (fake-bpy-module) to basedpyright stubPath.
-- Triggered by compound filetype "python.blender" (set in ftdetect/blender.lua).
-- Reinstall: uv pip install --target lsp_ext/blender-stubs fake-bpy-module-latest
--
-- Creates symlinks from python_stubs/<pkg>-stubs → blender-stubs/<pkg>-stubs
-- so basedpyright's existing stubPath resolves bpy, bmesh, mathutils, etc.
-- vim.fn.stdpath has no API equivalent
local config = vim.fn.stdpath("config") ---@diagnostic disable-line: param-type-mismatch
local src = config .. "/lsp_ext/blender-stubs"
local dst = config .. "/lsp_ext/python_stubs"

-- Symlink each *-stubs/ package into the shared stubPath directory.
-- Only runs once per session (idempotent, checks existing symlinks).
local dir = vim.uv.fs_scandir(src)
if dir then
    while true do
        local name, typ = vim.uv.fs_scandir_next(dir)
        if not name then break end
        if typ == "directory" and name:match("%-stubs$") then
            local link = dst .. "/" .. name
            if not vim.uv.fs_lstat(link) then
                vim.uv.fs_symlink(src .. "/" .. name, link)
            end
        end
    end
end
