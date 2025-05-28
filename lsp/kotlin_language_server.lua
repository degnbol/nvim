#!/usr/bin/env lua
return {
    -- assume project root is at git root.
    -- Some imports don't work if unset.
    -- If you have a project without git root at project root then 
    -- search upwards for something else with this function.
    root_dir = function() return vim.env.ROOT end
}
