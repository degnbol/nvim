-- author: "Thomas Leon Highbaugh thighbaugh@zoho.com"
require "base16"({
    base00 = "000000";
    base01 = "262a31";
    base02 = "3c3f4c";
    base03 = "454754";
    base04 = "555e70";
    base05 = "8b9cbe";
    base06 = "B2BFD9";
    base07 = "f4f4f7";
    base08 = "ff29a8";
    base09 = "85ffe0";
    base0A = "f0ffaa";
    base0B = "0badff";
    base0C = "8265ff";
    base0D = "00eaff";
    base0E = "00ffcc";
    base0F = "ff3d81";
}, true)
-- increased some grays.
-- Also just making comment brighter
vim.api.nvim_set_hl(0, "Comment", {fg="gray"})
