-- Taken from https://github.com/rktjmp/lush.nvim/blob/main/lua/lush/math.lua

local clamp = function(val, min, max)
  return math.min(max, math.max(min, val))
end

-- There is not math.round function so we make one.
-- round float, implementation rounds 0.5 upwards.
local round = function(val)
  return math.floor(val + 0.5)
end

return {
  round = round,
  clamp = clamp
}
