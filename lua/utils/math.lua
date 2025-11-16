local M = {}

---Clamp a value between a min and max.
---@param val number
---@param min number
---@param max number
---@return number
function M.clamp(val, min, max)
  return math.min(max, math.max(min, val))
end

---Rounds a float, implementation rounds 0.5 upwards.
---There is no math.round function so we make one.
---@param val number
---@return integer
function M.round(val)
  return math.floor(val + 0.5)
end

---Dot product between two vectors.
---Assumes same length.
---@param v1 table
---@param v2 table
---@return number
function M.dot(v1, v2)
    local ret = 0
    for i, v in ipairs(v1) do
        ret = ret + v * v2[i]
    end
    return ret
end

---Matrix multiplication onto a vector.
---@param mat table array of rows.
---@param vec table array of numbers.
---@return table array of numbers.
function M.mul_mat_vec(mat, vec)
    local ret = {}
    for i, row in ipairs(mat) do
        ret[i] = M.dot(row, vec)
    end
    return ret
end

return M
