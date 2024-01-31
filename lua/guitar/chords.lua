--- Instead of using vim.expand("<cword>") we find a word under cursor with any 
--- pattern of allowed chars.
--- line: string of the current line
--- c: 0-indexed column of cursor.
--- pattern: string of allowed chars, e.g. "%w()" for alphanumeric as well as parenthesis.
local function cword(line, c, pattern)
    return line:sub(1,c):match("["..pattern.."]*$") .. 
           line:sub(c+1):match("^["..pattern.."]*")
end

local function readjson(filename)
    file = io.open(filename, "r")
    content = vim.json.decode(file:read("*a"))
    file:close()
    return content
end

local rtp = vim.opt.runtimepath:get()[1]
local chordchart = readjson(rtp .. "/lua/guitar/chordchart.json")
local name2strings = {}
for _, chord in ipairs(chordchart) do
    name2strings[chord.name] = {chord.strings, chord.fret}
    -- add lowercase single char versions of minor chords
    if chord.name:match("^[A-Z]m$") then
        name2strings[chord.name:sub(1,1):lower()] = name2strings[chord.name]
    end
end

--- slash chords. As usual there are multiple ways to make one.
--- I go for the simplest one to code which may not be easy to play.
--- name: e.g. "C/G"
--- returns: {strings, fret}
function slashChord(name)
    root, bass = name:match("(.*)/(.*)")
    root = name2strings[root]
    bass = name2strings[bass]
    if root == nil or bass == nil then return end
    root, rootFret = unpack(root)
    bass, bassFret = unpack(bass)
    -- change to using the lowest note in bass and it should be the lowest 
    -- played note.
    lowest = tostring(bass):match("X*%d")
    -- account for any difference in offset
    if rootFret ~= bassFret then
        note = tonumber(lowest:sub(-1,-1)) - bassFret + rootFret
        lowest = lowest:sub(1,-2) .. tostring(note)
    end
    return lowest .. tostring(root):sub(#lowest+1), rootFret
end

--- Convert e.g. x00231 to pretty version in guitar notes:
--- ╳││││●
--- │││●││
--- ││││●│
local function prettyChord(strings, fret)
    strings = tostring(strings)
    -- default to 1, which isn't displayed for simplicity
    fret = fret or 1
    length = 4
    for d in strings:gmatch("%d") do
        length = math.max(length, d)
    end
    pretty = {}
    for i = 1, length do
        pretty[i] = {"│","│","│","│","│","│"}
    end
    for i = 1, 6 do
        d = tonumber(strings:sub(i,i))
        if d == nil then
            pretty[1][i] = "╳"
        elseif d ~= 0 then
            pretty[d][i] = "●"
        end
    end
    if fret ~= 1 then pretty[1][7] = fret end
    for i = 1, length do
        pretty[i] = table.concat(pretty[i])
    end
    return pretty
end

--- Detect a chord name or strings pattern (e.g. x00231) on current line and 
--- return {strings, fret}
local function detectChordLine()
    r, c = unpack(vim.api.nvim_win_get_cursor(0))
    line = vim.api.nvim_get_current_line()
    strings = line:match(("[%dxX-]"):rep(6))
    if strings ~= nil then
        return strings:gsub('-', '0'), 1
    else
        name = cword(line, c, "%w/()")
        print(name)
        if name:match("/") then
            return slashChord(name)
        end
        chord = name2strings[name]
        if chord == nil then
            -- may be written as single chars, e.g. "AeDA" for "A Em D A"
            char = line:sub(c+1,c+1)
            chord = name2strings[char]
            if chord == nil then return end
        end
        return unpack(chord)
    end
end

--- Convert e.g. "x00231" on current line OR "Am" as current word to pretty 
--- version in guitar notes:
--- ╳││││●
--- │││●││
--- ││││●│
function prettyChordLine()
    strings, fret = detectChordLine()
    if strings == nil then
        print("Chord code not understood.")
        return
    end
    pretty = prettyChord(strings, fret)
    -- place under
    vim.api.nvim_buf_set_lines(0, r, r, true, pretty)
end

vim.api.nvim_create_user_command("Chord", prettyChordLine, {})

