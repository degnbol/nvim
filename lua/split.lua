#!/usr/bin/env lua
-- from https://stackoverflow.com/questions/1426954/split-string-in-lua
-- string split that splits a,b,,c -> a,b,"",c rather than typical suggestion where -> a,b,c
-- USE: for field in s:split('\t') do ... end 
-- e.g. out = {} then a number of calls to table.insert(out, field) to put in a vector (table).
function string:split(pat)
    return (self .. pat):gmatch("([^"..pat.."]*)"..pat)
end
