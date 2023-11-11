#!/usr/bin/env julia
using GZip
using JSON

raw = GZip.open("plotlyjs.txt.gz") do io readlines(io) end

indices = findall(startswith.(raw, "Parent: ")) .- 1

names = raw[indices]
parents = [split(only(s[2:end]), '.') for s in split.(raw[indices .+ 1], ' ')]
types = [s[end] for s in split.(raw[indices .+ 2], ' '; limit=2)]
defaults = [s[1] == "Default:" ? s[end] : "" for s in split.(raw[indices .+ 3], ' '; limit=2)]
descIndices = indices .+ 3 .+ (defaults .!= "")
desc = raw[descIndices]
# parents have child entries right away instead of a line with desc
desc[descIndices .∈ Ref(indices)] .= ""

for (i, p) in enumerate(parents)
    m = match(r"data\[type=(\w+)\]", p[1])
    if m ≠ nothing
        parents[i][1] = m[1]
    end
end

parents = [rstrip.(p, Ref(['[', ']'])) for p in parents]

tree = Dict{String,Dict}()
for ps in unique(parents)
    t = tree
    for p in ps
        t = t[p] = Dict{String,Any}()
    end
end

for (i, name) in enumerate(names)
    t = tree
    for p in parents[i]
        if p ∉ keys(t) t[p] = Dict{String,Any}() end
        t = t[p]
    end
    if "items" ∉ keys(t) t["items"] = Dict{String,String}[] end
    push!(t["items"], Dict("label"=>name, "documentation"=>"$(desc[i])\nDefault: $(defaults[i])\nType: $(types[i])"))
end

open("plotlyjs.json", "w") do io
    JSON.print(io, tree, 2)
end
