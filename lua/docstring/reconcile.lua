-- Reconcile a parsed docstring against its signature: emit documented params in
-- signature order, add fresh entries for new params, and surface documented
-- params absent from the signature as `unmatched` orphans — never dropping a
-- hand-written description. Summary, returns, raises, and unmodeled sections
-- pass through untouched (returns/raises refresh is a later milestone).

local model_mod = require "docstring.model"

local M = {}

--- Refresh a matched entry's type from a signature annotation, honouring the
--- typeless convention: a type is only refreshed on an entry that *already*
--- carries one (code is the source of truth), and never injected into a typeless
--- entry — this config's Google style keeps types in the annotated signature. An
--- unannotated signature never blanks an existing hand-written type.
--- @param entry docstring.Param  the documented entry (mutated in place)
--- @param sig docstring.SigParam
local function refresh_type(entry, sig)
    if entry.type and sig.type then entry.type = sig.type end
end

--- Reconcile a DocModel's params against the ordered signature params, returning
--- a new model. Params are emitted in signature order; a documented entry
--- (whether currently a param or a surfaced orphan) is reused by name — keeping
--- its description and, subject to `refresh_type`, its type — while a signature
--- param with no documented entry gets a fresh typeless entry with an empty
--- description. Documented entries matching no signature param become
--- `unmatched` orphans, in their original order. Summary, returns, raises, and
--- other sections are carried through unchanged.
--- @param params docstring.SigParam[]  ordered signature params
--- @param model docstring.Model
--- @return docstring.Model
function M.reconcile(params, model)
    -- Pool of documented entries keyed by name, in first-seen order. Orphans
    -- re-enter the pool so a re-added param reclaims its old description rather
    -- than starting blank; a live param shadows an orphan of the same name.
    -- A docstring listing the same name twice is malformed; last wins.
    local pool = {}
    local order = {}
    for _, entry in ipairs(model.params) do
        if not pool[entry.name] then order[#order + 1] = entry.name end
        pool[entry.name] = entry
    end
    for _, entry in ipairs(model.unmatched) do
        if not pool[entry.name] then
            pool[entry.name] = entry
            order[#order + 1] = entry.name
        end
    end

    local result = model_mod.empty()
    result.summary = model.summary
    result.returns = model.returns
    result.raises = model.raises
    result.other = model.other

    local used = {}
    for _, sig in ipairs(params) do
        local entry = pool[sig.name]
        if entry then
            -- Copy before refreshing: the pooled table is the caller's, and a
            -- consumer may still hold the input model (rename matcher, hover).
            entry = { name = entry.name, type = entry.type, description = entry.description }
            refresh_type(entry, sig)
            used[sig.name] = true
        else
            --- @type docstring.Param
            entry = { name = sig.name, type = nil, description = "" }
        end
        result.params[#result.params + 1] = entry
    end
    for _, name in ipairs(order) do
        if not used[name] then
            result.unmatched[#result.unmatched + 1] = pool[name]
        end
    end
    return result
end

return M
