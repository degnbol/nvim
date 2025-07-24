local hi = require "utils/highlights"

-- Clear highlights for now.
-- Might be useful syntax captures.
hi.afterColorscheme(function()
    hi.clear("csvCol1")
    hi.clear("csvCol2")
    hi.link("csvCol3", "csvCol1")
    hi.clear("csvCol4")
    hi.link("csvCol5", "csvCol1")
    hi.clear("csvCol6")
    hi.link("csvCol7", "csvCol1")
    hi.clear("csvCol8")
    hi.link("csvCol9", "csvCol1")
    hi.clear("csvCol10")
    hi.link("csvCol11", "csvCol1")
end)
