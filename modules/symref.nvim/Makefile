.deps/plenary.nvim:
	git clone --depth 1 https://github.com/nvim-lua/plenary.nvim .deps/plenary.nvim

# Plenary starts each spec's nvim without --clean; its own XDG dirs keep user
# config, data and logs out.
test: export XDG_CONFIG_HOME = $(CURDIR)/.deps/xdg/config
test: export XDG_DATA_HOME = $(CURDIR)/.deps/xdg/data
test: export XDG_STATE_HOME = $(CURDIR)/.deps/xdg/state
test: .deps/plenary.nvim
	nvim --headless -u tests/minimal_init.lua -c "PlenaryBustedDirectory tests/ {minimal_init = 'tests/minimal_init.lua'}"

.PHONY: test
