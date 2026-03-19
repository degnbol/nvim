rtps/plenary.nvim:
	git clone --depth 1 https://github.com/nvim-lua/plenary.nvim rtps/plenary.nvim

test-deps: rtps/plenary.nvim

test: test-deps
	nvim --headless -u tests/minimal_init.lua -c "PlenaryBustedDirectory tests/plenary/ {minimal_init = 'tests/minimal_init.lua'}"

.PHONY: test-deps test
