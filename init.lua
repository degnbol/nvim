-- setup heavily inspired by NvChad (https://github.com/siduck76/NvChad)

require 'settings'
-- if you switch between package managers then remove everything from the other one with rm -rf ~/.local/share/nvim/site/pack/packer or paqs
-- require 'packages'
require 'plugins'
require 'keymappings'
require 'colors'
require 'TrueZen'
require 'whichkey'
require 'statusline'

-- Languages
require 'lsp'
require 'treesitter-nvim'
require'twilight'.setup{dimming={alpha=0.5}, context=30}
require('lspkind').init()
require 'gitsigns-nvim' -- git decoration
require('nvim_comment').setup()
require 'shebang-nvim'

require 'telescope-nvim'
require 'dashboard'
require 'tree'

require 'highlights'
require 'file-icons'

