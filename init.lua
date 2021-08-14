-- setup heavily inspired by NvChad (https://github.com/siduck76/NvChad)

require 'settings'
require 'packages'
require 'keymappings'
require 'colors'
require 'TrueZen'
require 'whichkey'
require 'top-bufferline'
require 'statusline'

-- Languages
require'nvim-lspconfig'
require 'nvim-compe'
require'lsp_signature'.on_attach()
require 'treesitter-nvim'
require'twilight'.setup{dimming={alpha=0.5}, context=30}
require('lspkind').init()
require 'gitsigns-nvim' -- git decoration
require('nvim_comment').setup()
-- require 'indent-blankline'
require 'shebang-nvim'
require 'cheat' -- :Cheat as alternative to StackOverflow

require 'iron-nvim'
-- require 'repl'
-- require 'ripple'

require 'telescope-nvim'
require 'dashboard'
require 'tree'

require 'highlights'
require 'file-icons'



