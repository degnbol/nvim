hi def link hyprCommand Statement
syn match Operator /=/ containedin=hyprCategory
syn match Operator /\$\ze[A-Za-z0-9]/
syn keyword HyprKey SHIFT ALT SUPER
hi def link HyprKey Constant
syn match Delimiter /,/ containedin=hyprValue
syn keyword Type rgba containedin=hyprCategory nextgroup=rgbaColor
syn region rgbaColor matchgroup=Delimiter start=/(/ end=/)/ containedin=hyprCategory
hi def link rgbaColor Special
syn keyword Boolean true false containedin=hyprCategory
syn match Number /[-+]\?\<[0-9]\+\.\?[0-9]*/ containedin=hyprCategory
