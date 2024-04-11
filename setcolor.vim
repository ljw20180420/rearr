:hi MyRedFg ctermfg=red
:hi MyGreenFg ctermfg=green
:hi MyYellowFg ctermfg=yellow
:hi MyBlueBg ctermbg=blue
:hi MyNormalBg
:syn region MyRedFg matchgroup=MyRedFgMat start=/\e\[31m/ end=/\e\[39m/ oneline concealends
:syn region MyGreenFg matchgroup=MyGreenFgMat start=/\e\[32m/ end=/\e\[39m/ oneline concealends containedin=ALL
:syn region MyYellowFg matchgroup=MyYellowFgMat start=/\e\[33m/ end=/\e\[39m/ oneline concealends containedin=ALL
:syn match MyBlueBg /([acgtnACGTN ])/ conceal cchar="n" containedin=ALL
:syn region MyNormalBg matchgroup=MyNormalBgMat start=/{/ skip=/[acgtnACGTN ]/ end=/}/ oneline concealends containedin=ALL

