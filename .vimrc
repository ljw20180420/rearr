set nocompatible
set foldcolumn=1
set foldmethod=marker
set foldmarker={{{,}}}
set foldtext=''
set nowrap
set scrollopt=hor
let s:dir = fnamemodify(expand('<sfile>:p'), ':h') . "/AnsiEsc.vim-12"
let &rtp = s:dir
1split
windo set scrollbind
set statusline=reference
set laststatus=0
windo set conceallevel=3
windo set concealcursor=nvic