all:mpnc_appendices.pdf

mpnc_appendices.pdf: mpnc_appendices.sty mpnc_appendices.md
	pandoc -H mpnc_appendices.sty mpnc_appendices.md -o mpnc_appendices.pdf 
	evince mpnc_appendices.pdf		
