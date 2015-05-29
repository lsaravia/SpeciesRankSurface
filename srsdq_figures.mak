all:srsdq_figures.pdf srsdq_figures.docx 

srsdq_figures.pdf: srsdq_figures.sty srsdq_figures.md
	pandoc -H srsdq_figures.sty srsdq_figures.md -o srsdq_figures.pdf 
	evince srsdq_figures.pdf		

srsdq_figures.docx: srsdq_figures.sty srsdq_figures_MEE.md
	pandoc -H srsdq_figures.sty srsdq_figures_MEE.md -o srsdq_figures.docx

