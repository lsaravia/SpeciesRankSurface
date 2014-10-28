all:srsdq_figures.pdf

srsdq_figures.pdf: srsdq_figures.sty srsdq_figures.md
	pandoc -H srsdq_figures.sty srsdq_figures.md -o srsdq_figures.pdf 
	evince srsdq_figures.pdf		
