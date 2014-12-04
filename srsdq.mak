all:srsdq.docx

srsdq.bib: srsdq.md
	cp "/home/leonardo/BibTeX/Manuscritos-SAD vs SRS-Dq.bib" srsdq.bib

srsdq.pdf: srsdq.md margins.sty srsdq.bib srsdq.mak
	pandoc -H margins.sty --bibliography srsdq.bib --csl=plos.csl srsdq.md -o srsdq.pdf 
	evince srsdq.pdf		

srsdq.docx: srsdq.md margins.sty srsdq.bib srsdq.mak
	pandoc -H margins.sty --bibliography srsdq.bib --csl=plos.csl srsdq.md -o srsdq.docx 
			
srsdq_MEE.pdf: srsdq_MEE.md margins.sty srsdq.bib
	pandoc -H margins.sty --bibliography srsdq.bib --csl=methods-in-ecology-and-evolution.csl srsdq_MEE.md -o srsdq.pdf 
	pdftk srsdq.pdf srsdq_figures.pdf output srsdq_MEE.pdf
	evince srsdq_MEE.pdf		
