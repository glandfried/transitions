pdflatex:
	pdflatex -interaction=nonstopmode article.tex

compile:
	pdflatex article.tex
	bibtex article.aux
	pdflatex article.tex
	pdflatex article.tex

summary:
	pdflatex summary.tex
	bibtex summary.aux
	pdflatex summary.tex
	pdflatex summary.tex

publish: clean pdfclean landfried-transiciones.pdf landfried-transitions.pdf 
	
landfried-transitions.pdf:
	sed -i 's/\\estrue/\\entrue/g' article.tex
	pdflatex article.tex
	bibtex article.aux
	pdflatex article.tex
	pdflatex article.tex
	cp article.pdf landfried-transitions.pdf
	
landfried-transiciones.pdf:
	sed -i 's/\\entrue/\\estrue/g' article.tex
	pdflatex article.tex
	bibtex article.aux
	pdflatex article.tex
	pdflatex article.tex
	mv article.pdf landfried-transiciones.pdf


pdfclean:
	- rm -f article*.pdf

clean:
	- rm -f *.log
	- rm -f *.soc
	- rm -f *.toc
	- rm -f *.aux
	- rm -f *.out
	- rm -f article.idx
	- rm -f *.bbl
	- rm -f *.bbg
	- rm -f *.dvi
	- rm -f *.blg
	- rm -f *.lof
	- rm -f *.nav
	- rm -f *.snm
	- rm -f *~

