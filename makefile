pdflatex:
	pdflatex -interaction=nonstopmode article.tex

compile:
	pdflatex article.tex
	bibtex article.aux
	pdflatex article.tex
	pdflatex article.tex

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

