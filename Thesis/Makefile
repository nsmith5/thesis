main.pdf: main.tex chapters/dft.tex chapters/ddft.tex chapters/statmech.tex references.bib
	pdflatex main.tex
	bibtex main1-blx
	bibtex main2-blx
	pdflatex main.tex
	pdflatex main.tex
	mv main.pdf crap.pdf
	mv main.tex crap.tex
	rm main*
	mv crap.tex main.tex
	mv crap.pdf main.pdf

.PHONY: clean

clean:
	mv main.pdf crap.pdf
	mv main.tex crap.tex
	rm main*
	mv crap.tex main.tex
	mv crap.pdf main.pdf	
