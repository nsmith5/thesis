SOURCES=main.tex chapters/dft_of_freezing.tex chapters/ddft.tex chapters/fundamentals.tex references.bib

main.pdf: $(SOURCES)
	latex main.tex
	bibtex main
	makeindex -s mcgillpage.ist abbr
	makeindex -s mcgillpage.ist keylist
	latex main.tex
	pdflatex main.tex
	rm *.{idx,ilg,ind,aux,bbl,blg,dvi,lof,lot,log,toc}
	
.PHONY: clean

clean:
	rm -f *.{idx,ilg,ind,aux,bbl,blg,dvi,lof,lot,log,toc}
	rm main.pdf
