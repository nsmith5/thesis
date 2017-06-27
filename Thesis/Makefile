SOURCES = main.tex chapters/dft_of_freezing.tex chapters/noise.tex\
		  chapters/fundamentals.tex chapters/binary_pfc.tex\
		  chapters/introduction.tex\
		  references.bib 

CRAP=idx,ilg,ind,aux,bbl,blg,dvi,lof,lot,log,toc,out

# &> /dev/null silences the noisey latex build process
# remove these to get your errors and what not back
all: main.pdf

debug:
	latex main.tex
	bibtex main
	latex main.tex
	pdflatex main.tex
	rm -rf *.{$(CRAP)}

main.pdf: $(SOURCES)
	latex main.tex &> /dev/null
	bibtex main &> /dev/null
	latex main.tex &> /dev/null
	pdflatex main.tex &>/dev/null
	rm -f *.{$(CRAP)}



.PHONY: clean, debug, all

clean:
	rm -f *.{$(CRAP)}
	rm -f main.pdf
