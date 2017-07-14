SOURCES = main.tex\
		  chapters/* \
		  appendices/* \
		  scraps/* \
		  references.bib

CRAP=idx,ilg,ind,aux,bbl,blg,dvi,lof,lot,log,toc,out,xml,bcf

# &> /dev/null silences the noisey latex build process
# remove these to get your errors and what not back
all: main.pdf

debug:
	pdflatex main.tex
	biber main.bcf
	pdflatex main.tex
	pdflatex main.tex
	rm -rf *.{$(CRAP)}

main.pdf: $(SOURCES)
	pdflatex main.tex &> /dev/null
	biber main.bcf &> /dev/null
	pdflatex main.tex &> /dev/null
	pdflatex main.tex &>/dev/null
	rm -f *.{$(CRAP)}



.PHONY: clean, debug, all

clean:
	rm -f *.{$(CRAP)}
	rm -f main.pdf
