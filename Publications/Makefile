TARGET = main

CRAP =bcf,xml,aux,log,out,bbl,blg

$(TARGET).pdf: $(TARGET).tex references.bib
	pdflatex $(TARGET).tex
	bibtex $(TARGET)
	pdflatex $(TARGET).tex
	pdflatex $(TARGET).tex
	rm -f *.{$(CRAP)} mainNotes.bib

.PHONY: clean

clean:
	rm -f main.pdf
	rm -f *.{$(CRAP)} mainNotes.bib
