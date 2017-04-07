TARGET = main

$(TARGET).pdf: $(TARGET).tex references/references.bib
	pdflatex -shell-escape $(TARGET).tex
	bibtex $(TARGET)
	pdflatex $(TARGET).tex
	pdflatex $(TARGET).tex
	mv $(TARGET).pdf crap.pdf
	mv $(TARGET).tex crap.tex
	rm $(TARGET)*
	mv crap.tex $(TARGET).tex
	mv crap.pdf $(TARGET).pdf

.PHONY: clean

clean:
	rm $(TARGET).pdf
