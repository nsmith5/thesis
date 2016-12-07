main.pdf: chapters/noise.tex chapters/fundamentals.tex
	pdflatex main.tex
	pdflatex main.tex
	rm main.bcf
	rm main.log
	rm main.out
	rm main.aux
	rm main.run.xml
	rm main.toc

.PHONY: clean

clean:
	rm main.pdf
