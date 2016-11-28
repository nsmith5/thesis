all:
	julia h5topng.jl

video: 
	ffmpeg -r 10 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p out.mp4	

.PHONY: clean

clean:
	rm *.png
