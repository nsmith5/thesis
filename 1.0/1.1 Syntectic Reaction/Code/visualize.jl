using PyPlot

export plot

function plot(s::State)
	imshow(s.c, cmap="coolwarm")
	clim(0, 1)
	colorbar()
	imshow(s.n, cmap="gray", alpha=0.7)
 	clim(0, 1)
	axis("off")
	return
end
