# h5topng.jl
#
# Usage:
# 	julia h5topng.jl <filename>
#
# Convert an HDF5 Simulation record  <filename> into a mountain 
# of png images

using HDF5
using PyPlot

function make_pictures(filename)
	println("Loading hdf5 datafile...")
	try
		fid = h5open(filename)
	catch 
		println("Couldn't open $filename, exiting now")
		return 1
	end

	println("Writing figures...")
	for time in names(fid)
		group = fid[time]
	
		density	= group["Density"] |> read
		concentration = group["Concentration"] |> read
	
		imshow(density, cmap="gray")
		clim(-1, 1)		

		imshow(concentration, cmap="coolwarm", alpha=0.5)	
		clim(0, 1)
		colorbar()	
		
		savefig(time*".png")
		println("Wrote figure $time.png")		

		clf()
	end

	close(fid)

	return nothing
end

make_pictures(ARGS[1])
