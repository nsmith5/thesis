println("Compiling Modules...")
using HDF5
using PyPlot

filename = "Data.h5"

function make_pictures(filename)
	println("Loading hdf5 datafile...")
	fid = h5open(filename)

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

make_pictures(filename)
