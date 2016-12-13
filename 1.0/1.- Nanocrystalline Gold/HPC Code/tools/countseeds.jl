# countseeds.jl
#
# Usage:
#   julia countseeds.jl <datafile> <outputfile>
#
# From a HDF5 Simulation record <datafile> read size distribution of seeds
# are append them to <outputfile>

using HDF5

function seedstatus(i, j, s)
    # Check if this index is in a seed
    # return 0 if no or 'id' of seed otherwise

    N, M = size(s)

    # Check all neighbour for seeds, the hashed positions below
    # show the places we need to check relative to our current co-ordinates
    #     j-1  j  j+1
    # i-1  # | # | O
    #     ----------
    #  i   # | x | O
    #     ----------
    # i+1  # | O | O

    seedid = s[i,j] | s[mod1(i-1, N), mod1(j-1, M)] |
                      s[i, mod1(j-1, M)] |
                      s[mod1(i+1, N), mod1(j-1, M)] |
                      s[mod1(i-1, N), j]
    return seedid

end

function firstseed(c, threshold)
    # Find the indices of the first seed we can find
    i = 1
    while c[i] < threshold
        i += 1
    end
    return ind2sub(c, i)
end

function findseeds(n, c, thresholdc=0.8, thresholdn=1.0)
    phases = String[]
    areas = Int64[]
    seedarr = zeros(Int64, size(c))
    N, M = size(c)

    i, j = firstseed(c, thresholdc)
    println("First seed at ($i, $j)")
    seedarr[i,j] = 1
    push!(areas, 1)
    n[i,j] > thresholdn ? push!(phases, "Solid") : push!(phases, "Liquid")

    for jj in j:M
        for ii in 1:N
            seedid = seedstatus(ii, jj, seedarr)
            inseed = seedid > 0
            if c[ii,jj] > thresholdc && inseed
                seedarr[ii,jj] = seedid
                try
                    areas[seedid] += 1
                catch
                    println("(ii, jj) = ($ii, $jj)")
                    println("areas = $areas")
                    println("phases = $phases")
                    return seedarr
                end
                if (phases[seedid] == "Liquid" && n[ii,jj] > thresholdn)
                    phases[seedid] = "Solid"
                end
            elseif c[ii,jj] > thresholdc && !inseed
                seedarr[ii,jj] = length(phases) + 1
                push!(areas, 1)
                n[ii,jj] > thresholdn ? push!(phases, "Solid") : push!(phases, "Liquid")
            end
        end
    end

    return seedarr, areas, phases
end

function get_pool_sizes(filename)
    fid = h5open(filename)
    t = Float64[]
    R = Float64[]
    for time in names(fid)
        group = fid[time]
        push!(t, read(attrs(group)["Time"])[1])
        Δx = read(attrs(group)["dx"])[1]
        c = read(group["Concentration"])
        push!(R, pool_size(c, Δx))
    end
    return t, R
end

function pool_size(c, Δx)
    # Read PHYSICAL REVIEW B 69, 081201(R) (2004)
    # for an explanation of this algorithm
    ising  = c |> φ -> φ - mean(c) |> sign
    counter = 0
    N, N = size(c)
    for i in 1:N
        for j in 1:N
            ising[i,j] ≠ ising[mod1(i+1,N), j] ?
                counter += 1 : nothing
            ising[i,j] ≠ ising[i, mod1(j+1,N)] ?
                counter += 1 : nothing
        end
    end
    L = counter*Δx
    A = (N*Δx)^2
    return A/L
end

function main(infile, outfile)
    t, R = get_pool_sizes(infile)
    writedlm(outfile, [t R])
    return
end

#main(ARGS[1], ARGS[2])
