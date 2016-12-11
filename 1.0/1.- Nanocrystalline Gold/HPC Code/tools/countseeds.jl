# countseeds.jl
#
# Usage:
#   julia countseeds.jl <datafile> <outputfile>
#
# From a HDF5 Simulation record <datafile> read size distribution of seeds
# are append them to <outputfile>

using HDF5

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

main(ARGS[1], ARGS[2])
