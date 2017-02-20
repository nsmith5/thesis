# countseeds.jl
#
# Usage:
#   julia countseeds.jl <datafile> <outputfile>
#
# From a HDF5 Simulation record <datafile> read size distribution of seeds
# are append them to <outputfile>

using HDF5

## Proprocessing

type Cluster
    size::Int
    solid::Bool
    function Cluster()
        new(0, false)
    end
end

function grow!(cluster::Cluster)
    cluster.size += 1
end

function solid!(cluster::Cluster)
    cluster.solid = true
end

function make_islands(n, c; threshold_n = 0.50, threshold_c = 0.60)
    return (n .> threshold_n, c .> threshold_c)
end

function neighbours(site, shape)
    i, j = site
    N, M = shape
    return [
        (mod1(i-1, N), mod1(j-1, M)), 
        (mod1(i-1, N), j),
        (mod1(i-1, N), mod1(j+1, M)), 
        (i, mod1(j-1, M)),
        (i, mod1(j+1, M)),
        (mod1(i+1, N), mod1(j-1, M)), 
        (mod1(i+1, N), j),
        (mod1(i+1, N), mod1(j+1, M))
    ]
end

function depth_first_search(isgold::AbstractArray{Bool, 2}, 
                            issolid::AbstractArray{Bool, 2}, 
                            isvisited::AbstractArray{Bool, 2}, 
                            cluster::Cluster,
                            i::Int, 
                            j::Int)
    stack = Tuple{Int, Int}[]
    push!(stack, (i, j))
    while !isempty(stack)
        site = pop!(stack)
        if !isvisited[site...] && isgold[site...]
            isvisited[site...] = true
            grow!(cluster)
            if issolid[site...]
                solid!(cluster)
            end
            for neighbour in neighbours(site, size(isvisited))
                push!(stack, neighbour)
            end
        end
    end
end

function find_clusters(n, c)
    clusters = Cluster[]
    N, M = size(n)
    isvisited = zeros(Bool, N, M)
    issolid, isgold = make_islands(n, c)
    for j in 1:M
        for i in 1:N
            if !isvisited[i,j] && isgold[i,j]
                cluster = Cluster()
                depth_first_search(isgold, issolid, isvisited, cluster, i, j)
                push!(clusters, cluster)
            end
        end
    end
    return clusters
end

function main(h5_filename, output_filename)
    h5file = h5open(h5_filename)
    for group in names(h5file)
        n = read(group["Density"])
        c = read(group["Concentration"])
        clusters = find_clusters(n, c)
        for cluster in clusters
            open(output_filename, "a+") do f
                write(f, "$group $cluster.size $Int(cluster.solid)")
            end
        end
    end
    close(h5file)
end

main(ARGS[1], ARGS[2])
