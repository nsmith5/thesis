function island(label, islands)
    island = 0
    for key ∈ keys(islands)
        label ∈ islands[key] ? island = key : nothing
    end
    island == 0 ? error("this label does not have an island") : nothing
    return island
end

function merge_islands!(currentlabel, connectlabel, islands)
    currentisland = island(currentlabel, islands)
    connectisland = island(connectlabel, islands)
    if currentisland != connectisland
        push!(islands[connectisland], islands[currentisland]...)
        delete!(islands, currentisland)
    end
end

function update_label!(label, solid)
    label[1] += 1                   # Increase count by 1
    label[2] = label[2] || solid    # Update the phase
    return nothing
end

function checkconnections(label_arr, i, j)
    N, M = size(label_arr)
    neighbours = [
        (mod1(i-1, N), mod1(j-1, M)),
        (mod1(i+1, N), mod1(j-1, M))
    ]
    return (label_arr[neighbours[1]...], label_arr[neighbours[2]...])
end

function resolve_connections!(neighbours, islands, currentlabel)
    #currentisland = island(currentlabel, islands)
    for neighbour in neighbours
        currentisland = island(currentlabel, islands)
        if neighbour == 0
            continue
        elseif neighbour ∉ islands[currentisland]
            merge_islands!(currentlabel, neighbour, islands)
        end
    end
    return nothing
end

function aggregate(island, labels)
    total = Any[0, false]
    for label in island
        total[1] += labels[label][1]
        total[2] = total[2] || labels[label][2]
    end
    return total
end

function aggregateall(islands, labels)
    out = Any[]
    for key in keys(islands)
        island = islands[key]
        push!(out, aggregate(island, labels))
    end
    return out
end

function find_droplets(n, c, thresholdc = 0.8, thresholdn = 1.0)
    N, M = size(c)
    label_arr = zeros(Int64, size(n))
    islands = Dict{Int64, Vector{Int64}}()
    labels = []
    currentlabel = 0
    currentisland = 0
    inisland = false

    for j ∈ 1:M
        for i ∈ 1:N
            solid = n[i,j] > thresholdn
            gold = c[i,j] > thresholdc

            if gold && !inisland

                # Flip the island switch on
                inisland = true

                # New label and new island
                currentisland += 1
                currentlabel += 1
                push!(labels, Any[0, false])
                islands[currentisland] = [currentlabel]
                # Typical update and resolve connections
                label_arr[i,j] = currentlabel
                update_label!(labels[currentlabel], solid)
                neighbours = checkconnections(label_arr, i, j)
                resolve_connections!(neighbours, islands, currentlabel)

            elseif gold && inisland

                # Typical update and resolve connections
                label_arr[i,j] = currentlabel
                update_label!(labels[currentlabel], solid)
                neighbours = checkconnections(label_arr, i, j)
                resolve_connections!(neighbours, islands, currentlabel)

            elseif !gold && inisland

                # Flip the island switch off
                inisland = false

            elseif !gold && !inisland

                # Do nothing
                nothing

            end
        end
    end

    # Aggregate the data for each island
    agg_data = aggregateall(islands, labels)

    solid_islands = Int64[]
    liquid_islands = Int64[]

    for island in agg_data
        if island[2] == false
            push!(liquid_islands, island[1])
        else
            push!(solid_islands, island[1])
        end
    end


    return liquid_islands, solid_islands
end
