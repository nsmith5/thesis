using IterTools
using PyPlot

function load(filename)
    raw = readdlm(filename)
    t = raw[:, 1]
    phase = raw[:, 3]
    return t, phase
end

function analysis(times, phases)
    output_times = Float64[]
    solid_fraction = Float64[]
   
    sorted = sort(collect(zip(times, phases)), by = pair -> pair[1])
    
    for droplets in groupby(pair -> pair[1], sorted)
        push!(output_times, droplets[1][1])
        push!(solid_fraction, mean(pair -> pair[2], droplets))
    end

    liquid_fraction = - (solid_fraction - 0.5) + 0.5 

    return output_times, liquid_fraction
end

function make_figure(time, fraction)
    loglog(time, fraction)
    xlim(minimum(time), 550_000)
    xlabel(L"Time ($\Delta t$)")
    ylabel(L"Fraction of Uncrystallized droplets $N_{liq} / N_{total}$")
    savefig("loglogincubation.svg")
end

function main(args)
    filename = args[1]
    times, phases = load(filename)
    time, fraction = analysis(times, phases)
    make_figure(time, fraction)
    return 0
end

main(ARGS)
