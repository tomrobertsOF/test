include("./siteanimals.jl")

##
max_tours = 100000;
max_size = 10;
delay = n -> max_tours/(1.5*max_size) * n
s, w, toursize = parallelpegarm(max_tours, max_size, getstartingstate());
show([w[n] / (s[1] - delay(n-1)) for n in 1:length(w)])

##
using Plots
plot(filter(>(0),toursize[trunc(Int,delay(max_size)+1):max_tours]))

## Plotting
using Plots;
gr(size = (500, 400));
using DelimitedFiles;
actuals = readdlm("siteanimalsexact.txt");
function plotsampledata(weights, samples, actuals, delay_function)
    estimate = [weights[n] / (samples[1] - delay_function(n-1)) for n in 1:length(w)]
    Plots.plot(0:(length(w)-1), log.(estimate), layout = 4, subplot = 1, label = "Estimate", title = "Log-plot of estimate vs actual", legend = :bottomright)
    Plots.plot!(0:(length(w)-1), samples, subplot = 2, label = "", title = "Samples for each n", ylims = [0, maximum(s)])
    Plots.plot!([estimate[n+1] / estimate[n] for n = 1:length(estimate)-1], subplot = 3, title = "Ratio between steps", legend = :none, ylims = [0, 4])
    Plots.plot!([actuals[n+1] / actuals[n] for n = 1:length(actuals)-1], subplot = 3)
    Plots.plot!(0:length(actuals)-1, log.(actuals), subplot = 1, label = "Actual")
    minlength = min(length.([actuals, estimate])...)
    Plots.plot!((@. (actuals[1:minlength] - estimate[1:minlength]) / actuals[1:minlength]), subplot = 4, title = "Relative Error", legend = :none)
end
plotsampledata(w, s, actuals, delay)

##
plot(toursize)

##
function drawtree(state::SiteTreeState)
    p = plot(minorgrids=true, ticks=(-100:5:100, -100:5:100), gridalpha=0.8, minorgridalpha=0.5, aspect_ratio = :equal, axis=(false,false))
    scatter!(p, Tuple.(state.occupied), markersize = 14, markercolor = :blue,label="occupied")
    scatter!(p, Tuple.(state.growth_candidates), markersize = 10, markercolor = :green, label="a+")
    scatter!(p, Tuple.(state.shrink_candidates), markersize = 8, markercolor = :red, label="a-")
end

## Stats gathering
global const PRINT_PROGRESS=false
num_runs = 1000;
estimates = fill([], num_runs)

max_tours = 10000;
max_size = 40;
delay = @inline n -> max_tours/(2*max_size) * n

Threads.@threads for i in 1:num_runs
    println("Runs complete: $i")
    s, w, toursize = progressivegarm(max_tours, max_size, getstartingstate());
    estimates[i] = [w[n] / (s[1] - delay(n-1)) for n in 1:length(w)]    
end

data = [[estimates[n][m] for n in 1:num_runs] for m in 1:length(estimates[1])]
using Statistics
var.(data)

## Equilibration
function toursizevariance(tours)
    s, w, toursize = parallelpegarm(tours, max_size, getstartingstate());
    return var(toursize)
end

toursizes = Vector{Int}[];
for t in 210000:10000:400000
    global toursizes;
    s, w, ws, ts = pegarm(t, max_size, getstartingstate());
    push!(toursizes, ts)
end
