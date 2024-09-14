using PTMC
using Plots
plotlyjs() #use plotlyjs backend

σ = 0.01
direction = [0.0, 0.0, -1.0] #source is pointed in the -̂z direction
center = [0.0, 0.0, 0.02]
source = GuassianSource(σ, direction, center)

inf_coord = 5000.0 # large enough to be considered infinite
air_box = Box((-inf_coord, -inf_coord, 0.0), (inf_coord, inf_coord, 0.04)) #z from 0 to 0.04
scatterer_box = Box((-inf_coord, -inf_coord, -inf_coord), (inf_coord, inf_coord, 0.0)) #z from -inf (practicaly) to 0

air = Solid(air_box, Turbid(); name="air")
scatterer = Solid(scatterer_box, Turbid(μₛ=90, μₐ=10, n=1.5, g=0.0); name="scatterer")
scene = Scene([air, scatterer])

#simulate 1000 particles and aggregate their history
aggregator = HistoryAggregator()
for i in 1:1000 simulate_particle(source, scene, aggregator) end

# plot the scene and the history of the first 10 particles
plt = plot(scene)
for i in 1:10
    plot!(aggregator.history[i]; markersize=1.)
end
xlims!(-0.3, 0.3)
ylims!(-0.3, 0.3)
zlims!(-0.05, 0.05)