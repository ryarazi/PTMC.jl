using PTMC
using Plots
plotlyjs() #use plotlyjs backend

σ = 0.01
direction = [0.,0.,-1.] #source is pointed in the -̂z direction
center = [0.,0.,0.02]
source = GuassianSource(σ, direction, center)

inf_coord = 5000. # large enough to be considered infinite
air_box = Box((-inf_coord,-inf_coord,0.), (inf_coord,inf_coord,0.04)) #z from 0 to 0.04
scatterer_box = Box((-inf_coord,-inf_coord,-inf_coord), (inf_coord,inf_coord,0.)) #z from -inf (practicaly) to 0

air = Solid(air_box, Turbid(); name="air")
scatterer = Solid(scatterer_box, Turbid(μₛ=90, μₐ=10, n=1.5, g=0.); name="scatterer")
scene = Scene([air, scatterer])

#simulate 10000 particles
end_particles = [simulate_particle(source, scene) for i in 1:10000]

plot(end_particles)
plot!(scene)
xlims!(-0.3, 0.3)
ylims!(-0.3, 0.3)
zlims!(-0.05, 0.05)