using PTMC
import PTMC: Ray, Box, Sphere, Photon
import PTMC: ⟂, sample_source, reflect, unpolarized_reflectance, refract, spin, roulette, rotate_vector, advance, rotation_a2b, find_container, exit_solid
using Test
using StaticArrays
using Statistics
using Random
using LinearAlgebra

@testset "PTMC.jl" begin
    @testset verbose = true "Unit Tests" begin
        @testset "linear rotations" begin    
            u = normalize(@SVector randn(3))
            cosθ = 2*rand() - 1
            ψ = 2pi*rand()
    
            @test norm(u) ≈ norm(rotate_vector(u, cosθ, ψ))
            @test u ≈ -rotate_vector(u, -1, ψ)
            @test dot(rotate_vector(u, cosθ, ψ), u) ≈ cosθ
    
            v = normalize(@SVector randn(3))
            rot = rotation_a2b(u,v)
            @test all(rot*u .≈ v)
            @test all(rot*cross(u,v) .≈ cross(u,v)) #axis of rotation stays in the same place
        end;
    
        @testset "geometric calculations" begin    
            ray = Ray(Tuple(rand(3)), Tuple(normalize(randn(3))))
    
            b = Box((0.,0.,0.), (1.,1.,1.))
    
            @test rand(3)∈b
            @test rand(3).+1∉b
            @test rand(3).-1∉b
    
            p = @SVector rand(3)
            direction = zeros(3)
            direction[rand((1,2,3))] = rand((-1,1))
            ray = Ray(p, direction)
    
            inter = ray ∩ b
            @test !(typeof(inter)<:Tuple)
    
            normal = isapprox.(inter, b.vmin; atol=eps(Float64)) - isapprox.(inter, b.vmax; atol=eps(Float64))
            @test all(normal .≈ -direction)
    
            R = 2. + 3*rand()
            sphere = Sphere(zeros(3), R)
            ray = Ray((@SVector rand(3)), (@SVector randn(3)))
            @test ray.origin∈sphere
            @test !isnothing(ray∩sphere)
            @test typeof(ray∩sphere)==SVector{3,Float64}
            @test norm(ray∩sphere) ≈ R
            
            direction = zeros(3)
            direction[rand((1,2,3))] = rand((-1,1))
            ray_axis = Ray(zeros(3), direction)
            @test ray_axis∩sphere ≈ R .* direction
    
            p_far =  10 .+ (@SVector rand(3))
            ray = Ray(p_far, -p_far)
            @test !isnothing(ray∩sphere)
            @test  length(ray∩sphere)==2
            ray_neg = Ray(p_far, p_far)
            @test isnothing(ray_neg∩sphere)
        end;
    
        @testset "laser sampling" begin    
            σ = rand() + 0.1
            direction = normalize(@SVector randn(3))
            center = @SVector rand(3)
            las = Laser(σ, direction, center) #laser is pointed to xy plane
            photons = [sample_source(las) for _ in 1:100000]
            
            photons_r = stack((p.r for p in photons))
            
            @test isapprox(vec(mean(photons_r, dims=2)), center, atol=5.e-2) #center is ok
            @test isapprox(mean(sum((photons_r .- center).^2, dims=1)), (√2σ)^2, atol=5.e-2) #σ is ok
            
            p1, p2 = rand(photons, 2)
            @test (p1.r - p2.r) ⟂ direction
        end
    
        @testset "scene creation" begin    
            mat1 = Turbid(μₛ=0.5, μₐ=0., n=1., g=0.)
    
            start_point = rand(3)
    
            s1 = Solid(Box(start_point, start_point + ones(3)), mat1)
    
            direction = 2*bitrand(3) .- 1.
            mat2 = Turbid(μₛ=1.5, μₐ=0., n=1., g=0.)
            s2 = Solid(Box(start_point+direction, start_point+direction+ones(3)), mat2)
    
            s3 = Solid(Box(min.(minimum(s1.geometry), minimum(s2.geometry))
                          ,max.(maximum(s1.geometry), maximum(s2.geometry))), mat1)
    
            scene = Scene([s1, s2, s3])
    
            p1 = Photon((start_point + 0.5*ones(3)), direction+1e-6*rand(3))
    
            @test p1∈scene
            @test p1∈s1
            @test p1∈s3
    
            p2 = advance(p1, 1.) 
    
            @test p2∈s2
    
            @test find_container(p1, scene) == s1
    
            rbound, normal = exit_solid(p1, find_container(p1, scene))
            @test rbound≈(p1.r+p2.r)/2 atol=1.e-3
        end;
    
        @testset "particle reflect" begin    
            #wikipedia example https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
            p = Photon([0, 0, 0.], [1/sqrt(2), -1/sqrt(2), 0.])
            n1 = 1.
            n2 = 1. / 0.9
            normal = [0, 1., 0.]
            reflect(p, normal)
    
            @test all(reflect(p, normal).u .≈ [1/sqrt(2), 1/sqrt(2), 0.0])
            @test norm(refract(p, normal, n1, n2).u - [0.636396, -0.771362, 0.0]) < 1.e-5
        end;
    
        @testset "particle fresnel reflectance" begin    
            p = Photon(zeros(3), normalize(randn(3)))
            normal = normalize(@SVector randn(3))
            normal = dot(normal, p.u)<0 ? normal : -normal
    
            @test unpolarized_reflectance(p, normal, 1., 1.) ≈ 0. atol=1.e-10
            @test unpolarized_reflectance(p, normal, 1., 1.e8) ≈ 1. atol=1.e-5
            @test unpolarized_reflectance(p, normal, 1.e8, 1.) ≈ 1. atol=1.e-5
    
            n1 = rand()*3 + 1.
            n2 = rand()*3 + 1.
            @test unpolarized_reflectance(p, -p.u, n1, n2) ≈ (n1-n2)^2/(n1+n2)^2 #https://en.wikipedia.org/wiki/Fresnel_equations#Normal_incidence
        end;
    
        @testset "particle snell law" begin
            p = Photon(zeros(3), normalize(randn(3)))
    
            normal = normalize(@SVector randn(3))
            normal = dot(normal, p.u)<0 ? normal : -normal
    
            n1 = rand()*3 + 1.
            n2 = rand()*3 + 1.
    
            reflected = reflect(p, normal)
            refracted = refract(p, normal, n1, n2)
            R = unpolarized_reflectance(p, normal, n1, n2)
    
            @test dot(p.u, normal) ≈ -dot(reflected.u, normal) # incident angle = reflection angle
    
            hit = [0., 0., 0.]
            start = hit - p.u
            hit2reflect = hit + reflected.u
            if !isnothing(refracted) hit2refract = hit + refracted.u end
    
            normal_vec = hit + normal
        end
    
        @testset "particle spin" begin    
            p = Photon(zeros(3), normalize(randn(3)))
            g = rand()*0.2+0.7
            cosθ_samp = [dot(spin(p, g).u, p.u) for i in 1:100000]
            @test g ≈ mean(cosθ_samp) atol=1.e-2
    
            g = 0.
            cosθ_samp = [dot(spin(p, g).u, p.u) for i in 1:100000]
            @test g ≈ mean(cosθ_samp) atol=1.e-2
        end;
    
        @testset "russian roulette" begin    
            roulette_chance = 0.1
    
            photons = vcat([Photon(zeros(3), normalize(randn(3)), rand()) for i in 1:1000000],
                        [Photon(zeros(3), normalize(randn(3)), 2*roulette_chance*rand()) for i in 1:1000000])
    
            energy = sum(x->x.weight, photons)
    
            roulette_photons = [roulette(p,roulette_chance) for p in photons]
            energy_roulette = sum(x->x.weight, roulette_photons)
    
            @test energy≈energy_roulette rtol=1.e-2 #energy conversion
        end
    end

    Nₚ = 100000 #number of particles
    @testset verbose = true "Integral Tests" begin
        @testset "PTMC benchmark - finite slab" begin
            las = Laser(0.04, [0.,0.,-1.], [0.,0.,0.02]) #laser is pointed in the -̂z direction

            air = Solid(Box((-5000.,-5000.,0.), (5000.,5000.,0.04)), Turbid(); name="air")
            scatterer = Solid(Box((-5000.,-5000.,-0.02), (5000.,5000.,0.)), Turbid(μₛ=90, μₐ=10, n=1., g=0.75); name="scatterer")
            air2 = Solid(Box((-5000.,-5000.,-0.04), (5000.,5000.,-0.02)), Turbid(); name="air2")
            scene = Scene([air, scatterer, air2])
            end_particles = [simulate_particle(las, scene) for i in 1:Nₚ]

            R = 0.
            T = 0.
            for p in end_particles
                if p.r[3] > 0.01 R += p.weight/Nₚ end
                if p.r[3] < -0.03 T += p.weight/Nₚ end
            end

            #see chapter 6 in https://omlc.org/classroom/ece532/class4/MCMan.pdf
            @test R ≈ 0.09739 rtol=5.e-1
            @test T ≈ 0.66096 rtol=5.e-1
        end;

        @testset "PTMC benchmark - infinite half-slab" begin
            las = Laser(0.01, [0.,0.,-1.], [0.,0.,0.02]) #laser is pointed in the -̂z direction

            air = Solid(Box((-5000.,-5000.,0.), (5000.,5000.,0.04)), Turbid(); name="air")
            scatterer = Solid(Box((-5000.,-5000.,-1000.), (5000.,5000.,0.)), Turbid(μₛ=90, μₐ=10, n=1.5, g=0.); name="scatterer")
            scene = Scene([air, scatterer])
            end_particles = [simulate_particle(las, scene) for i in 1:Nₚ]

            R = 0.
            for p in end_particles
                if p.r[3] > 0.01 R += p.weight/Nₚ end
            end

            @test R ≈ 0.26 rtol=5.e-1 #see chapter 6 in https://omlc.org/classroom/ece532/class4/MCMan.pdf
        end
    end
end;
