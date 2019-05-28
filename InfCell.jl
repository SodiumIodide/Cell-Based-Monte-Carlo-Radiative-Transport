#!/usr/bin/env julia
# InfCell.jl

include("PhysicsFunctions.jl")
using .PhysicsFunctions
using Random
using ProgressMeter
using CSV
using DataFrames

include("Constants.jl")
using .Constants
c = Constants

set_zero_subnormals(true)

mutable struct Particle
    # Weight value
    weight::Float64

    # Angle data
    u_val::Float64
    v_val::Float64
    w_val::Float64

    # Location data
    x_loc::Float64
    y_loc::Float64
    z_loc::Float64

    # Material properties
    mat_num::Int64
    opacity::Float64
    spec_heat::Float64
    dens::Float64
    chord::Float64

    # Time data
    t_remaining::Float64
end

@inline function isotropic_u(rand_num::Float64)::Float64
    # Random number between 0 and 1
    return @fastmath 2.0 * rand_num - 1.0
end

@inline function isotropic_v(rand_num::Float64, iso_u::Float64)::Float64
    # Random number between 0 and 1
    return @fastmath sqrt(1.0 - iso_u^2) * cos(2.0 * pi * rand_num)
end

@inline function isotropic_w(rand_num::Float64, iso_u::Float64)::Float64
    # Random number between 0 and 1
    return @fastmath sqrt(1.0 - iso_u^2) * sin(2.0 * pi * rand_num)
end

@inline function isotropic_uvw(gen::MersenneTwister)::NTuple{3, Float64}
    local u_value::Float64 = @fastmath isotropic_u(rand(gen))
    return @fastmath (u_value, isotropic_v(rand(gen), u_value), isotropic_w(rand(gen), u_value))
end

@inline function get_material_properties(mat_num::Int64)::NTuple{4, Float64}
    if (mat_num == 1)
        return (c.opacity_1, c.spec_heat_1, c.dens_1, c.chord_1)
    else
        return (c.opacity_2, c.spec_heat_2, c.dens_2, c.chord_2)
    end
end

@inline function generate_material(gen::MersenneTwister)::Tuple{Int32, Float64, Float64, Float64, Float64}
    if @fastmath (rand(gen) < c.volfrac_1)
        return (1, get_material_properties(1)...)
    else
        return (2, get_material_properties(2)...)
    end
end

function switch_material(particle::Particle)::Nothing
    if (particle.mat_num == 1)
        particle.mat_num = 2
        particle.opacity = c.opacity_2
        particle.spec_heat = c.spec_heat_2
        particle.dens = c.dens_2
        particle.chord = c.chord_2
    else
        particle.mat_num = 1
        particle.opacity = c.opacity_1
        particle.spec_heat = c.spec_heat_1
        particle.dens = c.dens_1
        particle.chord = c.chord_1
    end

    return nothing
end

@enum CollisionType begin
    x_bound = 1
    y_bound = 2
    z_bound = 3
end

function dist_to_boundary(particle::Particle)::Tuple{Float64, CollisionType}
    local (x_dist::Float64, y_dist::Float64, z_dist::Float64)
    if @fastmath (particle.u_val > 0.0)
        x_dist = @fastmath (c.x_len - particle.x_loc) / particle.u_val
    else
        x_dist = @fastmath abs(particle.x_loc / particle.u_val)
    end
    if @fastmath (particle.v_val > 0.0)
        y_dist = @fastmath (c.y_len - particle.y_loc) / particle.v_val
    else
        y_dist = @fastmath abs(particle.y_loc / particle.v_val)
    end
    if @fastmath (particle.w_val > 0.0)
        z_dist = @fastmath (c.z_len - particle.z_loc) / particle.w_val
    else
        z_dist = @fastmath abs(particle.z_loc / particle.w_val)
    end

    local min_col_dist::Float64 = @fastmath minimum([x_dist, y_dist, z_dist])
    if @fastmath (min_col_dist == x_dist)
        return (x_dist, x_bound)
    elseif @fastmath (min_col_dist == y_dist)
        return (y_dist, y_bound)
    else
        return (z_dist, z_bound)
    end
end

@inline function dist_to_transition(gen::MersenneTwister, chord::Float64)::Float64
    return @fastmath -(chord * c.sol) * log(rand(gen))
end

@inline function dist_to_census(particle::Particle)::Float64
    return @fastmath particle.t_remaining * c.sol
end

function move_particle(particle::Particle, distance::Float64)::Nothing
    @fastmath particle.x_loc += distance * particle.u_val
    @fastmath particle.y_loc += distance * particle.v_val
    @fastmath particle.z_loc += distance * particle.w_val

    return nothing
end

function energy_deposition(particle::Particle, dist::Float64, temp::Float64)::Float64
    local opacity::Float64 = sigma_a(particle.opacity, temp)
    local delta_intensity::Float64 = @fastmath particle.weight * (1.0 - exp(-opacity * dist))
    particle.weight = @fastmath particle.weight * exp(-opacity * dist)
    return delta_intensity * opacity
end

function deposit_energy(particle::Particle, accum_1::Float64, accum_2::Float64, distance::Float64, temp1::Float64, temp2::Float64)::NTuple{2, Float64}
    if @fastmath (particle.mat_num == 1)
        @fastmath accum_1 += energy_deposition(particle, distance, temp1)
    else
        @fastmath accum_2 += energy_deposition(particle, distance, temp2)
    end

    return (accum_1, accum_2)
end

function energy_kill(particle::Particle, temp::Float64)::Float64
    local opacity::Float64 = sigma_a(particle.opacity, temp)
    local delta_intensity::Float64 = particle.weight
    particle.weight = 0.0
    return delta_intensity * opacity
end

function boundary_hit(particle::Particle, bound::CollisionType)::Nothing
    # Change direction of travel (all surfaces reflective)
    if (bound == x_bound)
        particle.u_val = -particle.u_val
    elseif (bound == y_bound)
        particle.v_val = -particle.v_val
    else
        particle.w_val = -particle.w_val
    end

    return nothing
end

function main()::Nothing
    # Mersenne Twister random number generator
    local generator::MersenneTwister = MersenneTwister(1234)

    local source_location::NTuple{3, Float64} = (c.src_x, c.src_y, c.src_z)

    # All initial particles have the same weight in a homogeneous cell
    local init_weight::Float64 = @fastmath c.init_intensity / convert(Float64, c.num_particles)
    local particles::Vector{Particle} = [
        Particle(init_weight, isotropic_uvw(generator)..., source_location..., generate_material(generator)..., c.delta_t)
        for i=1:c.num_particles
    ]

    # Renormalize weights
    local num_1::Int64 = 0
    local num_2::Int64 = 0
    @simd for particle in particles
        if @fastmath (particle.mat_num == 1)
            @fastmath num_1 += 1
        else
            @fastmath num_2 += 1
        end
    end
    local init_weight_1::Float64 = @fastmath c.init_intensity / convert(Float64, num_1)
    local init_weight_2::Float64 = @fastmath c.init_intensity / convert(Float64, num_2)
    @simd for particle in particles
        if @fastmath (particle.mat_num == 1)
            particle.weight = init_weight_1
        else
            particle.weight = init_weight_2
        end
    end

    local times::Vector{Float64} = @fastmath [(c.t_init + c.delta_t * i) * c.sol for i=0:c.num_t]
    local intensity_1::Vector{Float64} = @fastmath fill(0.0, c.num_t + 1)
    local temperature_1::Vector{Float64} = @fastmath fill(0.0, c.num_t + 1)
    local intensity_2::Vector{Float64} = @fastmath fill(0.0, c.num_t + 1)
    local temperature_2::Vector{Float64} = @fastmath fill(0.0, c.num_t + 1)
    @inbounds intensity_1[1] = intensity_2[1] = c.init_intensity
    @inbounds temperature_1[1] = temperature_2[1] = c.init_temp

    # Outer loop - times
    @showprogress 1 for index=2:@fastmath(c.num_t + 1)
        local energy_dep_1::Float64 = 0.0
        local energy_dep_2::Float64 = 0.0
        # Inner loop - particles
        @simd for particle in particles
            particle.t_remaining = c.delta_t  # s
            # Per particle loop - propagate data in time-step
            while @fastmath (particle.t_remaining > 0.0)
                # Sample distances in time-step
                # Distance to boundary
                local (dist_b::Float64, bound::CollisionType) = dist_to_boundary(particle)  # cm
                # Distance to transition
                local dist_t::Float64 = dist_to_transition(generator, particle.chord)  # cm
                # Distance to census (end of time-step)
                local dist_c::Float64 = dist_to_census(particle)  # cm
                # Apply distance computations
                local min_distance::Float64 = @fastmath minimum([dist_b dist_c dist_t])  # cm
                (energy_dep_1, energy_dep_2) = deposit_energy(particle, energy_dep_1, energy_dep_2, min_distance, temperature_1[index - 1], temperature_2[index - 1])
                move_particle(particle, min_distance)
                # Particle meets boundary
                if @fastmath(min_distance == dist_b)
                    boundary_hit(particle, bound)
                    @fastmath particle.t_remaining -= (min_distance / c.sol)  # s
                # Material changes
                elseif @fastmath (min_distance == dist_t)
                    switch_material(particle)
                    @fastmath particle.t_remaining -= (min_distance / c.sol)  # s
                else
                    particle.t_remaining = 0.0  # s
                end
            end  # Particle "alive" loop
        end  # Particle loop

        # Kill particles if they fall below weight cutoff
        local particles_killed::Bool = false
        local particles_m1::Int64 = 0
        local particles_m2::Int64 = 0
        @simd for particle in particles
            # Check particle weight for small values and kill particles
            if @fastmath (particle.mat_num == 1)
                @fastmath particles_m1 += 1
                if @fastmath (particle.weight < c.weight_cutoff)
                    @fastmath energy_dep_1 += energy_kill(particle, temperature_1[index - 1])
                    particles_killed = true
                end
            else
                @fastmath particles_m2 += 1
                if @fastmath (particle.weight < c.weight_cutoff)
                    @fastmath energy_dep_2 += energy_kill(particle, temperature_2[index - 1])
                    particles_killed = true
                end
            end
        end

        # Compute material properties based on temperature
        local sigma_a_1::Float64 = @inbounds sigma_a(c.opacity_1, temperature_1[index - 1])
        local sigma_a_2::Float64 = @inbounds sigma_a(c.opacity_2, temperature_2[index - 1])
        local c_v_1::Float64 = @inbounds c_v(c.spec_heat_1, temperature_1[index - 1])
        local c_v_2::Float64 = @inbounds c_v(c.spec_heat_2, temperature_2[index - 1])

        @fastmath energy_dep_1 /= c.volfrac_1
        @fastmath energy_dep_2 /= c.volfrac_2

        # Calculate emission terms
        local energy_em_1::Float64 = @inbounds @fastmath c.sol * c.arad * sigma_a_1 * temperature_1[index - 1]^4
        local energy_em_2::Float64 = @inbounds @fastmath c.sol * c.arad * sigma_a_2 * temperature_2[index - 1]^4

        # Update intensity values based on deposition and emission
        @inbounds @fastmath intensity_1[index] = intensity_1[index - 1] + (energy_em_1 - energy_dep_1) * c.sol * c.delta_t
        @inbounds @fastmath intensity_2[index] = intensity_2[index - 1] + (energy_em_2 - energy_dep_2) * c.sol * c.delta_t

        # Calculate temperature change based on intensity values
        @inbounds temperature_1[index] = @fastmath temperature_1[index - 1] + (energy_dep_1 - energy_em_1) / (c_v_1 * c.dens_1) * c.delta_t
        @inbounds temperature_2[index] = @fastmath temperature_2[index - 1] + (energy_dep_2 - energy_em_2) / (c_v_2 * c.dens_2) * c.delta_t

        # Renormalize
        if @fastmath (particles_killed)
            local normal_intensity_1::Float64 = @inbounds @fastmath intensity_1[index] / convert(Float64, particles_m1)
            local normal_intensity_2::Float64 = @inbounds @fastmath intensity_2[index] / convert(Float64, particles_m2)
            @simd for particle in particles
                if @fastmath (particle.mat_num == 1)
                    particle.weight = normal_intensity_1
                else
                    particle.weight = normal_intensity_2
                end
            end
        end
    end  # Time loop

    local tabular::DataFrame = DataFrame(time=times[2:end], intensity1=intensity_1[2:end], temp1=temperature_1[2:end], intensity2=intensity_2[2:end], temp2=temperature_2[2:end])

    CSV.write("out/infcell.csv", tabular)
    println("File written")

    return nothing
end

main()
