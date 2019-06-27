#!/usr/bin/env julia
# InfCell.jl

set_zero_subnormals(true)

using Random
using ProgressMeter
using CSV
using DataFrames

include("Constants.jl")
const c = Constants

include("PhysicsFunctions.jl")
using .PhysicsFunctions
include("RunningStatistics.jl")
using .RunningStatistics
include("Common.jl")
using .Common

function energy_deposition!(particle::Particle, dist::Float64, temp::Float64)::Float64
    local opacity_mat::Float64 = sigma_a(particle.opacity, temp)
    local delta_energy::Float64 = @fastmath particle.weight * (1.0 - exp(-opacity_mat * dist))
    @fastmath particle.weight -= delta_energy

    # Time-Integrated
    return delta_energy
end

@inline function deposit_energy!(particle::Particle, accum::Vector{Float64}, dist::Float64, temp::Vector{Float64})::Nothing
    @inbounds @fastmath accum[particle.mat_num] += energy_deposition!(particle, dist, temp[particle.mat_num])

    return nothing
end

function main()::Nothing
    # Mersenne Twister random number generator
    local generator::MersenneTwister = MersenneTwister(123456)

    # Energy balance variables
    local energy_balance::RunningStat = RunningStat()
    local tot_eng::Float64 = @fastmath c.init_intensity * (c.vol_cell / c.sol) + sum(@.(c.init_temp * c.dens * c_v(c.spec_heat, c.init_temp) * c.vol))
    local prev_eng::Float64 = 0.0

    local init_energy::Float64 = energy_from_temp(c.init_temp)  # erg/cm^3

    # Generate the initial materials to map particles to
    local materials::Vector{Int64} = map_material.(rand(generator, c.num_particles))

    # Normalize weights when generating particles
    local num_1::Int64 = @fastmath count(m -> m == 1, materials)
    local num_mat::Vector{Int64} = @fastmath vec([
        num_1
        c.num_particles - num_1
    ])
    local init_weights::Vector{Float64} = @fastmath @. c.init_intensity * (c.vol / c.sol) / convert(Float64, num_mat)
    local particle_weights::Vector{Float64} = @inbounds [init_weights[n] for n in materials]

    # Vector of particles
    local particles::Vector{Particle} = [build_particle(generator, materials[i], particle_weights[i]) for i=1:c.num_particles]

    local times::Vector{Float64} = @fastmath [(c.t_init + c.delta_t * i) * c.sol for i=0:c.num_t]  # cm
    local energy_diff::Vector{Float64} = @fastmath fill(0.0, c.num_t + 1)
    local intensity::Array{Float64, 2} = @fastmath fill(0.0, 2, c.num_t + 1)  # erg/cm^2-s
    local energy::Array{Float64, 2} = @fastmath fill(0.0, 2, c.num_t + 1)  # eV
    @inbounds intensity[1, 1] = intensity[2, 1] = c.init_intensity  # erg/cm^2-s
    @inbounds energy[1, 1] = energy[2, 1] = init_energy  # eV

    # Outer loop - times
    @showprogress 1 for index=2:@fastmath(c.num_t + 1)
        # Track energy balance
        prev_eng = tot_eng

        # Update new particle times for existing particles
        next_time!.(particles)

        # Compute material properties based on temperature
        local sigma_a_mat::Vector{Float64} = c.opacity
        local c_v_mat::Vector{Float64} = c.spec_heat
        local rho_mat::Vector{Float64} = c.dens

        # Calculate emission terms
        local energy_em::Vector{Float64} = @inbounds @fastmath @. c.sol * sigma_a_mat * energy[:, index - 1] * c.vol * c.delta_t

        # Add more particles via emission
        local materials_add::Vector{Int64} = map_material.(rand(generator, c.num_particles_add))
        local num_1_add::Int64 = @fastmath count(m -> m == 1, materials_add)
        local num_add::Vector{Int64} = @fastmath vec([
            num_1_add
            c.num_particles_add - num_1_add
        ])
        local em_weights::Vector{Float64} = @fastmath @. energy_em / convert(Float64, num_add)
        local weights_add::Vector{Float64} = @inbounds [em_weights[n] for n in materials_add]
        local times_add::Vector{Float64} = produce_time.(rand(generator, c.num_particles_add))
        local particles_add::Vector{Particle} = @inbounds [build_particle(generator, materials_add[i], weights_add[i], times_add[i]) for i=1:c.num_particles_add]

        # Append newly-generated particles to existing vector
        particles = vcat(particles, particles_add)

        # Tally for deposited energy
        local energy_dep::Vector{Float64} = zeros(2)

        # Tally for radiative energy reaching boundary
        local energy_cen::Vector{Float64} = zeros(2)

        # Inner loop - particles
        for particle in particles
            # Error check
            if (isnan(particle.weight))
                error("NaN encountered")
            end

            # Per particle loop - propagate data in time-step
            while @fastmath (particle.t_remaining > 0.0)
                # Sample distances in time-step
                # Distance to boundary
                local (dist_b::Float64, bound::CollisionType) = dist_to_boundary(particle)  # cm
                # Distance to transition
                local dist_t::Float64 = dist_to_transition(generator, particle)  # cm
                # Distance to census (end of time-step)
                local dist_c::Float64 = dist_to_census(particle)  # cm
                # Apply distance computations
                local min_distance::Float64 = @fastmath minimum([dist_b dist_c dist_t])  # cm
                @inbounds @fastmath deposit_energy!(particle, energy_dep, min_distance, energy[:, index - 1])
                move_particle!(particle, min_distance)
                # Particle meets boundary
                if @fastmath (min_distance == dist_b)
                    boundary_hit!(particle, bound)
                # Material changes
                elseif @fastmath (min_distance == dist_t)
                    switch_material!(particle)
                end
            end  # Particle "alive" loop

            # Tally particles in each material at end of time-step
            @inbounds @fastmath energy_cen[particle.mat_num] += particle.weight
        end  # Particle loop

        # Splitting and Russian Roulette with particles
        local split_particles::Vector{Particle} = Particle[]
        local weight_target::Vector{Float64} = @fastmath @. energy_cen / (convert(Float64, c.num_particles) * c.volfrac)
        @simd for particle in particles
            local ratio::Float64 = @inbounds particle.weight / weight_target[particle.mat_num]
            if @fastmath (ratio > 1.0)
                local split::Int64 = trunc(Int64, ratio)
                if @fastmath (split > 1)
                    @simd for n = 1:split
                        push!(split_particles, build_particle(generator, particle.mat_num, particle.weight / convert(Float64, split)))
                    end
                    particle.weight = 0.0
                end
            else
                if @fastmath (ratio < rand(generator))
                    particle.weight = 0.0
                else
                    @fastmath particle.weight /= ratio
                end
            end
        end

        # Filter and append new particles
        filter!(p -> p.weight != 0.0, particles)
        particles = vcat(particles, split_particles)

        # Update intensity values based on energy remaining
        @inbounds intensity[:, index] = @fastmath @. energy_cen * (c.sol / c.vol)

        # Calculate temperature change based on deposition and emission
        @inbounds energy[:, index] = @fastmath @. energy[:, index - 1] + (energy_dep - energy_em) / (rho_mat * c_v_mat * c.vol) * 4.0 * c.arad

        # Update energy balance
        tot_eng = @inbounds @fastmath sum(@.(intensity[:, index] * (c.vol / c.sol))) + sum(@.(temp_from_energy(energy[:, index]) * rho_mat * c_v_mat * c.vol))
        @inbounds energy_diff[index] = @fastmath tot_eng - prev_eng
        @fastmath RunningStatistics.push(energy_balance, abs(tot_eng - prev_eng))
    end  # Time loop

    println("Average energy balance: ", mean(energy_balance))
    println("Maximum energy balance: ", greatest(energy_balance))

    local int_temp::Array{Float64, 2} = @. temp_from_energy(intensity / c.sol)

    local tabular::DataFrame = DataFrame(time=times[2:end], intensity1=intensity[1, 2:end], energy1=energy[1, 2:end], intensity2=intensity[2, 2:end], energy2=energy[2, 2:end], energydiff=energy_diff[2:end])

    CSV.write("out/infcelllinear.csv", tabular)
    println("File written")

    return nothing
end

main()
