#!/usr/bin/env julia
# InfCell.jl

using Random
using ProgressMeter
using CSV
using DataFrames
using Suppressor

@suppress include("Constants.jl")
using .Constants
c = Constants

@suppress include("PhysicsFunctions.jl")
using .PhysicsFunctions
@suppress include("RunningStatistics.jl")
using .RunningStatistics
@suppress include("Common.jl")
using .Common

set_zero_subnormals(true)

function energy_deposition!(particle::Particle, dist::Float64, temp::Float64)::NTuple{2, Float64}
    local opacity::Float64 = sigma_a(particle.opacity, temp)
    local delta_energy::Float64 = @fastmath particle.weight * (1.0 - exp(-opacity * dist))
    @fastmath particle.weight -= delta_energy

    local result_dep::Float64 = delta_energy
    local result_em::Float64 = @fastmath opacity * c.arad * temp^4 * c.vol * dist

    return (result_dep, result_em)
end

function deposit_energy!(particle::Particle, dep_1::Float64, dep_2::Float64, em_1::Float64, em_2::Float64, distance::Float64, temp1::Float64, temp2::Float64)::NTuple{4, Float64}
    local (result_dep::Float64, result_em::Float64)
    if @fastmath (particle.mat_num == 1)
        (result_dep, result_em) = energy_deposition!(particle, distance, temp1)
        @fastmath dep_1 += result_dep
        @fastmath em_1 += result_em
    else
        (result_dep, result_em) = energy_deposition!(particle, distance, temp2)
        @fastmath dep_2 += result_dep
        @fastmath em_2 += result_em
    end
    @fastmath particle.t_remaining -= distance / c.sol

    return (dep_1, dep_2, em_1, em_2)
end

function main()::Nothing
    # Mersenne Twister random number generator
    local generator::MersenneTwister = MersenneTwister(1234)

    # Energy balance variables
    local energy_balance::RunningStat = RunningStat()
    local tot_eng::Float64 = @fastmath c.init_intensity * (c.vol / c.sol) + (c.init_temp * c.dens_1 * c_v(c.spec_heat_1, c.init_temp) * c.vol * c.volfrac_1) + (c.init_temp * c.dens_2 * c_v(c.spec_heat_2, c.init_temp) * c.vol * c.volfrac_2)
    local prev_eng::Float64 = 0.0

    # All initial particles have the same weight in a homogeneous cell
    local init_weight::Float64 = @fastmath c.init_intensity * (c.vol / c.sol) / convert(Float64, c.num_particles)
    local particles::Vector{Particle} = [
        Particle(init_weight, isotropic_uvw(generator)..., generate_position(generator)..., generate_material(generator)..., c.delta_t)
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
    local init_weight_1::Float64 = @fastmath c.init_intensity * (c.vol / c.sol) / convert(Float64, num_1)
    local init_weight_2::Float64 = @fastmath c.init_intensity * (c.vol / c.sol) / convert(Float64, num_2)
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
        local energy_em_1::Float64 = 0.0
        local energy_em_2::Float64 = 0.0
        local particles_m1::Int64 = 0
        local particles_m2::Int64 = 0
        prev_eng = tot_eng
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
                (energy_dep_1, energy_dep_2, energy_em_1, energy_em_2) = deposit_energy!(particle, energy_dep_1, energy_dep_2, energy_em_1, energy_em_2, min_distance, temperature_1[index - 1], temperature_2[index - 1])
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
            if @fastmath (particle.mat_num == 1)
                @fastmath particles_m1 += 1
            else
                @fastmath particles_m2 += 1
            end
        end  # Particle loop

        # Compute material properties based on temperature
        local c_v_1::Float64 = @inbounds c_v(c.spec_heat_1, temperature_1[index - 1])
        local c_v_2::Float64 = @inbounds c_v(c.spec_heat_2, temperature_2[index - 1])

        #@fastmath energy_dep_1 /= c.volfrac_1
        #@fastmath energy_dep_2 /= c.volfrac_2
        #@fastmath energy_em_1 /= c.volfrac_1
        #@fastmath energy_em_2 /= c.volfrac_2

        # Update intensity values based on deposition and emission
        @inbounds @fastmath intensity_1[index] = intensity_1[index - 1] + (energy_em_1 - energy_dep_1) * (c.sol / c.vol)
        @inbounds @fastmath intensity_2[index] = intensity_2[index - 1] + (energy_em_2 - energy_dep_2) * (c.sol / c.vol)

        # Calculate temperature change based on intensity values
        @inbounds temperature_1[index] = @fastmath temperature_1[index - 1] + (energy_dep_1 - energy_em_1) / (c_v_1 * c.dens_1 * c.vol)
        @inbounds temperature_2[index] = @fastmath temperature_2[index - 1] + (energy_dep_2 - energy_em_2) / (c_v_2 * c.dens_2 * c.vol)

        # Track energy balance
        tot_eng = @fastmath (intensity_1[index] * c.vol / c.sol) * c.volfrac_1 + (intensity_2[index] * c.vol / c.sol) * c.volfrac_2 + (temperature_1[index] * c.dens_1 * c_v_1 * c.vol) * c.volfrac_1 + (temperature_2[index] * c.dens_2 * c_v_2 * c.vol) * c.volfrac_2
        push(energy_balance, tot_eng - prev_eng)

        # Renormalize
        local normal_energy_1::Float64 = @inbounds @fastmath intensity_1[index] * (c.sol / c.vol) / convert(Float64, particles_m1)
        local normal_energy_2::Float64 = @inbounds @fastmath intensity_2[index] * (c.sol / c.vol) / convert(Float64, particles_m2)
        @simd for particle in particles
            if @fastmath (particle.mat_num == 1)
                particle.weight = normal_energy_1
            else
                particle.weight = normal_energy_2
            end
        end
    end  # Time loop

    println("Average energy balance: ", mean(energy_balance))
    println("Maximum energy balance: ", greatest(energy_balance))

    local tabular::DataFrame = DataFrame(time=times[2:end], intensity1=intensity_1[2:end], temp1=temperature_1[2:end], intensity2=intensity_2[2:end], temp2=temperature_2[2:end])

    CSV.write("out/infcell_1.csv", tabular)
    println("File written")

    return nothing
end

main()