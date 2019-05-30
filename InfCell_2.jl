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
    local opacity::Float64 = sigma_a(particle.opacity, temp)
    local delta_energy::Float64 = @fastmath particle.weight * (1.0 - exp(-opacity * dist))
    @fastmath particle.weight -= delta_energy

    # Already Time-Integrated
    return delta_energy
end

function deposit_energy!(particle::Particle, accum::Float64, distance::Float64, temp::Float64)::Float64
    @fastmath accum += energy_deposition!(particle, distance, temp)

    return accum
end

function main()::Nothing
    # Mersenne Twister random number generator
    local generator::MersenneTwister = MersenneTwister(1234)

    # Energy balance variables
    local energy_balance::RunningStat = RunningStat()
    local tot_eng::Float64 = @fastmath c.init_intensity * (c.vol / c.sol) + (c.init_temp * c.dens_1 * c_v(c.spec_heat_1, c.init_temp) * c.vol_1) + (c.init_temp * c.dens_2 * c_v(c.spec_heat_2, c.init_temp) * c.vol_2)
    local prev_eng::Float64 = 0.0

    # All initial particles have the same weight in a homogeneous cell
    local init_weight::Float64 = c.init_intensity * (c.vol / c.sol) / convert(Float64, c.num_particles)
    local particles::Vector{Particle} = [
        Particle(init_weight, isotropic_uvw(generator)..., generate_position(generator)..., generate_material(generator)..., c.delta_t)
        for i=1:c.num_particles
    ]

    local times::Vector{Float64} = @fastmath [(c.t_init + c.delta_t * i) * c.sol for i=0:c.num_t]
    local intensity_1::Vector{Float64} = @fastmath fill(0.0, c.num_t + 1)
    local temperature_1::Vector{Float64} = @fastmath fill(0.0, c.num_t + 1)
    local intensity_2::Vector{Float64} = @fastmath fill(0.0, c.num_t + 1)
    local temperature_2::Vector{Float64} = @fastmath fill(0.0, c.num_t + 1)
    @inbounds intensity_1[1] = intensity_2[1] = c.init_intensity
    @inbounds temperature_1[1] = temperature_2[1] = c.init_temp

    # Outer loop - times
    @showprogress 1 for index=2:@fastmath(c.num_t + 1)
        local material_num::Int64 = random_material(generator)

        # Renormalize
        local normal_energy::Float64 = (convert(Float64, c.num_particles) * c.sol)^(-1)
        if @fastmath (material_num == 1)
            @inbounds @fastmath normal_energy *= intensity_1[index - 1] * c.vol_1
        else
            @inbounds @fastmath normal_energy *= intensity_2[index - 1] * c.vol_2
        end
        @simd for particle in particles
            particle.weight = normal_energy
        end

        local energy_dep::Float64 = 0.0
        prev_eng = tot_eng
        # Inner loop - particles
        @simd for particle in particles
            if (isnan(particle.weight))
                error("NaN encountered")
            end
            if (particle.mat_num != material_num)
                update_material!(particle, material_num)
            end
            particle.t_remaining = c.delta_t  # s
            # Per particle loop - propagate data in time-step
            while @fastmath (particle.t_remaining > 0.0)
                # Sample distances in time-step
                # Distance to boundary
                local (dist_b::Float64, bound::CollisionType) = dist_to_boundary(particle)  # cm
                # Distance to census (end of time-step)
                local dist_c::Float64 = dist_to_census(particle)  # cm
                # Apply distance computations
                local min_distance::Float64 = @fastmath minimum([dist_b dist_c])  # cm
                local temperature_mat::Float64
                if @fastmath (material_num == 1)
                    temperature_mat = temperature_1[index - 1]
                else
                    temperature_mat = temperature_2[index - 1]
                end
                energy_dep = deposit_energy!(particle, energy_dep, min_distance, temperature_mat)
                move_particle!(particle, min_distance)
                # Particle meets boundary
                if @fastmath(min_distance == dist_b)
                    boundary_hit!(particle, bound)
                end
            end  # Particle "alive" loop
        end  # Particle loop

        # Compute material properties based on temperature
        local (sigma_a_mat::Float64, c_v_mat::Float64, energy_em::Float64)
        if (material_num == 1)
            sigma_a_mat = @inbounds sigma_a(c.opacity_1, temperature_1[index - 1])
            c_v_mat = @inbounds c_v(c.spec_heat_1, temperature_1[index - 1])

            # Calculate emission term
            energy_em = @inbounds @fastmath c.sol * c.arad * sigma_a_mat * temperature_1[index - 1]^4 * c.vol_1 * c.delta_t

            # Update intensity value based on deposition and emission
            @inbounds intensity_1[index] = @fastmath intensity_1[index - 1] + (energy_em - energy_dep) * (c.sol / c.vol_1)
            @inbounds intensity_2[index] = @fastmath intensity_2[index - 1]

            # Calculate temperature change based on intensity values
            @inbounds temperature_1[index] = @fastmath temperature_1[index - 1] + (energy_dep - energy_em) / (c_v_mat * c.dens_1 * c.vol_1)
            @inbounds temperature_2[index] = @fastmath temperature_2[index - 1]
        else
            sigma_a_mat = @inbounds sigma_a(c.opacity_2, temperature_2[index - 1])
            c_v_mat = @inbounds c_v(c.spec_heat_2, temperature_2[index - 1])

            # Calculate emission term
            energy_em = @inbounds @fastmath c.sol * c.arad * sigma_a_mat * temperature_2[index - 1]^4 * c.vol_2 * c.delta_t

            # Update intensity value based on deposition and emission
            @inbounds intensity_1[index] = @fastmath intensity_1[index - 1]
            @inbounds intensity_2[index] = @fastmath intensity_2[index - 1] + (energy_em - energy_dep) * (c.sol / c.vol_2)

            # Calculate temperature change based on intensity values
            @inbounds temperature_1[index] = @fastmath temperature_1[index - 1]
            @inbounds temperature_2[index] = @fastmath temperature_2[index - 1] + (energy_dep - energy_em) / (c_v_mat * c.dens_2 * c.vol_2)
        end

        # Track energy balance
        tot_eng = @inbounds @fastmath (intensity_1[index] * c.vol_1 / c.sol) + (intensity_2[index] * c.vol_2 / c.sol) + (temperature_1[index] * c.dens_1 * c_v(c.spec_heat_1, temperature_1[index]) * c.vol_1) + (temperature_2[index] * c.dens_2 * c_v(c.spec_heat_2, temperature_2[index]) * c.vol_2)
        push(energy_balance, abs(tot_eng - prev_eng))
    end  # Time loop

    println("Average energy balance: ", mean(energy_balance))
    println("Maximum energy balance: ", greatest(energy_balance))

    local tabular::DataFrame = DataFrame(time=times[2:end], intensity1=intensity_1[2:end], temp1=energy_from_temp.(temperature_1[2:end]), intensity2=intensity_2[2:end], temp2=energy_from_temp.(temperature_2[2:end]))

    CSV.write("out/infcell_2.csv", tabular)
    println("File written")

    return nothing
end

main()
