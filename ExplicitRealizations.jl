#!/usr/bin/env julia

include("Constants.jl")
const c = Constants

include("GeometryGen.jl")
include("RunningStatistics.jl")
include("PhysicsFunctions.jl")
using .PhysicsFunctions
using Random
using LinearAlgebra
using DataFrames
using CSV
using ProgressMeter

function main()::Nothing
    set_zero_subnormals(true)

    # Iteration condition
    local generator::MersenneTwister = MersenneTwister(1234)

    # Computational values
    local times::Vector{Float64} = @fastmath [(x * c.delta_t + c.t_init) * c.sol for x in 1:c.num_t]

    # Arrays for online computation
    local stat_1_intensity::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:c.num_t]
    local stat_2_intensity::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:c.num_t]
    local stat_1_temp::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:c.num_t]
    local stat_2_temp::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:c.num_t]

    # Outer loop
    @showprogress 1 for iteration_number = 1:c.max_iterations
        # Initial values in the problem
        local intensity_value::Float64 = c.init_intensity  # erg/cm^2-s
        local temp_value::Vector{Float64} = vec([c.init_temp c.init_temp])  # eV
        local material_num::Int32 = @inbounds @fastmath (rand(generator) <= c.volfrac[1]) ? 1 : 2
        for time_index in 1:c.num_t
            local (t_delta::Vector{Float64}, t_arr::Vector{Float64}, materials::Vector{Int32}, num_cells::Int64) = GeometryGen.get_geometry(c.chord[1], c.chord[2], c.delta_t, c.num_divs, rng=generator, starting_material=material_num)

            local time_in::Vector{Float64} = zeros(2)
            local intensity_deposition::Vector{Float64} = zeros(2)

            local last_material_num::Int32 = material_num

            local opacity_mat::Float64 = @inbounds sigma_a(c.opacity[material_num], temp_value[material_num])

            for (t_delta_index, t_delta_value) in enumerate(t_delta)
                material_num = materials[t_delta_index]
                if @fastmath (material_num != last_material_num)
                    opacity_mat = @inbounds sigma_a(c.opacity[material_num], temp_value[material_num])
                end

                # Update intensity value
                @inbounds @fastmath intensity_deposition[material_num] += t_delta_value * opacity_mat * intensity_value

                @fastmath @inbounds intensity_value += (t_delta_value * opacity_mat * c.sol * (c.sol * c.arad * temp_value[material_num]^4 - intensity_value))

                @inbounds @fastmath time_in[material_num] += t_delta_value

                last_material_num = material_num
            end
            local frac::Vector{Float64} = @. @fastmath time_in / c.delta_t
            # Update temperature values
            @. @fastmath temp_value += (intensity_deposition - c.delta_t * sigma_a(c.opacity, temp_value) * c.arad * c.sol * temp_value^4) / (rho(c.dens, temp_value) * c_v(c.spec_heat, temp_value))

            RunningStatistics.push(stat_1_temp[time_index], temp_value[1])
            RunningStatistics.push(stat_2_temp[time_index], temp_value[2])

            # Push values statistically
            if @inbounds @fastmath(frac[1] > 0.0)
                @inbounds RunningStatistics.push(stat_1_intensity[time_index], intensity_value * frac[1])
            end
            if @inbounds @fastmath(frac[2] > 0.0)
                @inbounds RunningStatistics.push(stat_2_intensity[time_index], intensity_value * frac[2])
            end
        end
    end

    # Mean values
    local material_1_intensity::Vector{Float64} = RunningStatistics.mean.(stat_1_intensity)
    local material_2_intensity::Vector{Float64} = RunningStatistics.mean.(stat_2_intensity)
    local material_1_temp::Vector{Float64} = RunningStatistics.mean.(stat_1_temp)
    local material_2_temp::Vector{Float64} = RunningStatistics.mean.(stat_2_temp)

    # Variance values
    local variance_1_intensity::Vector{Float64} = RunningStatistics.variance.(stat_1_intensity)
    local variance_2_intensity::Vector{Float64} = RunningStatistics.variance.(stat_2_intensity)
    local variance_1_temp::Vector{Float64} = RunningStatistics.variance.(stat_1_temp)
    local variance_2_temp::Vector{Float64} = RunningStatistics.variance.(stat_2_temp)

    # Maximum and minimum values
    local max_intensity_1::Vector{Float64} = RunningStatistics.greatest.(stat_1_intensity)
    local min_intensity_1::Vector{Float64} = RunningStatistics.least.(stat_1_intensity)
    local max_intensity_2::Vector{Float64} = RunningStatistics.greatest.(stat_2_intensity)
    local min_intensity_2::Vector{Float64} = RunningStatistics.least.(stat_2_intensity)
    local max_temp_1::Vector{Float64} = RunningStatistics.greatest.(stat_1_temp)
    local min_temp_1::Vector{Float64} = RunningStatistics.least.(stat_1_temp)
    local max_temp_2::Vector{Float64} = RunningStatistics.greatest.(stat_2_temp)
    local min_temp_2::Vector{Float64} = RunningStatistics.least.(stat_2_temp)

    local tabular::DataFrame = DataFrame(time=times, intensity1=material_1_intensity, varintensity1=variance_1_intensity, maxintensity1=max_intensity_1, minintensity1=min_intensity_1, temperature1=energy_from_temp.(material_1_temp), vartemperature1=energy_from_temp.(variance_1_temp), maxtemperature1=energy_from_temp.(max_temp_1), mintemperature1=energy_from_temp.(min_temp_1), intensity2=material_2_intensity, varintensity2=variance_2_intensity, maxintensity2=max_intensity_2, minintensity2=min_intensity_2, temperature2=energy_from_temp.(material_2_temp), vartemperature2=energy_from_temp.(variance_2_temp), maxtemperature2=energy_from_temp.(max_temp_2), mintemperature2=energy_from_temp.(min_temp_2))

    CSV.write("out/nonlinear.csv", tabular)

    return nothing
end

main()
