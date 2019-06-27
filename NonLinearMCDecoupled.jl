#!/usr/bin/env julia

include("Constants.jl")
const c = Constants

include("RunningStatistics.jl")
include("PhysicsFunctions.jl")
using .PhysicsFunctions
using Random
using DataFrames
using CSV
using ProgressMeter

function main()::Nothing
    set_zero_subnormals(true)

    # Iteration condition
    local generator::MersenneTwister = MersenneTwister(1234)

    # Computational values
    local times::Vector{Float64} = @fastmath [(x * c.delta_t + c.t_init) * c.sol for x in 1:c.num_t]

    # Probability for material sampling
    local prob_1::Float64 = c.volfrac[1]
    local change_prob_1::Float64 = @fastmath 1.0 / c.chord[1] * c.delta_t
    local change_prob_2::Float64 = @fastmath 1.0 / c.chord[2] * c.delta_t
    if @fastmath((change_prob_1 > 1.0) || (change_prob_2 > 1.0))
        println("The value for delta_t is too large for sampling")
        return nothing
    end

    # Parallel arrays
    local stat_1_intensity::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:c.num_t]
    local stat_2_intensity::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:c.num_t]
    local stat_1_temp::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:c.num_t]
    local stat_2_temp::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:c.num_t]
    local stat_1_opacity::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:c.num_t]
    local stat_2_opacity::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:c.num_t]

    # Outer loop
    @showprogress 1 for iteration_number=1:c.max_iterations
        local rand_num::Float64

        # First loop uses initial conditions
        local (intensity_value::Float64, temp_value_1::Float64, temp_value_2::Float64) = (c.init_intensity, c.init_temp, c.init_temp)

        # Sample initial starting material
        rand_num = @fastmath rand(generator, Float64)
        local material_num::Int32 = @fastmath (rand_num < prob_1) ? 1 : 2

        # Inner loop
        for (index, time) in enumerate(times)
            local (opacity_term::Float64, spec_heat_term::Float64, dens_term::Float64) = @inbounds (c.opacity[material_num], c.spec_heat[material_num], c.dens[material_num])
            local (change_prob::Float64, temp_value::Float64) = (material_num == 1) ? (change_prob_1, temp_value_1) : (change_prob_2, temp_value_2)

            # Sample whether material changes
            rand_num = @fastmath rand(generator, Float64)

            if @fastmath(rand_num > change_prob)
                local opacity::Float64 = PhysicsFunctions.sigma_a(opacity_term, temp_value)
                local spec_heat::Float64 = PhysicsFunctions.c_v(spec_heat_term, temp_value)
                local dens::Float64 = PhysicsFunctions.rho(dens_term, temp_value)

                local new_intensity_value::Float64 = PhysicsFunctions.balance_intensity(opacity, c.delta_t, intensity_value, temp_value)
                local new_temp_value::Float64 = PhysicsFunctions.balance_temp(opacity, spec_heat, dens, c.delta_t, intensity_value, temp_value)

                (intensity_value, temp_value) = (new_intensity_value, new_temp_value)

                if @fastmath(material_num == 1)
                    temp_value_1 = temp_value  # eV
                    @inbounds RunningStatistics.push(stat_1_intensity[index], intensity_value)  # erg/cm^2-s
                    @inbounds RunningStatistics.push(stat_1_temp[index], temp_value)  # eV
                    @inbounds RunningStatistics.push(stat_1_opacity[index], opacity)  # cm^-1
                else
                    temp_value_2 = temp_value  # eV
                    @inbounds RunningStatistics.push(stat_2_intensity[index], intensity_value)  # erg/cm^2-s
                    @inbounds RunningStatistics.push(stat_2_temp[index], temp_value)  # eV
                    @inbounds RunningStatistics.push(stat_2_opacity[index], opacity)  # cm^-1
                end
            else
                material_num = @fastmath (material_num == 1) ? 2 : 1
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
    local max_opacity_1::Vector{Float64} = RunningStatistics.greatest.(stat_1_opacity)
    local min_opacity_1::Vector{Float64} = RunningStatistics.least.(stat_1_opacity)
    local max_opacity_2::Vector{Float64} = RunningStatistics.greatest.(stat_2_opacity)
    local min_opacity_2::Vector{Float64} = RunningStatistics.least.(stat_2_opacity)

    local tabular::DataFrame = DataFrame(time=times, intensity1=material_1_intensity, varintensity1=variance_1_intensity, maxintensity1=max_intensity_1, minintensity1=min_intensity_1, temperature1=energy_from_temp.(material_1_temp), vartemperature1=energy_from_temp.(variance_1_temp), maxtemperature1=energy_from_temp.(max_temp_1), mintemperature1=energy_from_temp.(min_temp_1), intensity2=material_2_intensity, varintensity2=variance_2_intensity, maxintensity2=max_intensity_2, minintensity2=min_intensity_2, temperature2=energy_from_temp.(material_2_temp), vartemperature2=energy_from_temp.(variance_2_temp), maxtemperature2=energy_from_temp.(max_temp_2), mintemperature2=energy_from_temp.(min_temp_2), maxopacity1=max_opacity_1, minopacity1=min_opacity_1, maxopacity2=max_opacity_2, minopacity2=min_opacity_2)

    CSV.write("out/nonlinearmc.csv", tabular)

    return nothing
end

main()
