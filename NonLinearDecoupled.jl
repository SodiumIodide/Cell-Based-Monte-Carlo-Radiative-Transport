#!/usr/bin/env julia

include("Constants.jl")
const c = Constants

include("GeometryGen.jl")
include("MeshMap.jl")
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
    @showprogress 1 for iteration_number=1:c.max_iterations
        local (t_delta::Vector{Float64}, t_arr::Vector{Float64}, materials::Vector{Int32}, num_cells::Int64) = GeometryGen.get_geometry(c.chord[1], c.chord[2], c.t_max, c.num_divs, rng=generator)

        local intensity::Vector{Float64} = zeros(num_cells)
        local temp::Vector{Float64} = zeros(num_cells)

        # First loop uses initial conditions
        local (intensity_value::Float64, temp_value_1::Float64, temp_value_2::Float64) = (c.init_intensity, c.init_temp, c.init_temp)

        # Inner loop
        for (index, material) in enumerate(materials)
            local delta_t_unstruct::Float64 = @inbounds t_delta[index]
            local (opacity_term::Float64, spec_heat_term::Float64, dens_term::Float64) = @inbounds (c.opacity[material], c.spec_heat[material], c.dens[material])
            local temp_value::Float64 = (material == 1) ? temp_value_1 : temp_value_2

            local original_terms::Vector{Float64} = vec([
                intensity_value
                temp_value
            ])
            local old_terms::Vector{Float64} = deepcopy(original_terms)
            local new_terms::Vector{Float64} = deepcopy(original_terms)

            local error::Float64 = 1.0

            # Newtonian loop
            while @fastmath(error >= c.tolerance)
                local opacity::Float64 = @inbounds PhysicsFunctions.sigma_a(opacity_term, old_terms[2])
                local spec_heat::Float64 = @inbounds PhysicsFunctions.c_v(spec_heat_term, old_terms[2])
                local dens::Float64 = @inbounds PhysicsFunctions.rho(dens_term, old_terms[2])

                local jacobian::Array{Float64, 2} = @inbounds PhysicsFunctions.complex_step_jacobian(old_terms[1], old_terms[2], delta_t_unstruct, opacity_term, dens_term, spec_heat_term, intensity_value, temp_value)
                local func_vector::Vector{Float64} = [
                    @inbounds PhysicsFunctions.balance_a(old_terms[1], old_terms[2], delta_t_unstruct, opacity, intensity_value)
                    @inbounds PhysicsFunctions.balance_b(old_terms[1], old_terms[2], delta_t_unstruct, opacity, spec_heat, dens, temp_value)
                ]

                local delta::Vector{Float64} = @inbounds @fastmath jacobian \ - func_vector

                @inbounds new_terms = @fastmath delta + old_terms
                old_terms = deepcopy(new_terms)

                error = PhysicsFunctions.relative_change(delta, original_terms)
            end
            intensity_value = @inbounds new_terms[1]
            @inbounds intensity[index] = intensity_value  # erg/cm^2-s
            if (material == 1)
                temp_value_1 = @inbounds new_terms[2]
                temp[index] = temp_value_1  # eV
            else
                temp_value_2 = @inbounds new_terms[2]
                temp[index] = temp_value_2  # eV
            end
        end

        local material_intensity_array::Array{Float64, 2} = MeshMap.material_calc(intensity, t_delta, num_cells, materials, c.delta_t, c.num_t, convert(Int32, 2))
        local material_temp_array::Array{Float64, 2} = MeshMap.material_calc(temp, t_delta, num_cells, materials, c.delta_t, c.num_t, convert(Int32, 2))

        @simd for k in 1:c.num_t
            if @inbounds @fastmath(material_intensity_array[k, 1] != 0.0)
                @inbounds RunningStatistics.push(stat_1_intensity[k], material_intensity_array[k, 1])  # erg/cm^2-s
                @inbounds RunningStatistics.push(stat_1_temp[k], material_temp_array[k, 1])  # eV
            end
            if @inbounds @fastmath(material_intensity_array[k, 2] != 0.0)
                @inbounds RunningStatistics.push(stat_2_intensity[k], material_intensity_array[k, 2])  # erg/cm^2-s
                @inbounds RunningStatistics.push(stat_2_temp[k], material_temp_array[k, 2])  # eV
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
