#!/usr/bin/env julia

# backward_euler.jl

using ProgressMeter
using DataFrames
using CSV

include("Constants.jl")
const c = Constants

function main()::Nothing
    # Problem constants
    local factor::Vector{Float64} = @fastmath @. 4.0 * c.arad / (c.dens * c.spec_heat)  # eV^-3

    # Initial conditions
    local init_energy::Float64 = @fastmath c.arad * c.init_temp^4  # erg/cm^3

    # Compute constant terms for matrix
    local term_1::Float64 = @fastmath @inbounds 1.0 + c.delta_t * c.sol * c.opacity[1] + c.delta_t / c.chord[1]
    local term_2::Float64 = @fastmath @inbounds - c.delta_t / c.chord[2]
    local term_3::Float64 = @fastmath @inbounds - c.delta_t * c.sol^2 * c.opacity[1]
    local term_4::Float64 = @fastmath @inbounds - c.delta_t / c.chord[1]
    local term_5::Float64 = @fastmath @inbounds 1.0 + c.delta_t * c.sol * c.opacity[2] + c.delta_t / c.chord[2]
    local term_6::Float64 = @fastmath @inbounds - c.delta_t * c.sol^2 * c.opacity[2]
    local term_7::Float64 = @fastmath @inbounds - c.delta_t * factor[1] * c.opacity[1]
    local term_8::Float64 = @fastmath @inbounds 1.0 + c.delta_t * factor[1] * c.opacity[1] * c.sol + c.delta_t / c.chord[1]
    local term_9::Float64 = @fastmath @inbounds - c.delta_t / c.chord[2]
    local term_10::Float64 = @fastmath @inbounds - c.delta_t * factor[2] * c.opacity[2]
    local term_11::Float64 = @fastmath @inbounds - c.delta_t / c.chord[1]
    local term_12::Float64 = @fastmath @inbounds 1.0 + c.delta_t * factor[2] * c.opacity[2] * c.sol + c.delta_t / c.chord[2]

    # Set up problem definitions
    # This matrix is a constant parameter
    local coupled_matrix::Array{Float64,2} = [
        term_1 term_2 term_3 0.0;
        term_4 term_5 0.0 term_6;
        term_7 0.0 term_8 term_9;
        0.0 term_10 term_11 term_12
    ]

    # Vectors in problem
    local vector::Vector{Float64} = Vector{Float64}(undef, 4)
    local previous_vector::Array{Float64,1} = @fastmath @inbounds [
        *(c.init_intensity, c.volfrac[1])
        *(c.init_intensity, c.volfrac[2])
        *(init_energy, c.volfrac[1])
        *(init_energy, c.volfrac[2])
    ]

    # Preallocate vectors
    local intensity_1::Vector{Float64} = Vector{Float64}(undef, c.num_t)
    local intensity_2::Vector{Float64} = Vector{Float64}(undef, c.num_t)
    local energy_1::Vector{Float64} = Vector{Float64}(undef, c.num_t)
    local energy_2::Vector{Float64} = Vector{Float64}(undef, c.num_t)

    @showprogress 1 for t = 1:c.num_t
        vector = @inbounds @fastmath coupled_matrix \ previous_vector
        previous_vector = deepcopy(vector)
        @inbounds intensity_1[t] = previous_vector[1]
        @inbounds intensity_2[t] = previous_vector[2]
        @inbounds energy_1[t] = previous_vector[3]
        @inbounds energy_2[t] = previous_vector[4]
    end

    local times::Vector{Float64} = @fastmath [(x * c.delta_t + c.t_init) * c.sol for x in 1:c.num_t]

    # Take the implicit volume fraction into account
    @fastmath @inbounds intensity_1 ./= c.volfrac[1]  # erg/cm^2-s
    @fastmath @inbounds intensity_2 ./= c.volfrac[2]  # erg/cm^2-s
    @fastmath @inbounds energy_1 ./= c.volfrac[1]  # erg/cm^3
    @fastmath @inbounds energy_2 ./= c.volfrac[2]  # erg/cm^3

    local tabular::DataFrame = DataFrame(time=times, intensity1=intensity_1, energy1=energy_1, intensity2=intensity_2, energy2=energy_2)

    CSV.write("out/analytic_linear.csv", tabular)

    return nothing
end

main()
