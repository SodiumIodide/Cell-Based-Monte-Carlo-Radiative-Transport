#!/usr/bin/env julia

# diagonalization.jl

using LinearAlgebra
using ProgressMeter
using DataFrames
using CSV

include("Constants.jl")
const c = Constants

function main()::Nothing
    # Problem constants
    local factor::Vector{Float64} = @fastmath @. 4.0 * c.arad / (c.dens * c.spec_heat)

    # Initial conditions
    local init_intensity::Vector{Float64} = @. c.init_intensity * c.volfrac
    local init_energy::Vector{Float64} = @. c.arad * c.init_temp^4 * c.volfrac

    # Compute constant terms for matrix
    local term_1::Float64 = @fastmath @inbounds c.sol * c.opacity[1] + 1.0 / c.chord[1]
    local term_2::Float64 = @fastmath @inbounds - 1.0 / c.chord[2]
    local term_3::Float64 = @fastmath @inbounds - c.sol^2 * c.opacity[1]
    local term_4::Float64 = @fastmath @inbounds - 1.0 / c.chord[1]
    local term_5::Float64 = @fastmath @inbounds c.sol * c.opacity[2] + 1.0 / c.chord[2]
    local term_6::Float64 = @fastmath @inbounds - c.sol^2 * c.opacity[2]
    local term_7::Float64 = @fastmath @inbounds - factor[1] * c.opacity[1]
    local term_8::Float64 = @fastmath @inbounds factor[1] * c.opacity[1] * c.sol + 1.0 / c.chord[1]
    local term_9::Float64 = @fastmath @inbounds - 1.0 / c.chord[2]
    local term_10::Float64 = @fastmath @inbounds - factor[2] * c.opacity[2]
    local term_11::Float64 = @fastmath @inbounds - 1.0 / c.chord[1]
    local term_12::Float64 = @fastmath @inbounds factor[2] * c.opacity[2] * c.sol + 1.0 / c.chord[2]

    local coupled_matrix::Array{Float64,2} = [
        term_1 term_2 term_3 0.0;
        term_4 term_5 0.0 term_6;
        term_7 0.0 term_8 term_9;
        0.0 term_10 term_11 term_12
    ]

    # Obtain initial matrices
    local initial_kappa::Vector{Float64} = @inbounds vec([
        init_intensity[1],
        init_intensity[2],
        init_energy[1],
        init_energy[2]
    ])

    function compute(t::Float64)::Vector{Float64}
        return @fastmath exp(- coupled_matrix * t) * initial_kappa
    end

    # Obtain eigenvectors and eigenvalues
    local eigenvectors::Array{Float64,2} = @fastmath eigvecs(coupled_matrix)
    local eigenvalues::Vector{Float64} = @fastmath eigvals(coupled_matrix)

    local initial_chi::Vector{Float64} = @inbounds @fastmath eigenvectors \ initial_kappa

    function compute2(t::Float64)::Vector{Float64}
        local chi::Vector{Float64} = Vector{Float64}(undef, 4)
        for index in eachindex(eigenvalues)
            @inbounds chi[index] = @fastmath exp(-eigenvalues[index] * t) * initial_chi[index]
        end
        return @inbounds @fastmath eigenvectors * chi
    end

    local kappa::Vector{Float64} = Vector{Float64}(undef, 4)
    local intensity_1::Vector{Float64} = Vector{Float64}(undef, c.num_t)
    local intensity_2::Vector{Float64} = Vector{Float64}(undef, c.num_t)
    local energy_1::Vector{Float64} = Vector{Float64}(undef, c.num_t)
    local energy_2::Vector{Float64} = Vector{Float64}(undef, c.num_t)

    @showprogress 1 for t = 1:c.num_t
        kappa = compute((c.t_init + c.delta_t * t))
        #println(kappa)
        intensity_1[t] = kappa[1]
        intensity_2[t] = kappa[2]
        energy_1[t] = kappa[3]
        energy_2[t] = kappa[4]
    end

    # Take the implicit volume fraction into account
    @inbounds @fastmath intensity_1 ./= c.volfrac[1]  # erg/cm^2-s
    @inbounds @fastmath intensity_2 ./= c.volfrac[2]  # erg/cm^2-s
    @inbounds @fastmath energy_1 ./= c.volfrac[1]  # eV
    @inbounds @fastmath energy_2 ./= c.volfrac[2]

    local times::Vector{Float64} = @fastmath [(x * c.delta_t + c.t_init) * c.sol for x in 1:c.num_t]

    local tabular::DataFrame = DataFrame(time=times, intensity1=intensity_1, energy1=energy_1, intensity2=intensity_2, energy2=energy_2)

    CSV.write("out/analytic_linear.csv", tabular)

    return nothing
end

main()
