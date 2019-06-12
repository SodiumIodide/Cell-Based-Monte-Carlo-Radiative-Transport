#!/usr/bin/env julia
# BuildMaterialPDF.jl

#=
Really only exists to ensure that the CDF for sampling is computed appropriately.

Relies on a comparison with known exponential distributions from the Distributions package, and the fitting methods therein.
=#

set_zero_subnormals(true)

using Random
using Distributions
const d = Distributions

include("Constants.jl")
using .Constants

include("Common.jl")
using .Common
const com = Common

function main()::Nothing
    generator = Random.MersenneTwister(1234)
    num_samples = convert(Int64, 1e6)
    particle_1 = com.build_particle(generator, 1, 0.0)
    particle_2 = com.build_particle(generator, 2, 0.0)

    samples_1 = [com.dist_to_transition(generator, particle_1) for p in 1:num_samples]
    samples_2 = [com.dist_to_transition(generator, particle_2) for p in 1:num_samples]

    lambda_1 = d.fit(d.Exponential, samples_1)
    lambda_2 = d.fit(d.Exponential, samples_2)

    println("Chord 1: ", d.params(lambda_1)[1])
    println("Chord 2: ", d.params(lambda_2)[1])

    return nothing
end

main()
