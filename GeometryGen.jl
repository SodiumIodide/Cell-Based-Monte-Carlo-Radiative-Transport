module GeometryGen
    export get_geometry

    using Random

    set_zero_subnormals(true)

    function get_geometry(chord_a::Float64, chord_b::Float64, end_time::Float64, num_divs::Int64; rng::MersenneTwister=MersenneTwister(1234))::Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Float64}
        # Computational utilities
        local material_num::Int32
        local rand_num::Float64
        local chord::Float64

        # Updateable values
        local cons_time::Float64 = 0.0  # s
        local time::Float64 = 0.0  # s

        # Return values
        local x_delta::Vector{Float64} = Float64[]
        local materials::Vector{Int32} = Int32[]
        local x_arr::Vector{Float64} = Float64[]
        local num_cells::Int64 = 0

        # Sample size
        chord = @fastmath min(chord_a, chord_b)
        local sample_time::Float64 = @fastmath chord * (-log(0.5))
        local est_size::Int64 = convert(Int64, ceil(@fastmath(end_time / sample_time * num_divs)))

        sizehint!(x_delta, est_size)
        sizehint!(materials, est_size)
        sizehint!(x_arr, est_size)

        # Determine first material to use
        local prob_a::Float64 = @fastmath chord_a / (chord_a + chord_b)
        rand_num = @fastmath rand(rng, Float64)
        material_num = @fastmath (rand_num < prob_a) ? 1 : 2

        # Loop to build geometry
        while (cons_time < end_time)
            # Generate a random number
            rand_num = @fastmath rand(rng, Float64)

            # Assign a chord length based on material number
            chord = @fastmath (material_num == 1) ? chord_a : chord_b  # s

            # Calculate and append the material length
            time = @fastmath chord * (-log(rand_num))  # s
            @fastmath cons_time += time  # s

            # Check on thickness to not overshoot the boundary
            if @fastmath(cons_time > end_time)
                @fastmath time += end_time - cons_time  # s
                cons_time = end_time  # s
            end

            # Further discretize geometry
            for i::Int64 = 1:num_divs
                push!(x_delta, @fastmath(time / convert(Float64, num_divs)))
                push!(x_arr, @fastmath(cons_time - time + (time / convert(Float64, num_divs) * convert(Float64, i))))
                push!(materials, material_num)
                num_cells += @fastmath 1
            end

            # Update material number
            material_num = @fastmath (material_num == 1) ? 2 : 1
        end

        return (x_delta, x_arr, materials, num_cells)
    end
end
