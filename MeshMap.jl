module MeshMap
    export struct_to_unstruct, unstruct_to_struct, material_calc

    set_zero_subnormals(true)

    function struct_to_unstruct(structured::Vector{Float64}, struct_delta::Float64, struct_size::Int64, unstruct_delta::Vector{Float64}, unstruct_size::Int64)::Vector{Float64}
        # Assign counter to zero due to pre-fetch increment
        local counter::Int64 = 0

        # Assign tallies
        local distance_tally::Float64 = 0.0  # s
        local unstruct_distance_tally::Float64 = 0.0  # s
        local struct_distance_tally::Float64 = 0.0  # s
        local leftover_distance = 0.0  # s
        local delta::Float64 = 0.0  # s

        # Create the new unstructured array map
        local unstructured::Vector{Float64} = Vector{Float64}(undef, unstruct_size)

        # Loop for mapping the structured to the unstructured mesh
        for i::Int64 = 1:unstruct_size
            # Start from the border with a zero weight for the current cell
            local distance_overlap::Bool = false
            local weight_tally::Float64 = 0.0  # unit*s

            # Increment unstructured distance tally
            @inbounds @fastmath unstruct_distance_tally += unstruct_delta[i]  # s

            # Carry over leftover distance
            if @fastmath(leftover_distance != 0.0)
                # If leftover_distance is still over-reaching the tally boundaries
                if @fastmath((distance_tally + leftover_distance) >= unstruct_distance_tally)
                    @inbounds delta = unstruct_delta[i]  # s
                    leftover_distance = @fastmath struct_distance_tally - unstruct_distance_tally  # s
                    distance_overlap = true
                else
                    delta = leftover_distance  # s
                    leftover_distance = 0.0  # s
                end
                @inbounds @fastmath weight_tally += delta * structured[counter]  # unit*s
                @fastmath distance_tally += delta  # s
            end

            while @fastmath((!distance_overlap) && (counter < struct_size))
                # Increment counter (structured index)
                @fastmath counter += 1

                # Increment structured distance tally
                @fastmath struct_distance_tally += struct_delta  # s

                # Check for boundary overlap
                if @fastmath((struct_distance_tally >= unstruct_distance_tally) || (counter == struct_size))
                    delta = @fastmath unstruct_distance_tally - distance_tally  # s
                    leftover_distance = @fastmath struct_distance_tally - unstruct_distance_tally  # s
                    distance_overlap = true
                else
                    delta = struct_delta  # s
                end

                # Increment the known distance tally
                distance_tally += delta  # s

                # Apply linear weighted tally
                @inbounds @fastmath weight_tally += delta * structured[counter]  # unit*s
            end  # Structured loop
            @inbounds @fastmath unstructured[i] = weight_tally / unstruct_delta[i]  # unit
        end  # Unstructured loop
        return unstructured
    end

    function unstruct_to_struct(unstructured::Vector{Float64}, unstruct_delta::Vector{Float64}, unstruct_size::Int64, struct_delta::Float64, struct_size::Int64)::Vector{Float64}
        # Assign counter to zero due to pre-fetch increment
        local counter::Int64 = 0

        # Assign tallies
        local distance_tally::Float64 = 0.0  # s
        local unstruct_distance_tally::Float64 = 0.0  # s
        local struct_distance_tally::Float64 = 0.0  # s
        local leftover_distance::Float64 = 0.0  # s
        local delta::Float64 = 0.0  # s

        # Create the new structured array map
        local structured::Vector{Float64} = Vector{Float64}(undef, struct_size)

        # Loop for mapping the unstructured to the structured mesh
        for i::Int64 = 1:struct_size
            # Start from the border with a zero weight for the current cell
            local distance_overlap::Bool = false
            local weight_tally::Float64 = 0.0  # unit*s

            # Increment structured distance tally
            @fastmath struct_distance_tally += struct_delta  # s

            # Carry over leftover distance
            if @fastmath(leftover_distance > 0.0)
                # If leftover_distance is still over-reaching the tally boundaries
                if @fastmath((distance_tally + leftover_distance) >= struct_distance_tally)
                    delta = struct_delta  # s
                    leftover_distance = @fastmath unstruct_distance_tally - struct_distance_tally  # s
                    distance_overlap = true
                else
                    delta = leftover_distance  # s
                    leftover_distance = 0.0  # s
                end
                @inbounds @fastmath weight_tally += delta * unstructured[counter]  # unit*s
                @fastmath distance_tally += delta  # s
            end

            while @fastmath((!distance_overlap) && (counter < unstruct_size))
                # Increment counter (unstructured index)
                @fastmath counter += 1

                # Increment unstructured distance tally
                @fastmath @inbounds unstruct_distance_tally += unstruct_delta[counter]  # s

                # Check for boundary overlap
                if @fastmath(unstruct_distance_tally >= struct_distance_tally || (counter == unstruct_size))
                    delta = @fastmath struct_distance_tally - distance_tally  # s
                    leftover_distance = @fastmath unstruct_distance_tally - struct_distance_tally  # s
                    distance_overlap = true
                else
                    @inbounds delta = unstruct_delta[counter]  # s
                end

                # Increment the known distance tally
                @fastmath distance_tally += delta  # s

                # Apply linear weighted tally
                @inbounds @fastmath weight_tally += delta * unstructured[counter]  # unit*s
            end  # Untructured loop
            @inbounds @fastmath structured[i] = weight_tally / struct_delta  # unit
        end  # Structured loop
        return structured
    end

    function material_calc(unstructured::Vector{Float64}, unstruct_delta::Vector{Float64}, unstruct_size::Int64, materials::Vector{Int32}, struct_delta::Float64, struct_size::Int64, num_materials::Int32)::Array{Float64, 2}
        # Create the new structured array map
        local material_struct::Array{Float64, 2} = zeros(Float64, struct_size, num_materials)
        for k::Int32 = 1:num_materials
            # Assign counter to zero due to pre-fetch increment
            local counter::Int64 = 0

            # Assign tallies
            local distance_tally::Float64 = 0.0  # s
            local unstruct_distance_tally::Float64 = 0.0  # s
            local struct_distance_tally::Float64 = 0.0  # s
            local leftover_distance::Float64 = 0.0  # s
            local switch::Float64 = 0.0
            local delta::Float64 = 0.0  # s

            # Loop for mapping the unstructured to the structured mesh
            for i::Int64 = 1:struct_size
                # Start from the border with a zero weight for the current cell
                local distance_overlap::Bool = false
                local weight_tally::Float64 = 0.0  # unit*s

                # Increment structured distance tally
                @fastmath struct_distance_tally += struct_delta  # s

                # Carry over leftover distance
                if @fastmath(leftover_distance > 0.0)
                    # If leftover distance is still over-reaching tally boundaries
                    if @fastmath((distance_tally + leftover_distance) >= struct_distance_tally)
                        delta = struct_delta  # s
                        leftover_distance = @fastmath unstruct_distance_tally - struct_distance_tally  # s
                        distance_overlap = true
                    else
                        delta = leftover_distance  # s
                        leftover_distance = 0.0  # s
                    end
                    @inbounds @fastmath weight_tally += delta * unstructured[counter] * switch  # unit*cm
                    @fastmath distance_tally += delta  # s
                end

                while @fastmath((!distance_overlap) && (counter < unstruct_size))
                    # Increment counter (unstructured index)
                    @fastmath counter += 1
                    # Increment unstructured distance tally
                    @fastmath @inbounds unstruct_distance_tally += unstruct_delta[counter]  # s

                    # Material number for calculations
                    # Tally switch for each material
                    if @inbounds @fastmath(materials[counter] == k)
                        switch = 1.0
                    else
                        switch = 0.0
                    end

                    # Check for boundary overlap
                    if @fastmath((unstruct_distance_tally >= struct_distance_tally) || (counter == unstruct_size))
                        delta = @fastmath struct_distance_tally - distance_tally  # s
                        leftover_distance = @fastmath unstruct_distance_tally - struct_distance_tally  # s
                        distance_overlap = true
                    else
                        @inbounds delta = unstruct_delta[counter]  # s
                    end

                    # Increment the known distance tally
                    @fastmath distance_tally += delta  # s

                    # Apply linear weighted tally
                    @inbounds @fastmath weight_tally += delta * unstructured[counter] * switch  # unit*s
                end  # Untructured loop

                # Average the results, or just append if no results previously
                @inbounds @fastmath material_struct[i, k] += (weight_tally / struct_delta)  # unit
            end  # Structured loop
        end  # Material loop
        return material_struct
    end
end
