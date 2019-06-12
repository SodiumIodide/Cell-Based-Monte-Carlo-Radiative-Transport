module Common
    export Particle, CollisionType, isotropic_uvw, produce_time, map_material, get_material_properties, switch_material!, next_time!, dist_to_boundary, dist_to_census, dist_to_transition, move_particle!, boundary_hit!, build_particle

    using ..Random

    using ..Constants
    c = Constants

    set_zero_subnormals(true)

    mutable struct Particle
        # Weight value
        weight::Float64

        # Angle data
        u_val::Float64
        v_val::Float64
        w_val::Float64

        # Location data
        x_loc::Float64
        y_loc::Float64
        z_loc::Float64

        # Material properties
        mat_num::Int64
        opacity::Float64
        spec_heat::Float64
        dens::Float64
        chord::Float64

        # Time data
        t_remaining::Float64
    end

    @enum CollisionType begin
        x_bound = 1
        y_bound = 2
        z_bound = 3
    end

    @inline function isotropic_u(rand_num::Float64)::Float64
        # Random number between 0 and 1
        return @fastmath 2.0 * rand_num - 1.0
    end

    @inline function isotropic_v(rand_num::Float64, iso_u::Float64)::Float64
        # Random number between 0 and 1
        return @fastmath sqrt(1.0 - iso_u^2) * cos(2.0 * pi * rand_num)
    end

    @inline function isotropic_w(rand_num::Float64, iso_u::Float64)::Float64
        # Random number between 0 and 1
        return @fastmath sqrt(1.0 - iso_u^2) * sin(2.0 * pi * rand_num)
    end

    @inline function isotropic_uvw(gen::MersenneTwister)::NTuple{3, Float64}
        local u_value::Float64 = @fastmath isotropic_u(rand(gen))
        local rand_val::Float64 = @fastmath rand(gen)
        return @fastmath (u_value, isotropic_v(rand_val, u_value), isotropic_w(rand_val, u_value))
    end

    @inline function produce_time(rand_num::Float64)::Float64
        return @fastmath rand_num * c.delta_t
    end

    @inline function get_material_properties(mat_num::Int64)::NTuple{4, Float64}
        return @inbounds (c.opacity[mat_num], c.spec_heat[mat_num], c.dens[mat_num], c.chord[mat_num])
    end

    @inline function map_material(rand_num::Float64)::Int64
        @fastmath return (rand_num < c.volfrac[1]) ? 1 : 2
    end

    @inline function generate_position(gen::MersenneTwister)::NTuple{3, Float64}
        return @fastmath (c.x_len * rand(gen), c.y_len * rand(gen), c.z_len * rand(gen))
    end

    function update_material!(particle::Particle, material::Int64)::Nothing
        if (material == 1)
            particle.mat_num = 1
            particle.opacity = @inbounds c.opacity[1]
            particle.spec_heat = @inbounds c.spec_heat[1]
            particle.dens = @inbounds c.dens[1]
            particle.chord = @inbounds c.chord[1]
        else
            particle.mat_num = 2
            particle.opacity = @inbounds c.opacity[2]
            particle.spec_heat = @inbounds c.spec_heat[2]
            particle.dens = @inbounds c.dens[2]
            particle.chord = @inbounds c.chord[2]
        end

        return nothing
    end

    @inline function switch_material!(particle::Particle)::Nothing
        @fastmath (particle.mat_num == 1) ? update_material!(particle, 2) : update_material!(particle, 1)

        return nothing
    end

    @inline function next_time!(particle::Particle)::Nothing
        particle.t_remaining = c.delta_t

        return nothing
    end

    function dist_to_boundary(particle::Particle)::Tuple{Float64, CollisionType}
        local (x_dist::Float64, y_dist::Float64, z_dist::Float64)
        if @fastmath (particle.u_val > 0.0)
            x_dist = @fastmath (c.x_len - particle.x_loc) / particle.u_val
        else
            x_dist = @fastmath abs(particle.x_loc / particle.u_val)
        end
        if @fastmath (particle.v_val > 0.0)
            y_dist = @fastmath (c.y_len - particle.y_loc) / particle.v_val
        else
            y_dist = @fastmath abs(particle.y_loc / particle.v_val)
        end
        if @fastmath (particle.w_val > 0.0)
            z_dist = @fastmath (c.z_len - particle.z_loc) / particle.w_val
        else
            z_dist = @fastmath abs(particle.z_loc / particle.w_val)
        end

        local min_col_dist::Float64 = @fastmath minimum([x_dist, y_dist, z_dist])
        if @fastmath (min_col_dist == x_dist)
            return (x_dist, x_bound)
        elseif @fastmath (min_col_dist == y_dist)
            return (y_dist, y_bound)
        else
            return (z_dist, z_bound)
        end
    end

    @inline function dist_to_transition(gen::MersenneTwister, particle::Particle)::Float64
        return @fastmath -(particle.chord * c.sol) * log(rand(gen))
    end

    @inline function dist_to_census(particle::Particle)::Float64
        return @fastmath particle.t_remaining * c.sol
    end

    function move_particle!(particle::Particle, distance::Float64)::Nothing
        @fastmath particle.x_loc += distance * particle.u_val
        @fastmath particle.y_loc += distance * particle.v_val
        @fastmath particle.z_loc += distance * particle.w_val
        @fastmath particle.t_remaining -= distance / c.sol

        return nothing
    end

    function boundary_hit!(particle::Particle, bound::CollisionType)::Nothing
        # Change direction of travel (all surfaces reflective)
        if @fastmath (bound == x_bound)
            particle.u_val = -particle.u_val
        elseif @fastmath (bound == y_bound)
            particle.v_val = -particle.v_val
        else
            particle.w_val = -particle.w_val
        end
    
        return nothing
    end

    @inline function build_particle(gen::MersenneTwister, mat_num::Int64, weight::Float64, time_left::Float64=c.delta_t)::Particle
        return Particle(weight, isotropic_uvw(gen)..., generate_position(gen)..., mat_num, get_material_properties(mat_num)..., time_left)
    end
end
