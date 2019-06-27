module PhysicsFunctions
    export sigma_a, c_v, rho, energy_from_temp, temp_from_energy

    using ..Constants
    const c = Constants

    set_zero_subnormals(true)

    # More generic functions to make use of complex step differentiation for the Jacobian term
    @inline function sigma_a(opacity_term::Float64, temp)
        return @fastmath opacity_term
    end

    @inline function c_v(spec_heat_term::Float64, temp)
        return @fastmath spec_heat_term * temp^3
    end

    @inline function rho(dens_term::Float64, temp)
        return @fastmath dens_term
    end

    # Conversion functions
    @inline function energy_from_temp(temp::Float64)::Float64
        return @fastmath c.arad * temp^4
    end

    @inline function temp_from_energy(eng::Float64)::Float64
        return @fastmath (eng / c.arad)^(0.25)
    end

    # Monte Carlo state functions
    function balance_intensity(opacity::Float64, dt::Float64, past_intensity::Float64, past_temp::Float64)::Float64
        set_zero_subnormals(true)
        local term_1::Float64 = @fastmath dt * c.sol^2 * opacity * c.arad * past_temp^4
        local term_2::Float64 = @fastmath (1.0 - dt * c.sol * opacity) * past_intensity

        return @fastmath term_1 + term_2
    end

    function balance_temp(opacity::Float64, spec_heat::Float64, density::Float64, dt::Float64, past_intensity::Float64, past_temp::Float64)::Float64
        set_zero_subnormals(true)
        local term_1::Float64 = @fastmath dt / (density * spec_heat) * opacity * past_intensity
        local term_2::Float64 = @fastmath - dt / (density * spec_heat) * c.sol * opacity * c.arad * past_temp^4
        local term_3::Float64 = past_temp

        return @fastmath term_1 + term_2 + term_3
    end

    # Realization functions
    # Realization numerical functions
    function balance_a(intensity::Float64, temp::Float64, delta_t::Float64, opacity::Float64, past_intensity::Float64)::Float64
        local term_1::Float64 = intensity
        local term_2::Float64 = @fastmath - delta_t * c.sol^2 * opacity * c.arad * temp^4
        local term_3::Float64 = @fastmath delta_t * c.sol * opacity * intensity
        local term_4::Float64 = @fastmath - past_intensity

        return @fastmath term_1 + term_2 + term_3 + term_4
    end

    function balance_b(intensity::Float64, temp::Float64, delta_t::Float64, opacity::Float64, spec_heat::Float64, density::Float64, past_temp::Float64)::Float64
        local term_1::Float64 = temp
        local term_2::Float64 = @fastmath - delta_t / (density * spec_heat) * opacity * intensity
        local term_3::Float64 = @fastmath delta_t / (density * spec_heat) * c.sol * opacity * c.arad * temp^4
        local term_4::Float64 = @fastmath - past_temp

        return @fastmath term_1 + term_2 + term_3 + term_4
    end

    # Jacobian computation for Newton iteration using complex step differentiation
    function complex_step_jacobian(intensity::Float64, temp::Float64, delta_t::Float64, opacity_term::Float64, density_term::Float64, spec_heat_term::Float64, past_intensity::Float64, past_temp::Float64)::Array{Float64, 2}
        local step::Float64 = 1e-8
        local complex_temp::Complex{Float64} = complex(temp, step)
        local complex_intensity::Complex{Float64} = complex(intensity, step)
        local r_opacity::Float64 = sigma_a(opacity_term, temp)
        local c_opacity::Complex{Float64} = sigma_a(opacity_term, complex_temp)
        local r_spec_heat::Float64 = c_v(spec_heat_term, temp)
        local c_spec_heat::Complex{Float64} = c_v(spec_heat_term, complex_temp)
        local r_density::Float64 = rho(density_term, temp)
        local c_density::Complex{Float64} = rho(density_term, complex_temp)

        local term_1::Float64 = @fastmath imag(complex_intensity - delta_t * c.sol^2 * r_opacity * c.arad * temp^4 + delta_t * c.sol * r_opacity * complex_intensity - past_intensity) / step
        local term_2::Float64 = @fastmath imag(intensity - delta_t * c.sol^2 * c_opacity * c.arad * complex_temp^4 + delta_t * c.sol * c_opacity * intensity - past_intensity) / step
        local term_3::Float64 = @fastmath imag(temp - delta_t / (r_density * r_spec_heat) * r_opacity * complex_intensity + delta_t / (r_density * r_spec_heat) * c.sol * r_opacity * c.arad * temp^4 - past_temp) / step
        local term_4::Float64 = @fastmath imag(complex_temp - delta_t / (c_density * c_spec_heat) * c_opacity * intensity + delta_t / (c_density * c_spec_heat) * c.sol * c_opacity * c.arad * complex_temp^4 - past_temp) / step

        return [
            term_1 term_2;
            term_3 term_4
        ]
    end

    function relative_change(delta::Vector{Float64}, original_terms::Vector{Float64})::Float64
        local current_norm::Float64 = @fastmath sqrt(sum([x^2 for x in delta]))
        local original_norm::Float64 = @fastmath sqrt(sum([x^2 for x in original_terms]))

        return @fastmath current_norm / original_norm
    end
end
