module PhysicsFunctions
    export sigma_a, c_v, energy_from_temp, temp_from_energy

    using ..Constants
    const c = Constants

    set_zero_subnormals(true)

    @inline function sigma_a(opacity_term::Float64, temp::Float64)::Float64
        return @fastmath opacity_term
    end

    @inline function c_v(spec_heat_term::Float64, temp::Float64)::Float64
        return @fastmath spec_heat_term * temp^3
    end

    @inline function energy_from_temp(temp::Float64)::Float64
        return @fastmath c.arad * temp^4
    end

    @inline function temp_from_energy(eng::Float64)::Float64
        return @fastmath (eng / c.arad)^(0.25)
    end
end
