module PhysicsFunctions
    export sigma_a, c_v

    set_zero_subnormals(true)

    @inline function sigma_a(opacity_term::Float64, temp::Float64)::Float64
        return @fastmath opacity_term / temp^3
    end

    @inline function c_v(spec_heat_term::Float64, temp::Float64)::Float64
        return spec_heat_term
    end
end
