module RunningStatistics
    export RunningStat, clear, push, num, mean, variance, standard_deviation, least, greatest

    using Base.Threads

    set_zero_subnormals(true)

    mutable struct RunningStat
        m_n::Int64
        m_oldM::Float64
        m_newM::Float64
        m_oldS::Float64
        m_newS::Float64
        m_min::Float64
        m_max::Float64

        # Internal constructor
        RunningStat() = new(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end

    function clear(r::RunningStat)::Nothing
        r.m_n = 0
        r.m_oldM = r.m_newM = r.m_oldS = r.m_newS = r.m_min = r.m_max = 0.0

        return nothing
    end

    function push(r::RunningStat, x::Float64)::Nothing
        @fastmath r.m_n += 1

        # See Knuth TAOCP vol 2, 3rd edition, page 232
        if (r.m_n == 1)
            r.m_oldM = r.m_newM = x
            r.m_oldS = 0.0
            r.m_min = r.m_max = x
        else
            r.m_newM = @fastmath r.m_oldM + (x - r.m_oldM) / r.m_n
            r.m_newS = @fastmath r.m_oldS + (x - r.m_oldM) * (x - r.m_newM)

            # Set up for next iteration
            r.m_oldM = r.m_newM
            r.m_oldS = r.m_newS

            # Min and max
            r.m_min = @fastmath min(x, r.m_min)
            r.m_max = @fastmath max(x, r.m_max)
        end

        return nothing
    end

    function num(r::RunningStat)::Int64
        return r.m_n
    end

    function mean(r::RunningStat)::Float64
        return @fastmath (r.m_n > 0) ? r.m_newM : 0.0
    end

    function variance(r::RunningStat)::Float64
        return @fastmath (r.m_n > 1) ? (r.m_newS / (r.m_n - 1)) : 0.0
    end

    function standard_deviation(r::RunningStat)::Float64
        return @fastmath sqrt(variance(r))
    end

    function least(r::RunningStat)::Float64
        return r.m_min
    end

    function greatest(r::RunningStat)::Float64
        return r.m_max
    end
end
