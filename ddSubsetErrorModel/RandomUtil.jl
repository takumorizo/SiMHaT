
module RandomUtil
    using Distributions

    export sample_multinomial, sample_beta

    function sample_multinomial(n::I, p::Array{R, 1})::Array{I, 1} where {I <: Integer, R <: Real}
        d = Multinomial((I)(n), convert.(Float64, p) ./ sum(convert.(Float64, p)) )
        return vec(rand(d, 1))
    end

    function sample_beta(α::R, β::R, epsilon::R = (R)(0.01))::R where {R <: Real}
        d = Beta((Float64)(α), (Float64)(β))
        return rand(d, 1)[1]
    end
    function sample_truncatedbeta(α::R, β::R, epsilon::R = (R)(0.01))::R where {R <: Real}
        d = Truncated(Distributions.Beta{R}((Float64)(α), (Float64)(β)), (Float64)(epsilon), (Float64)((1.0)-epsilon) )
        # d = Beta(α, β)
        return rand(d, 1)[1]
    end

end
