include("config64.jl")

module random
    using config64
    using Distributions

    export sampleMultiNomial, sampleBeta

    function sampleMultiNomial(n::I, p::Array{R, 1})::Array{I, 1} where {I <: Integer, R <: Real}
        d = Multinomial(n, p)
        return vec(rand(d,1))
    end

    function sampleBeta(α::R, β::R, epsilon::R = 0.01)::R where {R <: Real}
        d = Beta(α, β)
        return rand(d,1)[1]
    end
    function sampleTruncatedBeta(α::R, β::R, epsilon::R = 0.01)::R where {R <: Real}
        d = Truncated(Distributions.Beta{R}(α, β), epsilon, (R)(1.0)-epsilon )
        # d = Beta(α, β)
        return rand(d,1)[1]
    end

end
