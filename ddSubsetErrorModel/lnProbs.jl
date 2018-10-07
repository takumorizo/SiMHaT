Include("config.jl")

module lnProbs
    using ..config
    using Distributions
    using SpecialFunctions
    # x ∈ {0,1}, x = 1 w.p. p
    function ln_P_ber(x::INT, p::REAL)
        return (REAL)(x * log(e, p) + (1-x) * log(e, (1-p)))
    end

    function ln_P_beta(p::REAL, α::REAL, β::REAL)
        return (REAL)((α-1)*log(e, p) + (β-1)*log(e,(1-p)) - SpecialFunctions.lbeta(α, β))
    end
end
