include("config64.jl")

module lnProbs
    using config64
    using Distributions
    # x ∈ {0,1}, x = 1 w.p. p
    function ln_P_ber(x::INT, p::REAL)
        return x * log(e, p) + (1-x) * log(e, (1-p))
    end

    function ln_P_beta(p::REAL, α::REAL, β::REAL)
        return (α-1)*log(e, p) + (β-1)*log(e,(1-p)) - lbeta(α, β)
    end
end
