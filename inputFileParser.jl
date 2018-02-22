include("config64.jl")

module inputParser
    using ConfParser
    using config64

    export Parameters
    type Parameters{I <: Integer, R <: Real}
        α_s::R
        α_v::R # dirichlet process hyper parameter for sample/variant

        minXIn::I
        minYIn::I
        minXOut::I
        minYOut::I

        λs::Array{R,1} # []λ_e::R, λ_m::R, λ_t::R] frequency of error/merged/tree for each block
        βs::Array{R,2}
        # β_e_0::R, β_e_1::R # beta hyper parameter for sample-with-variant frequency in error block
        # β_m_0::R, β_m_1::R # beta hyper parameter for sample-with-variant frequency in merged block
        # β_t_0::R, β_t_1::R # beta hyper parameter for sample-with-variant frequency in tree block
        ln_G_t::Array{R,1}
        ln_G_B::Array{R,1}
    end

    function parseConfigFile(configFile::String)::Parameters{INT, REAL}
        conf = ConfParse(configFile)
        parse_conf!(conf)
        println(conf)
        α_s   = parse(REAL, String(retrieve(conf, "model", "alpha_s")) )
        α_v   = parse(REAL, String(retrieve(conf, "model", "alpha_v")) )

        minXIn   = parse(INT, String(retrieve(conf, "model", "minXIn")) )
        minYIn   = parse(INT, String(retrieve(conf, "model", "minYIn")) )
        minXOut  = parse(INT, String(retrieve(conf, "model", "minXOut")) )
        minYOut  = parse(INT, String(retrieve(conf, "model", "minYOut")) )

        λ_e   = parse(REAL, String(retrieve(conf, "model", "lambda_e")) )
        λ_m   = parse(REAL, String(retrieve(conf, "model", "lambda_m")) )
        λ_t   = parse(REAL, String(retrieve(conf, "model", "lambda_t")) )
        β_e_0 = parse(REAL, String(retrieve(conf, "model", "beta_e_0")) )
        β_e_1 = parse(REAL, String(retrieve(conf, "model", "beta_e_1")) )
        β_m_0 = parse(REAL, String(retrieve(conf, "model", "beta_m_0")) )
        β_m_1 = parse(REAL, String(retrieve(conf, "model", "beta_m_1")) )
        β_t_0 = parse(REAL, String(retrieve(conf, "model", "beta_t_0")) )
        β_t_1 = parse(REAL, String(retrieve(conf, "model", "beta_t_1")) )

        ln_1m_g1 = parse(REAL, String(retrieve(conf, "model", "ln_1m_g_t")) )
        ln_g1 = log(e, (REAL)(1.0) - exp(ln_1m_g1))

        ln_1m_g2 = parse(REAL, String(retrieve(conf, "model", "ln_1m_g_B")) )
        ln_g2 = log(e, (REAL)(1.0) - exp(ln_1m_g2))

        return Parameters{INT,REAL}(α_s, α_v,
                                minXIn, minYIn, minXOut, minYOut,
                                [λ_e, λ_m, λ_t],
                                [β_e_0 β_e_1; β_m_0 β_m_1; β_t_0 β_t_1;],
                                [ln_g1, ln_1m_g1], [ln_g2, ln_1m_g2])
    end

    function parseInputSummary(summaryPath::String)::Array{REAL, 2}
        ans = []
        isFirst = true
        open(summaryPath, "r") do f
            for line::String in eachline(f)
                if isFirst
                    line = strip(line)
                    isFirst = false
                    continue
                end
                line = strip(line)
                lineCol::Array{String,1} = split(line, '\t')
                nodeID   = parse(INT, lineCol[1])
                sampleID = parse(INT, lineCol[2])
                push!(ans, map(x -> parse(REAL, x), lineCol[3:length(lineCol)]) )
            end
        end
        S = length(ans)
        M = 0
        if S > 0
            M = length(ans[1])
        end
        mat = zeros(REAL, S, M)
        for (s,m) in Iterators.product(1:S,1:M)
            mat[s,m] = ans[s][m]
        end
        return mat
    end
end
