Include("config.jl")

module InputParser
    using ConfParser
    using ..config

    export Parameters
    export Annealer

    mutable struct Annealer{I <: Integer, R <: Real}
        period::I
        ln_p_ladders::Array{R, 2}
    end

    mutable struct Parameters{I <: Integer, R <: Real}
        α_s::R
        α_v::R # dirichlet process hyper parameter for sample/variant

        δ_s::Array{R, 1} # freq vector for founder existence for B[c,m] = 1, 2, 3. (shared/merge/unique)
        λ_s::Array{R, 2}
        # λ_smu_0::R, λ_s_1::R # beta hyper parameter for variant frequency in shared/merge/unique block
        # λ_e_0::R,   λ_e_1::R # beta hyper parameter for variant frequency in error block
        β_s::Array{R, 1}
        # β_smu::R,  # ber hyper parameter for haplotype in shared/merge/unique block
        # β_e::R     # ber hyper parameter for haplotype in error block
        γ_e::Array{R, 1}
        # γ_e,  # beta hyper parameter for false positive mutation rate.

        p_merge::R
        p_hap::R
        p_back::R
        p_unique::R
        ln_p_v::Array{R, 1} # error block state penalty, [prob_valid, prob_invalid]
    end

    function parse_config_file(config_file::String)::Tuple{Parameters{INT, REAL}, Annealer{INT, REAL}}
        conf = ConfParse(config_file)
        parse_conf!(conf)
        println(conf)
        α_s   = parse(REAL, String(retrieve(conf, "model", "alpha_s")) )
        α_v   = parse(REAL, String(retrieve(conf, "model", "alpha_v")) )

        δ_s   = parse(REAL, String(retrieve(conf, "model", "delta_s")) )
        δ_m   = parse(REAL, String(retrieve(conf, "model", "delta_m")) )
        δ_u   = parse(REAL, String(retrieve(conf, "model", "delta_u")) )

        λ_smu_0 = parse(REAL, String(retrieve(conf, "model", "lambda_smu_0")) )
        λ_smu_1 = parse(REAL, String(retrieve(conf, "model", "lambda_smu_1")) )
        λ_e_0   = parse(REAL, String(retrieve(conf, "model", "lambda_e_0")) )
        λ_e_1   = parse(REAL, String(retrieve(conf, "model", "lambda_e_1")) )

        β_smu = parse(REAL, String(retrieve(conf, "model", "beta_smu")) )
        β_e   = parse(REAL, String(retrieve(conf, "model", "beta_e")) )

        γ_e_0 = parse(REAL, String(retrieve(conf, "model", "gamma_e_0")) )
        γ_e_1 = parse(REAL, String(retrieve(conf, "model", "gamma_e_1")) )

        p_merge   = parse(REAL, String(retrieve(conf, "model", "p_merge")) )
        p_hap     = parse(REAL, String(retrieve(conf, "model", "p_hap")) )
        p_back    = parse(REAL, String(retrieve(conf, "model", "p_back")) )
        p_unique  = parse(REAL, String(retrieve(conf, "model", "p_unique")) )

        ln_1m_p_v = parse(REAL, String(retrieve(conf, "model", "ln_1m_p_v")) )
        ln_p_v    = log(ℯ, (REAL)(1.0) - exp(ln_1m_p_v))

        period    = parse(REAL, String(retrieve(conf, "annealer", "period")) )

        ln_p_ladders_strings = (retrieve(conf, "annealer", "ln_p_ladders"))
        ln_p_ladders = zeros(REAL, length(ln_p_ladders_strings), 2)

        for i in (INT)(1):(INT)(length(ln_p_ladders_strings))
            step = String(ln_p_ladders_strings[i])
            print("step: "); println(step)
            ln_1m_p_step = parse(REAL, step)
            ln_p_step    = log(ℯ, (REAL)(1.0) - exp(ln_1m_p_step))
            ln_p_ladders[i, 1] = ln_p_step
            ln_p_ladders[i, 2] = ln_1m_p_step
        end

        return (Parameters{INT,REAL}(α_s, α_v,
                                    [δ_s, δ_m, δ_u],
                                    [λ_smu_0 λ_smu_1; λ_e_0 λ_e_1;],
                                    [β_smu, β_e], [γ_e_0, γ_e_1],
                                    p_merge, p_hap, p_back, p_unique,
                                    [ln_p_v, ln_1m_p_v]),
               Annealer{INT, REAL}(period, ln_p_ladders))
    end

    function parse_input_summary(summary_path::String)::Array{REAL, 2}
        ans = []
        isFirst = true
        open(summary_path, "r") do f
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
