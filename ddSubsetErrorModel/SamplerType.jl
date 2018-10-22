@Include "InputParser.jl"
@Include "RandomUtil.jl"
@Include "BuffPhyloMatrixType.jl"
@Include "TableGraphType.jl"
@Include "DistanceParser.jl"

"""
ddSubsetErrorModel sampling script
"""
module SamplerType
    using ..InputParser
    using ..RandomUtil
    using ..BuffPhyloMatrixType
    using ..TableGraphType
    using ..DistanceParser
    using Random
    using Distributions
    using SpecialFunctions


    # DONE: L, S, usage_s ->  L^(1), L^(2), S^(1)_s S^(2)_s, usage_s^(1), usage_s^(2)
    mutable struct Sampler{I,R}
        Z::Array{I, 2} # Z == 1 err, Z == 2 mutation
        H::Array{I, 2} # H == 1 mat, H == 2 pat

        L_s::TableGraph{I} # Table Link graph for sample subset
        L_n::TableGraph{I} # Table Link graph for each identical sample

        s_s::Array{I, 1}
        s_n::Array{I, 1}
        s_v::Array{I, 1}
        usage_s::Dict{I, Array{I, 1}}
        usage_n::Dict{I, Array{I, 1}}
        usage_v::Dict{I, Array{I, 1}}
        unused_s::Set{I}
        unused_n::Set{I}
        unused_v::Set{I}

        a::Array{I, 1}  # a == 1 not merged, a == 2 merged
        f::Array{R, 1}  # freq of mutation for each col
        g::Array{I, 1}  # haplotype for each col
        u::Array{I, 1}  # u[j] : an unique sample having mutation j

        p_err::R        # A frequency in which mutation is false positive.
        er::Array{I, 1} # er[j] == 1 : not error, er[j] == 2 : error
        erset::Set{I}   # {j | er[j] == 2}

        B::Dict{Tuple{I,I}, I}
        # B[c,m] == 1 : shared
        # B[c,m] == 2 : merged
        # B[c,m] == 3 : unique
        # B[c,m] == 4 : error

        tree_cache::BuffPhyloMatrix{I}
        ln_p_data::Array{R, 3}
        param::Parameters{I, R}
        anneal::Annealer{I, R}
        ln_prob::R
    end
    export Sampler

    function isvalid(samp::Sampler{I, R}; debug::Bool = false)::Bool where {I <: Integer, R <: Real}
        return  _is_main_valid(samp.usage_s, samp.usage_v, samp.B) &&
                _is_error_valid(samp.B, samp.usage_s, samp.usage_v, samp.erset, samp.tree_cache, blocktype = (I)(1))
    end

    # DONE: L_n, s_n, usage_n, unused_n, valid checker for debug.
    function assert_nested_cluster_validity(samp::Sampler{I, R}, comment::String = "")::Nothing where {I <: Integer, R <: Real}
        for c in keys(samp.usage_n)
            for i in samp.usage_n[c]
                if samp.s_n[i] != c
                    println("invalid cluster setting now")
                    println(comment)
                    println(samp.L_n)
                    println(samp.usage_n)
                    println(samp.unused_n)
                    println(samp.s_n)
                    @assert samp.s_n[i] == c
                end
            end
        end
        for c in samp.unused_n
            if c ∈ keys(samp.usage_n)
                println("invalid cluster setting now")
                println(comment)
                println(samp.L_n)
                println(samp.usage_n)
                println(samp.unused_n)
                println(samp.s_n)
                @assert c ∉ keys(samp.usage_n)
            end
        end
    end

    #= DONE:
        add:    ln_P_Link(samp.L^(2), samp.L^(1)) in _and ln_p_all
        modify: ln_P_Y to model L^(2)
    =#
    function ln_p_all(samp::Sampler{I,R}, debug::Bool = false)::R where {I<:Integer, R<:Real}
        ans::R = (R)(0.0)
        S::I, M::I = size(samp.Z)
        ans += _ln_p_Link(samp.L_s)
        ans += _ln_p_Link(samp.L_n, samp.usage_s, samp.s_s)
        ans += _ln_p_beta(samp.p_err, samp.param.γ_e[1], samp.param.γ_e[2])
        ans += _ln_p_er(samp.er, samp.erset, samp.p_err)
        ans += _ln_p_crp(samp.usage_v, samp.param.α_v)
        ans += _ln_p_b(samp.B, samp.usage_s, samp.usage_v, samp.param.δ_s)
        ans += _ln_p_v(samp.Z, samp.usage_s, samp.usage_v, samp.B,
                                samp.erset, samp.tree_cache, samp.param.ln_p_v)
        ans += _ln_p_a(samp.a, samp.param.p_merge)
        ans += _ln_p_f(samp.f, samp.param.λ_s, samp.er)
        ans += _ln_p_g(samp.g, samp.param.β_s, samp.er)


        ans += _ln_p_y(samp.Z, samp.H,
                       samp.s_s, samp.s_n ,samp.s_v, samp.usage_n, samp.B,
                       samp.a, samp.f, samp.g, samp.u, samp.er, samp.param)

        ans += _ln_p_data_y(samp.ln_p_data, samp.Z, samp.H)
        if debug
            print("_ln_p_Link(samp.L): ");println(_ln_p_Link(samp.L_s))
            print("_ln_p_er(samp.er, samp.erset, samp.p_err): ");println(_ln_p_er(samp.er, samp.erset, samp.p_err))
            print("_ln_p_crp(samp.usage_v, samp.param.α_v): ");println(_ln_p_crp(samp.usage_v, samp.param.α_v))
            print("_ln_p_b(samp.B, samp.usage_s, samp.usage_v, samp.param.δ_s): ");println(_ln_p_b(samp.B, samp.usage_s, samp.usage_v, samp.param.δ_s))
            print("_ln_p_v(samp.Z, samp.usage_s, samp.usage_v, samp.B, samp.erset, samp.tree_cache, samp.param.ln_p_v): ");println(_ln_p_v(samp.Z, samp.usage_s, samp.usage_v, samp.B, samp.erset, samp.tree_cache, samp.param.ln_p_v))
            print("_ln_p_a(samp.a, samp.param.p_merge): ");println(_ln_p_a(samp.a, samp.param.p_merge))
            print("_ln_p_f(samp.f, samp.param.λ_s, samp.er): ");println(_ln_p_f(samp.f, samp.param.λ_s, samp.er))
            print("_ln_p_g(samp.g, samp.param.β_s, samp.er): ");println(_ln_p_g(samp.g, samp.param.β_s, samp.er))
            print("_ln_p_y(samp.Z, samp.H, samp.s_s, samp.s_v, samp.B, samp.a, samp.f, samp.g, samp.u, samp.er, samp.param): ");println(_ln_p_y(samp.Z, samp.H, samp.s_s, samp.s_v, samp.B, samp.a, samp.f, samp.g, samp.u, samp.er, samp.param))
            print("_ln_p_data_y(samp.ln_p_data, samp.Z, samp.H): ");println(_ln_p_data_y(samp.ln_p_data, samp.Z, samp.H))
        end
        return ans
    end

    #= DONE:
        samp.u[j] == i => samp.u[j] ∈ usage_n[s_n[i]]
        samp.u[j] != i => samp.u[j] ∉ usage_n[s_n[i]]
      TODO: refactoring all sampler funciton to return ratio of ln_p_next / ln_p_prev
    =#
    function _sample_z!(samp::Sampler{I, R}, debug::Bool = false)::Nothing where {I<:Integer, R <: Real}
        S::I, M::I = size(samp.Z)
        t::I = (I)(0)
        ln_p::Array{R,1} = [0.0, 0.0]
        ln_f::Array{R,1} = [0.0, 0.0]
        now_z::I = (I)(0)

        # (debug) && (debug = debug && (R == Float64))
        # (debug) && (ln_p_prev_true::R = (R)(0.0); ln_p_next_true::R = (R)(0.0); ln_p_prev::R = (R)(0.0); ln_p_next::R = (R)(0.0);
        #             ln_p_ratio_true::R = (R)(0.0); ln_p_ratio::R = (R)(0.0); eps::R =(R)(0.0))
        # (debug) && ( if (R==Float64); eps = (R)(1e-2); end;)
        for j in (I)(1):M
            for i in (I)(1):S
                # (debug) && (ln_p_prev_true = ln_p_all(samp))
                now_z = samp.Z[i,j]
                t = 1 + samp.H[i,j]
                ln_p[1] = samp.ln_p_data[i,j,1]; ln_p[2] = samp.ln_p_data[i,j,t];
                ln_f[1] = log(ℯ, 1.0 - samp.f[j]); ln_f[2] = log(ℯ, samp.f[j]);
                if     samp.er[j] == 1 && samp.B[(samp.s_s[i], samp.s_v[j])] == 1
                    ln_p[1] += ln_f[1]
                    ln_p[2] += ln_f[2]
                elseif samp.er[j] == 1 && samp.B[(samp.s_s[i], samp.s_v[j])] == 2
                    ln_p[1] += (samp.a[i] == 2) * ln_f[1] + (samp.a[i] == 1) * (R)(log(ℯ, 1.0-samp.param.p_back))
                    ln_p[2] += (samp.a[i] == 2) * ln_f[2] + (samp.a[i] == 1) * (R)(log(ℯ, samp.param.p_back))
                elseif samp.er[j] == 1 && samp.B[(samp.s_s[i], samp.s_v[j])] == 3
                    isIn::Bool    = (i ∈ samp.usage_n[samp.s_n[samp.u[j]]])
                    isNotIn::Bool = (i ∉ samp.usage_n[samp.s_n[samp.u[j]]])
                    ln_p[1] += (isIn) * ln_f[1] + (isNotIn) * (R)(log(ℯ, 1.0-samp.param.p_unique))
                    ln_p[2] += (isIn) * ln_f[2] + (isNotIn) * (R)(log(ℯ, samp.param.p_unique))
                elseif samp.er[j] == 2
                    ln_p[1] += ln_f[1]
                    ln_p[2] += ln_f[2]
                    isv::Bool        = isvalid(samp)
                    ln_p[now_z]     += (R)(isv) * samp.param.ln_p_v[1] + (R)(!isv) * samp.param.ln_p_v[2]
                    BuffPhyloMatrixType.edit!(samp.tree_cache, i, j, (I)((3 - now_z) - 1))
                    isv              = isvalid(samp)
                    ln_p[3 - now_z] += (R)(isv) * samp.param.ln_p_v[1] + (R)(!isv) * samp.param.ln_p_v[2]
                end
                _exp_normalize!(ln_p)
                samp.Z[i,j] = (I)(argmax( RandomUtil.sample_multinomial((I)(1), ln_p) ))
                if samp.er[j] == 2
                    BuffPhyloMatrixType.edit!(samp.tree_cache, i, j, samp.Z[i,j]-(I)(1))
                end
                # (debug) && (ln_p_next_true = ln_p_all(samp))
                # (debug) && (ln_p_prev = log(ℯ, ln_p[now_z]))
                # (debug) && (ln_p_next = log(ℯ, ln_p[samp.Z[i,j]]))
                # (debug) && (ln_p_ratio_true = ln_p_next_true - ln_p_prev_true;
                #             ln_p_ratio      = ln_p_next      - ln_p_prev;)
                # (debug) && ( if abs(ln_p_ratio_true) > 1e-5 && abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) > eps;
                #                 println((ln_p_ratio_true, ln_p_ratio));
                #                 println((ln_p_next_true, ln_p_next));
                #                 println((ln_p_prev_true, ln_p_prev));
                #                 @assert abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) <= eps;
                #             end;
                # )
            end
        end
        return nothing
    end

    function _sample_h!(samp::Sampler{I, R}, debug::Bool = false)::Nothing where {I<:Integer, R <: Real}
        S::I, M::I = size(samp.Z)
        ln_p::Array{R, 1}   = [0.0, 0.0]
        ln_p_hap::Array{R, 1} = log.(ℯ, [ 1.0 - samp.param.p_hap, samp.param.p_hap ])
        ln_p_rn::Array{R, 1}  = convert.(R, log.(ℯ, [0.50, 0.50])) #convert.(R, log.(ℯ, [0.50, 0.50]))
        # (debug) && (debug = debug && (R == Float64))
        # (debug) && (ln_p_prev_true::R = (R)(0.0); ln_p_next_true::R = (R)(0.0); ln_p_prev::R = (R)(0.0); ln_p_next::R = (R)(0.0);
        #             ln_p_ratio_true::R = (R)(0.0); ln_p_ratio::R = (R)(0.0); eps::R = (R)(0.0);)
        # (debug) && ( if (R==Float64); eps = (R)(1e-2); end;)
        for j in (I)(1):M
            for i in (I)(1):S
                # (debug) && (ln_p_prev_true = ln_p_all(samp))
                # (debug) && (now_H::I = samp.H[i,j])
                ln_p[1] = (samp.Z[i,j] == 2) * samp.ln_p_data[i,j, 1+1]
                ln_p[2] = (samp.Z[i,j] == 2) * samp.ln_p_data[i,j, 1+2]
                if     samp.er[j] == 1
                    ln_p[1] += (samp.g[j] == 1) * ln_p_hap[2] + (samp.g[j] == 2) * ln_p_hap[1]
                    ln_p[2] += (samp.g[j] == 2) * ln_p_hap[2] + (samp.g[j] == 1) * ln_p_hap[1]
                elseif samp.er[j] == 2
                    ln_p[1] += (samp.g[j] == 1) * ln_p_rn[2] + (samp.g[j] == 2) * ln_p_rn[1]
                    ln_p[2] += (samp.g[j] == 2) * ln_p_rn[2] + (samp.g[j] == 1) * ln_p_rn[1]
                end
                _exp_normalize!(ln_p)
                samp.H[i,j] = argmax( RandomUtil.sample_multinomial((I)(1), ln_p) )

                # (debug) && (ln_p_next_true = ln_p_all(samp))
                # (debug) && (ln_p_prev = log(ℯ, ln_p[now_H]))
                # (debug) && (ln_p_next = log(ℯ, ln_p[samp.H[i,j]]))
                # (debug) && (ln_p_ratio_true = ln_p_next_true - ln_p_prev_true;
                #             ln_p_ratio      = ln_p_next      - ln_p_prev;)
                # (debug) && ( if abs(ln_p_ratio_true) > 1e-5 && abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) > eps;
                #                 println((ln_p_ratio_true, ln_p_ratio));
                #                 println((ln_p_next_true, ln_p_next));
                #                 println((ln_p_prev_true, ln_p_prev));
                #                 @assert abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) <= eps;
                #             end;
                # )
            end
        end
        return nothing
    end

    function _sample_a!(samp::Sampler{I, R}, debug::Bool = false)::Nothing where {I <: Integer, R <: Real}
        ln_p::Array{R, 1}   = convert.(R, [0.0, 0.0])
        ln_p_merge::Array{R, 1} = convert.(R, log.(ℯ, [1.0 - samp.param.p_merge, samp.param.p_merge] ))
        ln_p_back::Array{R, 1}  = convert.(R, log.(ℯ, [1.0 - samp.param.p_back, samp.param.p_back] ))
        S::I, M::I = size(samp.Z)
        c::I = (I)(0)
        # (debug) && (debug = debug && (R == Float64))
        # (debug) && (ln_p_prev_true::R = (R)(0.0); ln_p_next_true::R = (R)(0.0); ln_p_prev::R = (R)(0.0); ln_p_next::R = (R)(0.0);
        #             ln_p_ratio_true::R = (R)(0.0); ln_p_ratio::R = (R)(0.0); eps::R = (R)(0.0);)
        # (debug) && ( if (R==Float64); eps = (R)(1e-2); end;)
        for i in (I)(1):S
            # (debug) && (ln_p_prev_true = ln_p_all(samp))
            # (debug) && (now_A::I = samp.a[i])
            ln_p .= ln_p_merge
            c = samp.s_s[i]
            for m in keys(samp.usage_v)
                if samp.B[c,m] == 2 # block merge
                    for j in samp.usage_v[m]
                        ln_p[2] += (samp.Z[i,j] == 2) * log(ℯ, samp.f[j])
                        ln_p[2] += (samp.Z[i,j] == 1) * log(ℯ, 1.0 - samp.f[j])
                        ln_p[1] += (samp.Z[i,j] == 2) * ln_p_back[2]
                        ln_p[1] += (samp.Z[i,j] == 1) * ln_p_back[1]
                    end
                end
            end
            _exp_normalize!(ln_p)
            samp.a[i] = argmax( RandomUtil.sample_multinomial((I)(1), ln_p) )
            # (debug) && (ln_p_next_true = ln_p_all(samp))
            # (debug) && (ln_p_prev = log(ℯ, ln_p[now_A]))
            # (debug) && (ln_p_next = log(ℯ, ln_p[samp.a[i]]))
            # (debug) && (ln_p_ratio_true = ln_p_next_true - ln_p_prev_true;
            #             ln_p_ratio      = ln_p_next      - ln_p_prev;)
            # (debug) && ( if abs(ln_p_ratio_true) > 1e-5 && abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) > eps;
            #                 println((ln_p_ratio_true, ln_p_ratio));
            #                 println((ln_p_next_true, ln_p_next));
            #                 println((ln_p_prev_true, ln_p_prev));
            #                 @assert abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) <= eps;
            #             end;
            # )
        end
        return nothing
    end

    #= DONE:
        samp.u[j] == i => samp.u[j] ∈ usage_n[s_n[i]]
       # TODO: refactoring all sampler funciton to return ratio of ln_p_next / ln_p_prev
    =#
    function _sample_f!(samp::Sampler{I, R}, debug::Bool = false)::Nothing where {I <: Integer, R <: Real}
        m::I = (I)(0)
        λ::Array{R ,1} = [0.0, 0.0]
        S::I, M::I = size(samp.Z)
        # (debug) && (debug = debug && (R == Float64))
        # (debug) && (ln_p_prev_true::R = (R)(0.0); ln_p_next_true::R = (R)(0.0); ln_p_prev::R = (R)(0.0); ln_p_next::R = (R)(0.0);
        #             ln_p_ratio_true::R = (R)(0.0); ln_p_ratio::R = (R)(0.0); eps::R = (R)(0.0);)
        # (debug) && ( if (R==Float64); eps = (R)(1e-2); end;)

        for j in (I)(1):M
            # (debug) && (ln_p_prev_true = ln_p_all(samp))
            # (debug) && (now_f::R = samp.f[j])

            m = samp.s_v[j]
            λ[1] = (samp.er[j] == 1) * samp.param.λ_s[1, 1] + (samp.er[j] == 2) * samp.param.λ_s[2, 1]
            λ[2] = (samp.er[j] == 1) * samp.param.λ_s[1, 2] + (samp.er[j] == 2) * samp.param.λ_s[2, 2]
            for c in keys(samp.usage_s)
                for i in samp.usage_s[c]
                    isadd::Bool = (samp.B[c,m] == 1 || samp.B[c,m] == 4) ||
                                  (samp.B[c,m] == 2 && samp.a[i] == 2)   ||
                                  (samp.B[c,m] == 3 && (i ∈ samp.usage_n[samp.s_n[samp.u[j]]]) )
                    λ[1] += (R)(isadd) * (samp.Z[i,j] == 2)
                    λ[2] += (R)(isadd) * (samp.Z[i,j] == 1)
                end
            end
            samp.f[j] = RandomUtil.sample_beta(λ[1], λ[2])

            # (debug) && (ln_p_next_true = ln_p_all(samp))
            # (debug) && (ln_p_prev = _ln_p_beta(now_f, λ[1], λ[2]))
            # (debug) && (ln_p_next = _ln_p_beta(samp.f[j], λ[1], λ[2]))
            # (debug) && (ln_p_ratio_true = ln_p_next_true - ln_p_prev_true;
            #             ln_p_ratio      = ln_p_next      - ln_p_prev;)
            # (debug) && ( if abs(ln_p_ratio_true) > 1e-5 && abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) > eps;
            #                 println((ln_p_ratio_true, ln_p_ratio));
            #                 println((ln_p_next_true, ln_p_next));
            #                 println((ln_p_prev_true, ln_p_prev));
            #                 @assert abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) <= eps;
            #             end;
            # )
        end
        return nothing
    end

    function _sample_g!(samp::Sampler{I, R}, debug::Bool = false)::Nothing where {I <: Integer, R <: Real}
        S::I, M::I = size(samp.Z)
        m::I = (I)(0)
        ln_p::Array{R, 1} = convert.(R,log.(ℯ, [0.5, 0.5]))
        ln_p_hap::Array{R, 1} = log.(ℯ, [ 1.0 - samp.param.p_hap, samp.param.p_hap ])
        # (debug) && (debug = debug && (R == Float64))
        # (debug) && (ln_p_prev_true::R = (R)(0.0); ln_p_next_true::R = (R)(0.0); ln_p_prev::R = (R)(0.0); ln_p_next::R = (R)(0.0);
        #             ln_p_ratio_true::R = (R)(0.0); ln_p_ratio::R = (R)(0.0); eps::R = (R)(0.0);)
        # (debug) && ( if (R==Float64); eps = (R)(1e-2); end;)
        for j in (I)(1):M
            # (debug) && (ln_p_prev_true = ln_p_all(samp))
            # (debug) && (now_g::I = samp.g[j])

            m = samp.s_v[j]
            ln_p .= samp.param.β_s
            for c in keys(samp.usage_s)
                for i in samp.usage_s[c]
                    (samp.B[c,m] != 4) && (
                        ln_p[1] += ln_p_hap[ 1 + (samp.H[i,j] == 1) ];
                        ln_p[2] += ln_p_hap[ 1 + (samp.H[i,j] == 2) ];
                    )
                end
            end
            # ln_p_temp::Array{R, 1} = deepcopy(ln_p)
            _exp_normalize!(ln_p)
            samp.g[j] = argmax( RandomUtil.sample_multinomial((I)(1), ln_p) )

            # (debug) && (ln_p_next_true = ln_p_all(samp))
            # (debug) && (ln_p_prev = log(ℯ, ln_p[now_g]))
            # (debug) && (ln_p_next = log(ℯ, ln_p[samp.g[j]]))
            # (debug) && (ln_p_ratio_true = ln_p_next_true - ln_p_prev_true;
            #             ln_p_ratio      = ln_p_next      - ln_p_prev;)
            # (debug) && ( if abs(ln_p_ratio_true) > 1e-5 && abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) > eps;
            #                 println((ln_p_ratio_true, ln_p_ratio));
            #                 println((ln_p_next_true, ln_p_next));
            #                 println((ln_p_prev_true, ln_p_prev));
            #                 @assert abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) <= eps;
            #             end;
            # )
        end
        return nothing
    end

    #=
    # TODO: refactoring all sampler funciton to return ratio of ln_p_next / ln_p_prev
    DONE: write a code like this.
    ln_p_cons = 0.0
    for c in filter(c->(samp.B[c,m] == 3), keys(samp.usage_s))
        for u in usage_s[c]
            for i in usage_n[s_n[u]]
                consAt::R = (samp.Z[i,j] == 2) * ln_p_unique[2] + (samp.Z[i,j] == 1) * ln_p_unique[1]
                prosAt::R = (samp.Z[i,j] == 2) * ln_p_f[2]      + (samp.Z[i,j] == 1) * ln_p_f[1]
                ln_p[u] += (prosAt-consAt)
                ln_p_cons += consAt
            end
        end
    end
    ln_p .+= ln_p_cons
    =#
    function _sample_u!(samp::Sampler{I, R}, debug::Bool = false)::Nothing where {I <: Integer, R <: Real}
        S::I, M::I = size(samp.Z)
        m::I = (I)(0)
        ln_p::Array{R, 1} = zeros(R, S)
        ln_p_unique::Array{R, 1} = convert.(R, log.(ℯ, [1.0 - samp.param.p_unique, samp.param.p_unique]) )
        ln_p_f::Array{R, 1} = [0.0, 0.0]

        # (debug) && (debug = debug && (R == Float64))
        # (debug) && (ln_p_prev_true::R = (R)(0.0); ln_p_next_true::R = (R)(0.0); ln_p_prev::R = (R)(0.0); ln_p_next::R = (R)(0.0);
        #             ln_p_ratio_true::R = (R)(0.0); ln_p_ratio::R = (R)(0.0); eps::R = (R)(0.0);)
        # (debug) && ( if (R==Float64); eps = (R)(1e-2); end;)
        for j in (I)(1):M
            # (debug) && (ln_p_prev_true = ln_p_all(samp))
            # (debug) && (now_u::I = samp.u[j])

            m = samp.s_v[j]
            ln_p .= (R)(0.0)
            ln_p_f[1] = log(ℯ, 1.0 - samp.f[j]); ln_p_f[2] = log(ℯ, samp.f[j]);

            ln_p_cons = (R)(0.0)
            for c in filter(c->(samp.B[c,m] == 3), keys(samp.usage_s))
                for u in samp.usage_s[c]
                    for i in samp.usage_n[samp.s_n[u]]
                        consAt::R = (samp.Z[i,j] == 2) * ln_p_unique[2] + (samp.Z[i,j] == 1) * ln_p_unique[1]
                        prosAt::R = (samp.Z[i,j] == 2) * ln_p_f[2]      + (samp.Z[i,j] == 1) * ln_p_f[1]
                        ln_p[u] += (prosAt-consAt)
                        ln_p_cons += consAt
                    end
                end
            end
            ln_p .+= ln_p_cons

            # ln_p_temp::Array{R, 1} = deepcopy(ln_p)
            _exp_normalize!(ln_p)
            samp.u[j] = argmax( RandomUtil.sample_multinomial((I)(1), ln_p) )

            # (debug) && (ln_p_next_true = ln_p_all(samp))
            # (debug) && (ln_p_prev = log(ℯ, ln_p[now_u]))
            # (debug) && (ln_p_next = log(ℯ, ln_p[samp.u[j]]))
            # (debug) && (ln_p_ratio_true = ln_p_next_true - ln_p_prev_true;
            #             ln_p_ratio      = ln_p_next      - ln_p_prev;)
            # (debug) && ( if abs(ln_p_ratio_true) > 1e-5 && abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) > eps;
            #                 println((ln_p_ratio_true, ln_p_ratio));
            #                 println((ln_p_next_true, ln_p_next));
            #                 println((ln_p_prev_true, ln_p_prev));
            #                 @assert abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) <= eps;
            #             end;
            # )
        end
        return nothing
    end

    function _sample_s_s!(samp::Sampler{I, R}, i::I)::Nothing where {I <: Integer, R <: Real}
        c::I = samp.s_s[i]
        prev_cluster::I = samp.s_s[i]
        prev_B::Dict{Tuple{I,I}, I} = deepcopy(samp.B)
        next_cluster::I = _sample_crp(samp.s_s[i], samp.usage_s, samp.unused_s, samp.param.α_s)
        next_B::Dict{Tuple{I,I}, I} = deepcopy(samp.B)

        (   # previous cluster disappears
            if length(samp.usage_s[prev_cluster]) == 1;
                for m in keys(samp.usage_v);
                    pop!(next_B, (prev_cluster, m));
                end;
            end;
            # novel cluster appears
            if ( next_cluster ∈ samp.unused_s ) ||
               ( length(samp.unused_s) == 0  && next_cluster == prev_cluster);
                for m in keys(samp.usage_v);
                    (m != 0) && (next_B[next_cluster, m] = argmax( RandomUtil.sample_multinomial((I)(1), samp.param.δ_s) ));
                    (m == 0) && (next_B[next_cluster, m] = 4);
                end;
            end;
        )
        # acc from prev s_s, B related
        acc_prev::R = (R)(0.0)
        (
            acc_prev += _ln_p_v(samp.Z, samp.usage_s, samp.usage_v, samp.B, samp.erset,
                                samp.tree_cache, samp.param.ln_p_v);
            acc_prev += _ln_p_y(samp.Z, samp.H, samp.s_s, samp.s_v, samp.B, samp.a, samp.f,
                                samp.g, samp.u, samp.er, samp.param, range_s = i:i);
        )

        acc_next::R = (R)(0.0)
        (
            _update_cluster!(prev_cluster, next_cluster, samp.s_s, samp.usage_s, samp.unused_s, i);
            samp.B = deepcopy(next_B);
            BuffPhyloMatrixType.update!(samp.tree_cache, samp.B, samp.usage_s, samp.usage_v);
            acc_next += _ln_p_v(samp.Z, samp.usage_s, samp.usage_v, samp.B, samp.erset,
                                samp.tree_cache, samp.param.ln_p_v);
            acc_next += _ln_p_y(samp.Z, samp.H, samp.s_s, samp.s_v, samp.B, samp.a, samp.f,
                                samp.g, samp.u, samp.er, samp.param, range_s = i:i);
        )
        # select prev/next
        acc_rate::R = min(1.0, exp(acc_next - acc_prev))
        if rand() > acc_rate # rejected
            # revert to a previous state
            _update_cluster!(next_cluster, prev_cluster, samp.s_s, samp.usage_s, samp.unused_s, i)
            samp.B = deepcopy(prev_B)
            BuffPhyloMatrixType.update!(samp.tree_cache, samp.B, samp.usage_s, samp.usage_v)
        end
        return nothing
    end

    function _sample_p_err!(samp::Sampler{I, R}, debug::Bool = false)::Nothing where {I <: Integer, R <: Real}
        # prev_p_err::R = samp.p_err      # debug
        # prev_p_all::R = ln_p_all(samp)  # debug
        # (debug) && (debug = debug && (R == Float64))
        # (debug) && (ln_p_prev_true::R = (R)(0.0); ln_p_next_true::R = (R)(0.0); ln_p_prev::R = (R)(0.0); ln_p_next::R = (R)(0.0);
        #             ln_p_ratio_true::R = (R)(0.0); ln_p_ratio::R = (R)(0.0); eps::R = (R)(0.0);)
        # (debug) && ( if (R==Float64); eps = (R)(1e-2); end;)
        # (debug) && (ln_p_prev_true = ln_p_all(samp))
        # (debug) && (now_P::R = samp.p_err)

        γ_e = [0.0, 0.0]
        γ_e .= samp.param.γ_e
        γ_e[1] += (R)(length(samp.erset))
        γ_e[2] += (R)(length(samp.er) - length(samp.erset))
        samp.p_err = RandomUtil.sample_beta(γ_e[1], γ_e[2])

        # (debug) && (ln_p_next_true = ln_p_all(samp))
        # (debug) && (ln_p_prev = _ln_p_beta(now_P, γ_e[1], γ_e[2]))
        # (debug) && (ln_p_next = _ln_p_beta(samp.p_err, γ_e[1], γ_e[2]))
        # (debug) && (ln_p_ratio_true = ln_p_next_true - ln_p_prev_true;
        #             ln_p_ratio      = ln_p_next      - ln_p_prev;)
        # (debug) && ( if abs(ln_p_ratio_true) > 1e-5 && abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) > eps;
        #                 println((ln_p_ratio_true, ln_p_ratio));
        #                 println((ln_p_next_true, ln_p_next));
        #                 println((ln_p_prev_true, ln_p_prev));
        #                 @assert abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) <= eps;
        #             end;
        # )
        return nothing
    end

    # TODO: refactoring all sampler funciton to return ratio of ln_p_next / ln_p_prev
    # DONE: change argment of ln_P_Y, ln_P_V, and s_s -> S^(1)_s
    function _sample_s_v!(samp::Sampler{I, R}, j::I, debug::Bool = false)::Nothing where {I <: Integer, R <: Real}
        # (debug) && (debug = debug && (R == Float64))
        # (debug) && (ln_p_prev_true::R = (R)(0.0); ln_p_next_true::R = (R)(0.0); ln_q_prev::R = (R)(0.0); ln_q_next::R = (R)(0.0);
        #             ln_p_ratio_true::R = (R)(0.0); ln_p_ratio::R = (R)(0.0); eps::R = (R)(0.0);)
        # (debug) && ( if (R==Float64); eps = (R)(1e-2); end;)
        # (debug) && (ln_p_prev_true = ln_p_all(samp))
        # (debug) && (ln_q_prev += _ln_p_er(samp.er, samp.erset, samp.p_err) )
        # (debug) && (ln_q_prev += _ln_p_crp(samp.usage_v, samp.param.α_v) )
        # (debug) && (ln_q_prev += _ln_p_b(samp.B, samp.usage_s, samp.usage_v, samp.param.δ_s) )
        # (debug) && (ln_q_prev += _ln_p_f(samp.f, samp.param.λ_s, samp.er) )

        prev_cluster::I = samp.s_v[j]
        prev_er::I = samp.er[j]
        prev_f::R  = samp.f[j]
        prev_B::Dict{Tuple{I,I}, I} = deepcopy(samp.B)

        next_er::I = argmax( RandomUtil.sample_multinomial((I)(1), [1.0 - samp.p_err, samp.p_err]) )
        next_cluster::I = (I)(0)
        (next_er==1) && (next_cluster = _sample_crp(samp.s_v[j], samp.usage_v, samp.unused_v, samp.param.α_v) )
        (next_er==2) && (next_cluster = 0)
        next_f::R  = RandomUtil.sample_beta( samp.param.λ_s[next_er, 1], samp.param.λ_s[next_er, 2] )
        next_B::Dict{Tuple{I,I}, I} = deepcopy(samp.B)
        (   # init next_B
            # previous cluster disappears
            if length(samp.usage_v[prev_cluster]) == 1;
                for c in keys(samp.usage_s);
                    pop!(next_B, (c, prev_cluster));
                end;
            end;
            if next_cluster == 0;
                for c in keys(samp.usage_s);
                    next_B[c, next_cluster] = 4;
                end;
            elseif next_cluster != 0; # novel cluster appears in non error cluster
                if ( next_cluster ∈ samp.unused_v ) ||
                   ( length(samp.unused_v) == 0  && next_cluster == prev_cluster);
                    for c in keys(samp.usage_s);
                        next_B[c, next_cluster] = argmax( RandomUtil.sample_multinomial((I)(1), samp.param.δ_s) );
                    end;
                end;
            end;
        )
        acc_prev::R = (R)(0.0)  # acc from prev s_s, B related
        (
            acc_prev += _ln_p_v(samp.Z, samp.usage_s, samp.usage_v, samp.B, samp.erset,
                                samp.tree_cache, samp.param.ln_p_v);
            acc_prev += _ln_p_y(samp.Z, samp.H, samp.s_s, samp.s_n, samp.s_v, samp.usage_n, samp.B, samp.a, samp.f,
                                samp.g, samp.u, samp.er, samp.param, range_v = j:j);
        )
        acc_next::R = (R)(0.0) # update prev to next state
        (
            _update_subset!(j, prev_er, next_er, samp.erset, samp.er);
            _update_cache_in_error!(j, prev_er, next_er, view(samp.Z, :, j).-(I)(1), samp.tree_cache);
            _update_cluster!(prev_cluster, next_cluster, samp.s_v, samp.usage_v, samp.unused_v, j);
            samp.B = deepcopy(next_B);
            samp.f[j] = next_f;
            BuffPhyloMatrixType.update!(samp.tree_cache, samp.B, samp.usage_s, samp.usage_v);
            acc_next += _ln_p_v(samp.Z, samp.usage_s, samp.usage_v, samp.B, samp.erset,
                                samp.tree_cache, samp.param.ln_p_v);
            acc_next += _ln_p_y(samp.Z, samp.H, samp.s_s, samp.s_n, samp.s_v, samp.usage_n, samp.B, samp.a, samp.f,
                                samp.g, samp.u, samp.er, samp.param, range_v = j:j);
        )
        # (debug) && (ln_p_next_true = ln_p_all(samp))
        # (debug) && (ln_q_next += _ln_p_er(samp.er, samp.erset, samp.p_err) )
        # (debug) && (ln_q_next += _ln_p_crp(samp.usage_v, samp.param.α_v) )
        # (debug) && (ln_q_next += _ln_p_b(samp.B, samp.usage_s, samp.usage_v, samp.param.δ_s) )
        # (debug) && (ln_q_next += _ln_p_f(samp.f, samp.param.λ_s, samp.er) )
        #
        # (debug) && (ln_p_ratio_true = ln_p_next_true - ln_p_prev_true;
        #             ln_p_ratio      = acc_next - acc_prev - ln_q_prev + ln_q_next;)
        # (debug) && ( if abs(ln_p_ratio_true) > 1e-5 && abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) > eps;
        #                 println(j);
        #                 println((ln_p_ratio_true, ln_p_ratio));
        #                 println((ln_p_next_true, ln_q_next, acc_next));
        #                 println((ln_p_prev_true, ln_q_prev, acc_prev));
        #                 @assert abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) <= eps;
        #             end;
        # )
        # select prev/next
        acc_rate::R = min((R)(1.0), exp(acc_next - acc_prev))
        if rand() > acc_rate # rejected and revert to a previous state
            _update_subset!(j, next_er, prev_er, samp.erset, samp.er)
            _update_cache_in_error!(j, next_er, prev_er, view(samp.Z, :, j).-(I)(1), samp.tree_cache)
            _update_cluster!(next_cluster, prev_cluster, samp.s_v, samp.usage_v, samp.unused_v, j)
            samp.B = deepcopy(prev_B);
            samp.f[j]  = prev_f
            BuffPhyloMatrixType.update!(samp.tree_cache, samp.B, samp.usage_s, samp.usage_v)
        end
        return nothing
    end

    # DONE: Check again whether your M-H algorithm update is correct.
    #=
    DONE: ln_P_Y: add argument
    =#
    function _ln_p_y_marginal!(samp::Sampler{I, R}, c::I)::R where {I <: Integer, R <: Real}
        ans::R = (R)(0.0)
        ln_p::Array{R, 1} = [0.0, 0.0, 0.0]
        for m in keys(samp.usage_v)
            ln_p .= (R)(0.0)
            now_b::I = samp.B[c,m]
            if now_b == 4
                continue
            end
            for b in (I)(1):(I)(3)
                ln_p[b] += log(ℯ, samp.param.δ_s[b])
                samp.B[c,m] = b
                ln_p[b] += _ln_p_y(samp.Z, samp.H, samp.s_s, samp.s_n, samp.s_v, samp.usage_n, samp.B, samp.a, samp.f,
                                   samp.g, samp.u, samp.er, samp.param,
                                   range_s = samp.usage_s[c], range_v = samp.usage_v[m])
            end
            samp.B[c,m] = now_b
            ans += _log_sum_exp(ln_p)
        end
        return ans
    end

    function _ln_p_b_y!(samp::Sampler{I, R})::R where {I <: Integer, R <: Real}
        ans::R = (R)(0.0)
        ln_p::Array{R, 1} = convert.(R, [0.0, 0.0, 0.0])
        for c in keys(samp.usage_s)
            for m in keys(samp.usage_v)
                now_b::I = samp.B[c, m]
                if (samp.B[c, m] == 4)
                    continue
                end
                ln_p .= (R)(0.0)
                for t in (I)(1):(I)(3)
                    samp.B[c, m] = t
                    ln_p[t] += log(ℯ, samp.param.δ_s[t])
                    ln_p[t] += _ln_p_y(samp.Z, samp.H, samp.s_s, samp.s_n, samp.s_v, samp.usage_n, samp.B, samp.a, samp.f,
                                       samp.g, samp.u, samp.er, samp.param,
                                       range_s = samp.usage_s[c], range_v = samp.usage_v[m])
                end
                samp.B[c, m] = now_b
                _exp_normalize!(ln_p)
                ans += log(ℯ, ln_p[now_b])
            end
        end
        return ans
    end

    # DONE: samp.usage_s -> samp.usage_s^(1), samp.s_s -> samp.S^(1)_s
    # DONE: add argument to ln_P_Y
    function _sample_b_patially_gibbs!(samp::Sampler{I,R}, c::I, m::I)::Nothing where {I <: Integer, R <: Real }
        (samp.B[c, m] == 4) && (return nothing)
        ln_p::Array{R, 1} = convert.(R, [0.0, 0.0, 0.0])
        for t in (I)(1):(I)(3)
            samp.B[c, m] = t
            ln_p[t] += log(ℯ, samp.param.δ_s[t])
            ln_p[t] += _ln_p_y(samp.Z, samp.H, samp.s_s, samp.s_n, samp.s_v, samp.usage_n, samp.B, samp.a, samp.f,
                               samp.g, samp.u, samp.er, samp.param,
                               range_s = samp.usage_s[c], range_v = samp.usage_v[m])
        end
        # ln_p_temp::Array{R, 1} = deepcopy(ln_p)
        _exp_normalize!(ln_p)
        samp.B[c, m] = argmax( RandomUtil.sample_multinomial((I)(1), ln_p) )
        BuffPhyloMatrixType.update!(samp.tree_cache, samp.B, samp.usage_s, samp.usage_v)
        return nothing
    end

    #=
    DONE: add proposing L_n part
        add: _propose_nested_link(smap.L.W, i, usage_s, s_s)::tableGraph{I, R}
        add: _set_cluster_by_link(L_n, s_n, usage_n, unused_n)
        modify: update L_n state
    # TODO: refactoring all sampler funciton to return ratio of ln_p_next / ln_p_prev
    =#
    function _sample_l!(samp::Sampler{I, R}, i::I, debug::Bool = false)::Nothing where {I <: Integer, R <: Real}
        # (debug) && (debug = debug && (R == Float64))
        # (debug) && (ln_p_prev_true::R = (R)(0.0); ln_p_next_true::R = (R)(0.0); ln_q_prev::R = (R)(0.0); ln_q_next::R = (R)(0.0);
        #             ln_p_ratio_true::R = (R)(0.0); ln_p_ratio::R = (R)(0.0); eps::R = (R)(0.0);)
        # (debug) && ( if (R==Float64); eps = (R)(1e-2); end;)
        # (debug) && (ln_p_prev_true = ln_p_all(samp))
        # (debug) && (ln_q_prev += _ln_p_Link(samp.L_s) )
        # (debug) && (ln_q_prev += _ln_p_Link(samp.L_n, samp.usage_s, samp.s_s) )
        # (debug) && (ln_q_prev += _ln_p_b_y!(samp) )

        prev_to::I = TableGraphType.getlink(samp.L_s, i)
        next_to::I = _sample_link(samp.L_s.W, i)
        s_pp::Set{I} = Set{I}([])
        s_pn::Set{I} = Set{I}([])
        s_np::Set{I} = Set{I}([])
        s_nn::Set{I} = Set{I}([])
        TableGraphType.diffgroup!(samp.L_s, s_pp, s_pn, s_np, s_nn, i, prev_to, next_to, both = true)
        prev_L_n::TableGraph{I, R} = deepcopy(samp.L_n)

        prev_usage_s::Dict{I, Array{I, 1}} = deepcopy(samp.usage_s)
        prev_unused_s::Set{I}              = deepcopy(samp.unused_s)
        prevs_s::Array{I}                  = deepcopy(samp.s_s)
        prev_usage_n::Dict{I, Array{I, 1}} = deepcopy(samp.usage_n)
        prev_unused_n::Set{I}              = deepcopy(samp.unused_n)
        prevs_n::Array{I}                  = deepcopy(samp.s_n)
        prev_B::Dict{Tuple{I,I}, I}        = deepcopy(samp.B)
        acc_prev::R = (R)(0.0)
        (
            acc_prev += _ln_p_y_marginal!(samp, samp.s_s[prev_to]);
            if samp.s_s[prev_to] != samp.s_s[next_to];
                acc_prev += _ln_p_y_marginal!(samp, samp.s_s[next_to]);
            end;
            acc_prev += _ln_p_v(samp.Z, samp.usage_s, samp.usage_v, samp.B, samp.erset,
                                samp.tree_cache, samp.param.ln_p_v);
        )
        acc_next::R = (R)(0.0)
        (# update to proposal state
            TableGraphType.rmedge!(samp.L_s, i, prev_to);
            TableGraphType.addedge!(samp.L_s, i, next_to);
            (removed::Array{I, 1}, added::Array{I, 1}) = _update_s_s_by_group_diff!(s_pp, s_pn, s_np, s_nn,
                                                                                        samp.usage_s, samp.unused_s, samp.s_s);
            err_m::Set{I} = Set{I}([]);
            for c in removed; for m in keys(samp.usage_v);
                (samp.B[c,m] == 4) && (push!(err_m, m));
                pop!(samp.B, (c,m));
            end; end;
            for c in added; for m in keys(samp.usage_v);
                samp.B[c,m] = 1;
                (m ∈ err_m) && (samp.B[c,m] = 4);
            end; end;
            _propose_nested_link!(samp.L_n, samp.usage_s, samp.s_s, Set{I}([prev_to, next_to]));
            _set_cluster_by_link!(samp.L_n, samp.s_n, samp.usage_n, samp.unused_n);
            for c in added; for m in keys(samp.usage_v);
                _sample_b_patially_gibbs!(samp, c, m);
            end; end;
            acc_next += _ln_p_y_marginal!(samp, samp.s_s[prev_to]);
            if samp.s_s[prev_to] != samp.s_s[next_to];
                acc_next += _ln_p_y_marginal!(samp, samp.s_s[next_to]);
            end;
            acc_next += _ln_p_v(samp.Z, samp.usage_s, samp.usage_v, samp.B, samp.erset,
                                samp.tree_cache, samp.param.ln_p_v);
        )
        # select prev/next
        acc_rate::R = min((R)(1.0), exp(acc_next - acc_prev))

        # (debug) && (ln_p_next_true = ln_p_all(samp))
        # (debug) && (ln_q_next += _ln_p_Link(samp.L_s) )
        # (debug) && (ln_q_next += _ln_p_Link(samp.L_n, samp.usage_s, samp.s_s) )
        # (debug) && (ln_q_next += _ln_p_b_y!(samp) )
        # (debug) && (ln_p_ratio_true = ln_p_next_true - ln_p_prev_true;
        #             ln_p_ratio      = acc_next - acc_prev - ln_q_prev + ln_q_next;)

        # (debug) && ( if abs(ln_p_ratio_true) > 1e-5 && abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) > eps;
        #                 println("==== NG =====");
        #                 println(i);
        #                 println(samp.usage_s);
        #                 println(prev_usage_s);
        #                 println(samp.usage_n);
        #                 println(prev_usage_n);
        #                 println((ln_p_ratio_true, ln_p_ratio));
        #                 println((ln_p_next_true, ln_q_next, acc_next));
        #                 println((ln_p_prev_true, ln_q_prev, acc_prev));
        #                 @assert abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) <= eps;
        #             end;
        # )

        if rand() > acc_rate # rejected
            TableGraphType.rmedge!(samp.L_s, i, next_to)
            TableGraphType.addedge!(samp.L_s, i, prev_to)
            samp.L_n      = deepcopy(prev_L_n)
            samp.usage_n  = deepcopy(prev_usage_n)
            samp.unused_n = deepcopy(prev_unused_n)
            samp.s_n      = deepcopy(prevs_n)
            samp.usage_s  = deepcopy(prev_usage_s)
            samp.unused_s = deepcopy(prev_unused_s)
            samp.s_s      = deepcopy(prevs_s)
            samp.B        = deepcopy(prev_B)
        end
        return nothing
    end

    # DONE: ln_P_Y add argument
    # TODO: refactoring all sampler funciton to return ratio of ln_p_next / ln_p_prev
    function _sample_b!(samp::Sampler{I, R}, c::I, m::I, debug::Bool = false) where {I <: Integer, R <: Real}
        (samp.B[c, m] == 4) && (return nothing)
        # ln_p_correct::Array{R, 1} = [0.0, 0.0, 0.0]
        ln_p::Array{R, 1} = convert.(R, [0.0, 0.0, 0.0])
        range_b::AbstractRange{I} = (I)(1):(I)(3)
        (samp.B[c,m] != 1) && (range_b = (I)(3):(I)(-1):(I)(1))

        # (debug) && (debug = debug && (R == Float64))
        # (debug) && (ln_p_prev_true::R = (R)(0.0); ln_p_next_true::R = (R)(0.0); ln_p_prev::R = (R)(0.0); ln_p_next::R = (R)(0.0);
        #             ln_p_ratio_true::R = (R)(0.0); ln_p_ratio::R = (R)(0.0); eps::R = (R)(0.0);)
        # (debug) && ( if (R==Float64); eps = (R)(1e-2); end;)
        # (debug) && (ln_p_prev_true = ln_p_all(samp))
        # (debug) && (now_b::I = samp.B[c,m])

        for t in (I)(1):(I)(3)
            samp.B[c, m] = t
            BuffPhyloMatrixType.update!(samp.tree_cache, samp.B, samp.usage_s, samp.usage_v)
            ln_p[t] += log(ℯ, samp.param.δ_s[t])
            ln_p[t] += _ln_p_y(samp.Z, samp.H, samp.s_s, samp.s_n, samp.s_v, samp.usage_n, samp.B, samp.a, samp.f,
                               samp.g, samp.u, samp.er, samp.param,
                               range_s = samp.usage_s[c], range_v = samp.usage_v[m])
            ln_p[t] += _ln_p_v(samp.Z, samp.usage_s, samp.usage_v, samp.B, samp.erset,
                               samp.tree_cache, samp.param.ln_p_v)
        end
        # ln_p_temp::Array{R, 1} = deepcopy(ln_p)
        _exp_normalize!(ln_p)
        samp.B[c, m] = argmax( RandomUtil.sample_multinomial((I)(1), ln_p) )
        BuffPhyloMatrixType.update!(samp.tree_cache, samp.B, samp.usage_s, samp.usage_v)

        # (debug) && (ln_p_next_true = ln_p_all(samp))
        # (debug) && (ln_p_prev = log(ℯ, ln_p[now_b]))
        # (debug) && (ln_p_next = log(ℯ, ln_p[samp.B[c,m]]))
        # (debug) && (ln_p_ratio_true = ln_p_next_true - ln_p_prev_true;
        #             ln_p_ratio      = ln_p_next      - ln_p_prev;)
        # (debug) && ( if abs(ln_p_ratio_true) > 1e-5 && abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) > eps;
        #                 println((ln_p_ratio_true, ln_p_ratio));
        #                 println((ln_p_next_true, ln_p_next));
        #                 println((ln_p_prev_true, ln_p_prev));
        #                 @assert abs((ln_p_ratio_true - ln_p_ratio)/ln_p_ratio) <= eps;
        #             end;
        # )

        return nothing
    end

    # DONE: _sample_l! -> sampleLs!
    # DONE: modify _sample_s_v!, _sample_b!, _sample_f!, _sample_z!
    # TODO: refactoring all sampler funciton to return ratio of ln_p_next / ln_p_prev
    function sample_map!(samp::Sampler{I, R};
                         seed::I = (I)(0),
                         iter::I = (I)(100000),
                         thin::I = (I)(1),
                         burnin::I = (I)(0),
                         progress_count::I = (I)(100))::Tuple{Sampler{I, R}, Array{R, 1}} where {I <:Integer, R <: Real}
        Random.seed!(seed)
        ln_p_v_true::Array{R, 1} = samp.param.ln_p_v
        map_state::Sampler{I,R} = deepcopy(samp)
        ln_probs::Array{R, 1} = []
        max_ln_prob::R = (R)(-Inf)
        for count in (I)(1):(iter+burnin)
            S::I, M::I = size(samp.Z)
            _update_penalty!((I)(count),
                             samp.param.ln_p_v,
                             samp.anneal.ln_p_ladders,
                             samp.anneal.period)
            _sample_z!(samp)
            _sample_h!(samp)

            for i in (I)(1):S
                _sample_l!(samp, (I)(i))
            end

            _sample_p_err!(samp)
            for j in (I)(1):M
                _sample_s_v!(samp,  (I)(j))
            end


            _sample_a!(samp)
            _sample_f!(samp)
            _sample_g!(samp)
            _sample_u!(samp)

            for (c,m) in Iterators.product(keys(samp.usage_s), keys(samp.usage_v))
                _sample_b!(samp,  (I)(c), (I)(m))
            end

            if count > burnin && count % thin == 0
                samp.param.ln_p_v .= ln_p_v_true
                samp.ln_prob = ln_p_all(samp)
                push!(ln_probs, samp.ln_prob)
                if max_ln_prob < samp.ln_prob
                    max_ln_prob = samp.ln_prob
                    map_state = deepcopy(samp)
                end
            end
            if count % progress_count == 0
                println(count)
            end
        end
        return (map_state, ln_probs)
    end

    # DONE: _sample_l! -> sampleLs!
    # DONE: modify _sample_s_v!, _sample_b!, _sample_f!, _sample_z!
    # return the (state of MAP, iterCount, ln_prob) with in this iterations
    function sample_all!(samp::Sampler{I, R};
                         seed::I = (I)(0),
                         iter::I = (I)(5000),
                         thin::I = (I)(10),
                         burnin::I = (I)(1000))::Array{Tuple{Sampler{I, R}, I}, 1} where {I <:Integer, R <: Real }
        # setting the given random seed
        Random.seed!(seed)
        ln_p_v_true::Array{R, 1} = samp.param.ln_p_v
        sampled::Array{ Tuple{Sampler{I,R}, I}, 1} = []
        for count in (I)(1):(iter+burnin)
            S::I, M::I = size(samp.Z)
            _update_penalty!((I)(count),
                             samp.param.ln_p_v,
                             samp.anneal.ln_p_ladders,
                             samp.anneal.period)
            _sample_z!(samp)
            _sample_h!(samp)

            for i in (I)(1):S
                _sample_l!(samp, i)
            end

            for j in (I)(1):M
                _sample_s_v!(samp, j)
            end

            _sample_a!(samp)
            _sample_f!(samp)
            _sample_g!(samp)
            _sample_u!(samp)

            for (c,m) in Iterators.product(keys(samp.usage_s), keys(samp.usage_v))
                _sample_b!(samp, c, m)
            end

            if count > burnin && count % thin == 0
                samp.param.ln_p_v .= ln_p_v_true
                samp.ln_prob = ln_p_all(samp)
                now::Sampler{I,R} = deepcopy(samp)
                push!(sampled, (now, count))
            end

            if count % 100 == 0
                println(count)
            end
        end
        return sampled
    end

    # DONE: L, S, usage_s ->  L^(1), L^(2), S^(1)_s S^(2)_s, usage_s^(1), usage_s^(2)
    # DONE: Use different distance matix for L_s and L_n
    function init(err_score_path::String, pat_score_path::String, mat_score_path::String, param_path::String,
                  INT::Type{<:Integer} = Int32, REAL::Type{<:Real} = Float32)
        ln_p_d::Array{REAL, 3}   = _parseData(err_score_path, pat_score_path, mat_score_path, INT, REAL)
        param::Parameters{INT,REAL}, anneal::Annealer{INT, REAL} = InputParser.parse_config_file(param_path, INT, REAL)
        println("=========== data matrix ==========")
        println(ln_p_d)
        println("==================================")

        S::INT = size(ln_p_d, 1)
        M::INT = size(ln_p_d, 2)

        Z::Array{INT, 2}   = convert.(INT, fill(1,S,M)) # init ℤ, Z[i,j] ∈ {1,2}, 1: error, 2: tumor
        H::Array{INT, 2}   = convert.(INT, fill(1,S,M)) # init H, H[i,j] ∈ {1,2}, 1: mat,   2: pat

        D_s::Array{REAL, 2}  = DistanceParser.parse_bf_hamming_distance(err_score_path, pat_score_path, mat_score_path, param_path, "alpha_s", INT, REAL)
        decayfunction_s = DistanceParser.parse_decayfunction(param_path, INT = INT, REAL = REAL)
        D_n::Array{REAL, 2}  = DistanceParser.parse_bf_hamming_distance(err_score_path, pat_score_path, mat_score_path, param_path, "alpha_n")
        decayfunction_n = DistanceParser.parse_decayfunction(param_path, distance_tag = "distance",
                                                             decay_function_tag = "decayFunctionNest", decay_rate_tag = "decayRateNest",
                                                             INT = INT, REAL = REAL)
        for (i,j) in Iterators.product((INT)(1):S,(INT)(1):S)
            (i != j) && ( D_s[i,j] = decayfunction_s(D_s[i,j]))
            (i != j) && ( D_n[i,j] = decayfunction_n(D_n[i,j]))
        end
        println("f(distance)")
        println(D_s)
        println(D_n)
        L_s::TableGraph{INT, REAL} = TableGraphType.init(S, D_s)
        L_n::TableGraph{INT, REAL} = TableGraphType.init(S, D_n)

        s_s::Array{INT, 1} = convert.(INT, collect(1:S)) # init sample wise cluster
        s_n::Array{INT, 1} = convert.(INT, collect(1:S)) # init sample wise cluster
        s_v::Array{INT, 1} = convert.(INT, collect(1:M)) # init mutation wise cluster
        usage_s::Dict{INT,Array{INT,1}} = Dict{INT,Array{INT,1}}()
        usage_n::Dict{INT,Array{INT,1}} = Dict{INT,Array{INT,1}}()
        usage_v::Dict{INT,Array{INT,1}} = Dict{INT,Array{INT,1}}()
        for s in (INT)(1):S; usage_s[s] = [s]; end
        for s in (INT)(1):S; usage_n[s] = [s]; end
        for m in (INT)(1):M; usage_v[m] = [m]; end

        a::Array{INT, 1} = convert.(INT, fill(1, S))    # init merge sample indicator
        f::Array{REAL,1} = convert.(REAL,fill(0.90, M)) # init mutation freq for each mutation
        g::Array{INT, 1} = convert.(INT, fill(1, M))    # init haplotype for each mutation
        u::Array{INT, 1} = convert.(INT, fill(1, M))    # init mutation to unique sample indicator
        p_err::REAL      = (REAL)(0.05)                      # init false positive mutation rate
        er::Array{INT,1} = convert.(INT, fill(1, M))    # init error indicator for each mutation

        # init block cluster ∈ {1,2,3,4}, 1:shared, 2:meged, 3:unique, 4:error
        B::Dict{Tuple{INT,INT},INT} = Dict{Tuple{INT,INT},INT}()
        for (s,m) in Iterators.product((INT)(1):S,(INT)(1):M); B[(s,m)] = 1; end

        tree_cache::BuffPhyloMatrix{INT} = BuffPhyloMatrixType.init(S, M, buffer_size = M)

        samp::Sampler{INT, REAL} = Sampler{INT, REAL}(Z, H, L_s, L_n, s_s, s_n, s_v, usage_s, usage_n, usage_v, Set{INT}(), Set{INT}(), Set{INT}(),
                                                      a, f, g, u, p_err, er, Set{INT}(),  B, tree_cache, ln_p_d, param, anneal, (REAL)(0.0))
        samp.ln_prob = ln_p_all(samp)
        return samp
    end

    function exec_map(err_scores::String, mat_scores::String, pat_scores::String, ini_file::String;
                      seed::INT = (I)(0), iter::INT = (I)(100000), thin::INT = (I)(10), burnin::INT = (I)(10),
                      REAL::Type{<:Real} = Float32) where {INT <: Integer}
        samp = SamplerType.init(err_scores, mat_scores, pat_scores, ini_file, INT, REAL)
        map, ln_probs = SamplerType.sample_map!(samp, seed = seed, iter = iter, thin = thin, burnin = burnin)
        return (map, ln_probs)
    end


    function _ping_sampler(err_scores::String, mat_scores::String, pat_scores::String, ini_file::String)
        samp = SamplerType.init(err_scores, mat_scores, pat_scores, ini_file)
        sampled = SamplerType.sample_all!(samp)
        return sampled
    end

    function __ping_sampler()
        samp = SamplerType.init("../../simulationTree/err.score.txt", "../../simulationTree/mat.score.txt",
                            "../../simulationTree/pat.score.txt", "./simpleModel.ini")
        sampled = SamplerType.sample_all!(samp)
        return sampled
    end

    #
    # private funcitons
    #
    function _parseData(err_score_path::String,
                        pat_score_path::String,
                        mat_score_path::String,
                        INT::Type{<:Integer} = Int32, REAL::Type{<:Real}=Float32)::Array{REAL, 3}
        err_score::Array{REAL, 2} = InputParser.parse_input_summary(err_score_path)
        mat_score::Array{REAL, 2} = InputParser.parse_input_summary(mat_score_path)
        pat_score::Array{REAL, 2} = InputParser.parse_input_summary(pat_score_path)
        # param::Parameters{REAL}  = InputParser.parse_config_file(param_path)
        S::INT, M::INT = size(err_score)
        ln_p::Array{REAL, 3} = convert.(REAL, zeros(S,M,3))
        for (s,m) in Iterators.product((INT)(1):S,(INT)(1):M)
            ln_p[s,m,1] = err_score[s,m]
            ln_p[s,m,2] = mat_score[s,m]
            ln_p[s,m,3] = pat_score[s,m]
        end
        return ln_p
    end

    function _exp_normalize!(ln_p::Array{R,1})::Nothing where {R <: Real}
        ln_p .= (ln_p .- maximum(ln_p))
        ln_p .= exp.(ln_p)
        ln_p .= ln_p ./ sum(ln_p)
        return nothing
    end

    function _exp_normalize!(ln_p::AbstractArray{R,1})::Nothing where {R <: Real}
        ln_p .= (ln_p .- maximum(ln_p))
        ln_p .= exp.(ln_p)
        ln_p .= ln_p ./ sum(ln_p)
        return nothing
    end

    function _log_sum_exp(ln_p::Array{R,1})::R where {R <: Real}
        maxVal::R = maximum(ln_p)
        return log(ℯ, sum(exp.( (ln_p .- maxVal) ))) + maxVal
    end

    #= DONE:
        add: nested ln_P_Link(L1, L2) as another function.
    =#
    function _ln_p_Link(L::TableGraph{I, R})::R where {I <: Integer, R <: Real}
        ans::R = (R)(0.0)
        ln_p = zeros(R, L.V)
        for from in (I)(1):L.V
            ln_p .= (R)(0.0)
            ln_p .= log.(ℯ, L.W[:, from])
            _exp_normalize!(ln_p)
            ln_p .= log.(ℯ, ln_p)
            for to in (I)(1):L.V
                ans += L.E[to, from] * ln_p[to]
            end
        end
        return ans
    end

    function _ln_p_Link(L::TableGraph{I, R}, usage_s::Dict{I,Array{I,1}}, s_s::Array{I, 1} )::R where {I <: Integer, R <: Real}
        ans::R = (R)(0.0)
        ln_p = zeros(R, L.V)
        for from in (I)(1):L.V
            ln_p .= (R)(0.0)
            ln_p .= log.(ℯ, L.W[:, from])
            _exp_normalize!(view(ln_p, usage_s[s_s[from]]))
            view(ln_p, usage_s[s_s[from]]) .= log.(ℯ, ln_p[usage_s[s_s[from]]])
            for to in usage_s[s_s[from]]
                ans += L.E[to, from] * ln_p[to]
            end
        end
        return ans
    end

    function _ln_p_crp(usage::Dict{I,Array{I,1}}, α::R)::R where {R <: Real, I <: Integer}
        ans::R = (R)(0.0)
        (0 ∉ keys(usage)) && (ans += (length(keys(usage))) * log(ℯ, α))
        (0 ∈ keys(usage)) && (ans += (length(keys(usage))-(I)(1)) * log(ℯ, α))

        totalNum::I = 0
        for (cluster, indexes) in usage
            (cluster != 0) && (
                number::I = (I)(length(indexes));
                ans += lgamma(number);
                totalNum += number;
            )
        end
        ans += lgamma(α) - lgamma(α+totalNum) #lnAF(α, totalNum)
        return ans
    end

    # x ∈ {0,1}, x = 1 w.p. p
    function _ln_p_ber(x::I, p::R)::R where {R <: Real, I <: Integer}
        return x * log(ℯ, p) + (1-x) * log(ℯ, (1.0-p))
    end

    function _ln_p_beta(p::R, α::R, β::R)::R where {R <: Real}
        return (α-(R)(1.0))*log(ℯ, p) + (β-(R)(1.0))*log(ℯ,((R)(1.0)-p)) - SpecialFunctions.lbeta(α, β)
    end

    function _ln_p_er(er::Array{I,1}, erset::Set{I}, p_err::R)::R where {R <: Real, I <: Integer}
        ln_p_err::R = log(ℯ, p_err)
        ln_1m_p_err::R = log(ℯ, (R)(1.0) - p_err)
        return (R)(length(erset)) * ln_p_err + (R)(length(er) - length(erset)) * ln_1m_p_err
    end

    function _ln_p_b(B::Dict{Tuple{I,I},I},
                     usage_s::Dict{I, Array{I,1}},
                     usage_v::Dict{I, Array{I,1}},
                     δ::Array{R, 1})::R where {I <: Integer, R <: Real}
        ans::R = 0.0
        for c in keys(usage_s)
            for m in keys(usage_v)
                (c != 0 && m != 0) &&
                (value::I = B[c,m]; ans += log(ℯ, δ[value]);)
            end
        end
        return ans
    end

    function _is_error_valid(B::Dict{Tuple{I,I},I},
                             usage_s::Dict{I, Array{I,1}},
                             usage_v::Dict{I, Array{I,1}},
                             erset::Set{I},
                             tree_cache::BuffPhyloMatrix{I};
                             blocktype::I = (I)(1), debug = false)::Bool where {I <: Integer}
        BuffPhyloMatrixType.update!(tree_cache, B, usage_s, usage_v, blocktype = (I)(1))
        # (debug) && (print("erset: "); println(erset))
        # (debug) && (print("istree: "); println(BuffPhyloMatrixType.istree(tree_cache)))
        return (length(erset) == 0) || (!BuffPhyloMatrixType.istree(tree_cache))
    end

    function _is_main_valid(usage_s::Dict{I, Array{I,1}},
                            usage_v::Dict{I, Array{I,1}},
                            B::Dict{Tuple{I,I},I})::Bool where {I <: Integer}
        # return true
        ans::Bool = true
        for m in keys(usage_v)
            num_nonunique::I = 0
            num_unique::I = 0
            for c in keys(usage_s)
                num_nonunique += (I)(B[c,m] == 1) + (I)(B[c,m] == 2)
                num_unique += (I)(B[c,m] == 3)
            end
            ans = ans && (!(num_nonunique > 0 && num_unique > 0))
        end
        for m in keys(usage_v)
            num_shared::I = 0
            num_merged::I = 0
            for c in keys(usage_s)
                num_shared += (I)(B[c,m] == 1)
                num_merged += (I)(B[c,m] == 2)
            end
            ans = ans && ( (num_merged == 0) || (num_shared > 0 &&  num_merged > 0) )
        end
        return ans
    end

    function _ln_p_v(Z::Array{I, 2},
                     usage_s::Dict{I, Array{I,1}},
                     usage_v::Dict{I, Array{I,1}},
                     B::Dict{Tuple{I,I},I},
                     erset::Set{I},
                     tree_cache::BuffPhyloMatrix{I},
                     ln_p_penalties::Array{R, 1})::R where {I <: Integer, R <: Real}
        isvalid::Bool = _is_main_valid(usage_s, usage_v, B) && _is_error_valid(B, usage_s, usage_v, erset, tree_cache, blocktype = (I)(1))
        return (R)(isvalid) * ln_p_penalties[1] + (R)(!isvalid) * ln_p_penalties[2]
    end

    function _ln_p_a(a::Array{I, 1}, p_merge::R)::R where {I <: Integer, R <: Real}
        ans::R = (R)(0.0)
        ln_p_merge::R    = log(ℯ, p_merge)
        ln_1m_p_merge::R = log(ℯ, 1.0 - p_merge)
        for x in a
            ans += (x == 2) * ln_p_merge + (x == 1) * ln_1m_p_merge
        end
        return ans
    end

    function _ln_p_f(f::Array{R, 1}, λ::Array{R, 2}, er::Array{I, 1})::R where {I <: Integer, R <: Real}
        ans::R = (R)(0.0)
        for j in 1:length(f)
            ans += _ln_p_beta(f[j], λ[er[j],1], λ[er[j],2])
        end
        return ans
    end

    function _ln_p_g(g::Array{I, 1}, β::Array{R, 1}, er::Array{I, 1})::R where {I <: Integer, R <: Real}
        ans::R = (R)(0.0)
        for j in 1:length(g)
            ans += _ln_p_ber(g[j]-(I)(1), β[er[j]])
        end
        return ans
    end

    # DONE: add function argumant of S^(1)_s, usage_s^(1), S^(2)_s, usage_s^(2)
    # samp.u[j] == i => samp.u[j] ∈ usage_n[s_n[i]
    # samp.u[j] != i => samp.u[j] ∉ usage_n[s_n[i]
    function _ln_p_y(Z::Array{I, 2},
                     H::Array{I, 2},
                     s_s::Array{I, 1},
                     s_n::Array{I, 1},
                     s_v::Array{I, 1},
                     usage_n::Dict{I, Array{I,1}},
                     B::Dict{Tuple{I,I},I},
                     a::Array{I, 1},
                     f::Array{R, 1},
                     g::Array{I, 1},
                     u::Array{I, 1},
                     er::Array{I, 1},
                     param::Parameters{I, R};
                     range_s::AbstractArray{I, 1} = (I)(1):(I)(0),
                     range_v::AbstractArray{I, 1} = (I)(1):(I)(0))::R where {I <: Integer, R <: Real}
        ans::R = (R)(0.0)
        S::I, M::I = size(Z)
        if length(range_s) == 0
            range_s = (I)(1):S
        end
        if length(range_v) == 0
            range_v = (I)(1):M
        end
        for j in range_v
            for i in range_s
                c::I, m::I = (s_s[i], s_v[j])
                if     B[c,m] == 1 # shared
                    ans += (R)(g[j] == 2) * _ln_p_ber(H[i,j]-(I)(1), param.p_hap)
                    ans += (R)(g[j] == 1) * _ln_p_ber(H[i,j]-(I)(1), (R)(1.0) - param.p_hap)
                    ans += _ln_p_ber(Z[i,j]-(I)(1), f[j])
                elseif B[c,m] == 2 # merged
                    ans += (R)(g[j] == 2) * _ln_p_ber(H[i,j]-(I)(1), param.p_hap)
                    ans += (R)(g[j] == 1) * _ln_p_ber(H[i,j]-(I)(1), (R)(1.0) - param.p_hap)
                    ans += (R)(a[i] == 1) * _ln_p_ber(Z[i,j]-(I)(1), param.p_back)
                    ans += (R)(a[i] == 2) * _ln_p_ber(Z[i,j]-(I)(1), f[j])
                elseif B[c,m] == 3 # unique
                    ans += (R)(g[j] == 2) * _ln_p_ber(H[i,j]-(I)(1), param.p_hap)
                    ans += (R)(g[j] == 1) * _ln_p_ber(H[i,j]-(I)(1), (R)(1.0) - param.p_hap)
                    ans += (R)(i ∉ usage_n[s_n[u[j]]]) * _ln_p_ber(Z[i,j]-(I)(1), param.p_unique)
                    ans += (R)(i ∈ usage_n[s_n[u[j]]]) * _ln_p_ber(Z[i,j]-(I)(1), f[j])
                elseif B[c,m] == 4 # error
                    ans += _ln_p_ber(H[i,j]-(I)(1), (R)(0.50))
                    ans += _ln_p_ber(Z[i,j]-(I)(1), f[j])
                else
                    throw("unexpected Block type @ ln_P_Y")
                end
            end
        end
        return ans
    end

    function _ln_p_data_y(ln_p_data::Array{R, 3},
                          Z::Array{I, 2},
                          H::Array{I, 2})::R where {I <: Integer, R <: Real}
        S::I, M::I = size(Z)
        ans::R = (R)(0.0)
        for (i,j) in Iterators.product(1:S, 1:M)
            t::I = (I)(Z[i,j] == 2) * ( (I)(1) + H[i,j] ) + (I)(Z[i,j] == 1) * 1
            ans += ln_p_data[i,j,t]
        end
        return ans
    end

    function _setnums_candidates(usage::Dict{I, Array{I, 1}}, unused::Set{I},
                                 now_cluster::I, α::R;
                                 rm_irrerevant::Bool = true)::Tuple{Array{R, 1}, Array{I, 1}} where {I <: Integer, R <: Real}
        candidates::Array{I, 1} = collect(keys(usage))
        nums::Array{R, 1} = [ (R)(length(usage[k])) for k in keys(usage)]

        if rm_irrerevant && (0 ∈ keys(usage))
            # irrerevant_at::I = findfirst(candidates, (I)(0))
            irrerevant_at::I = something(findfirst(isequal((I)(0)),candidates), 0)
            deleteat!(candidates, irrerevant_at)
            deleteat!(nums, irrerevant_at)
        end

        if now_cluster != 0
            # now_cluster_at::I = findfirst(candidates, now_cluster)
            now_cluster_at::I = something(findfirst(isequal((I)(now_cluster)),candidates), 0)
            if nums[now_cluster_at] == (R)(1.0)
                deleteat!(candidates, now_cluster_at)
                deleteat!(nums, now_cluster_at)
            end
        end

        new_cluster::I = 0
        (length(unused) >  0) && (new_cluster = first(unused))
        (length(unused) == 0) && (new_cluster = now_cluster)
        push!(candidates, new_cluster)
        push!(nums, α)

        return (nums, candidates)
    end

    function _sample_crp(now_cluster::I, usage::Dict{I, Array{I, 1}}, unused::Set{I},
                         alpha::R)::I where {I <: Integer, R <: Real}
        @assert now_cluster ∈ keys(usage)
        nums::Array{R, 1}, candidates::Array{I, 1} = _setnums_candidates(usage, unused, now_cluster, alpha)
        ln_p::Array{R, 1} = convert.(R, log.(ℯ, nums))
        _exp_normalize!(ln_p)

        sampled_cluster::I = candidates[ argmax(RandomUtil.sample_multinomial((I)(1), ln_p)) ]
        return sampled_cluster
    end

    function _update_cluster!(prev_cluster::I, next_cluster::I, array::Array{I, 1}, usage::Dict{I, Array{I, 1}}, unused::Set{I},
                              data_index::I)::Nothing where {I <: Integer}
        @assert 0 ∉ unused
        array[data_index] = next_cluster
        # prev_data_index::I = findfirst(usage[prev_cluster], data_index)
        prev_data_index::I = something(findfirst(isequal((I)(data_index)),usage[prev_cluster]), 0)
        deleteat!(usage[prev_cluster], prev_data_index)
        if length(usage[prev_cluster]) == 0
            pop!(usage, prev_cluster)
            (prev_cluster != 0) && push!(unused, prev_cluster)
        end

        if next_cluster ∈ keys(usage)
            push!(usage[next_cluster], data_index)
        else
            usage[next_cluster] = [data_index]
            (next_cluster != 0) && pop!(unused, next_cluster)
        end

        return nothing
    end

    function _update_subset!(element::I, prev_state::I, next_state::I, set::Set{I}, array::Array{I, 1};
                             positive::I = (I)(2), negative::I = (I)(1))::Nothing where {I <: Integer}
        (prev_state == next_state) && (return nothing)

        if next_state == positive
            push!(set, element)
            array[element] = positive
        elseif next_state == negative
            pop!(set, element)
            array[element] = negative
        end
        return nothing
    end

    function _update_cache_in_error!(element::I, prev_state::I, next_state::I,
                                     v::AbstractArray{I, 1},
                                     tree_cache::BuffPhyloMatrix{I};
                                     positive::I = (I)(2), negative::I = (I)(1))::Nothing where {I <: Integer}
        (prev_state == next_state) && (return nothing)

        if next_state == positive
            BuffPhyloMatrixType.add!(tree_cache, element, v)
        elseif next_state == negative
            BuffPhyloMatrixType.rm!(tree_cache, element)
        end
        return nothing
    end

    function _update_penalty!(iter::I,
                              ln_p_penalty::Array{R, 1},
                              ln_p_ladders::Array{R, 2},
                              period::I)::Nothing where {I <: Integer, R <: Real}
        ln_p_penalty .= ln_p_ladders[(div(iter, period) % length(ln_p_ladders[:,1]))+1, :]
        return nothing
    end

    function _sample_link(D::Array{R, 2}, i::I) where{I <: Integer, R <: Real}
        @assert size(D)[1] == size(D)[2]
        ln_p::Array{R, 1} = log.(ℯ, D[:, i])
        _exp_normalize!(ln_p)
        return argmax( RandomUtil.sample_multinomial((I)(1), ln_p) )
    end

    # DONE: modify: D::AbstractArray{R, 2}
    function _sample_link(D::AbstractArray{R, 2}, i::I) where{I <: Integer, R <: Real}
        @assert size(D)[1] == size(D)[2]
        ln_p::Array{R, 1} = log.(ℯ, D[:, i])
        _exp_normalize!(ln_p)
        return argmax( RandomUtil.sample_multinomial((I)(1), ln_p) )
    end

    function _rm_cluster!(S::Set{I}, usage::Dict{I, Array{I, 1}}, unused::Set{I}, s_s::Array{I, 1})::I where {I <: Integer}
        if length(S) > 0
            i::I = first(S)
            k::I = s_s[i]
            for n in usage[k]
                s_s[n] = -1
            end
            pop!(usage, k)
            push!(unused, k)
            return k
        else
            return -1
        end
    end

    function _add_cluster!(S::Set{I}, usage::Dict{I, Array{I, 1}}, unused::Set{I}, s_s::Array{I, 1})::I where {I <: Integer}
        if length(S) > 0
            i::I = first(S)
            k::I = pop!(unused)
            usage[k] = collect(S)
            for n in S
                s_s[n] = k
            end
            return k
        else
            return -1
        end
    end

    # DONE: add: _propose_nested_link(smap.L.W, i, usage_s, s_s)::tableGraph{I, R}
    function _propose_nested_link!(L::TableGraph{I, R},
                                   usage_s::Dict{I, Array{I, 1}},
                                   s_s::Array{I, 1}, focus::Set{I})::Nothing where {I <: Integer, R <: Real}
        for v in focus
            up_cluster::I = s_s[v]
            ln_p::Array{R, 1} = log.(ℯ, L.W[usage_s[up_cluster], v])
            _exp_normalize!(ln_p)
            propose::I = usage_s[up_cluster][argmax( RandomUtil.sample_multinomial((I)(1), ln_p) )]
            for to in (I)(1):L.V
                (to == propose) && (TableGraphType.addedge!(L, v, to))
                (to != propose) && (TableGraphType.rmedge!(L, v, to))
            end
        end
    end

    # DONE: add: _set_cluster_by_link!(L_n, s_n, usage_n, unused_n)
    function _set_cluster_by_link!(L::TableGraph{I, R},
                                   s_n::Array{I, 1},
                                   usage_n::Dict{I, Array{I, 1}},
                                   unused_n::Set{I})::Nothing where {I <: Integer, R <: Real}
        N::I = size(s_n)[1]
        unused_n = deepcopy(Set{I}())
        s_n    .= (I)(0)
        for c in keys(usage_n)
            delete!(usage_n, c)
        end
        c_now::I = (I)(1)
        for n in (I)(1):N
            if s_n[n] == (I)(0)
                tempSet::Set{I} = Set{I}()
                TableGraphType.dfs!(L, n, tempSet)
                for nn in tempSet
                    s_n[nn] = c_now
                end
                c_now += (I)(1)
            end
        end
        for n in (I)(1):N
            if s_n[n] ∈ keys(usage_n)
                push!(usage_n[s_n[n]], n)
            else
                usage_n[s_n[n]] = [n]
            end
        end
        for c in c_now:N
            push!(unused_n, c)
        end
    end

    function _update_s_s_by_group_diff!(s_pp::Set{I}, s_pn::Set{I}, s_np::Set{I}, s_nn::Set{I},
                                        usage_s::Dict{I, Array{I, 1}},
                                        unused_s::Set{I}, s_s::Array{I, 1})::Tuple{Array{I, 1}, Array{I, 1}} where {I <: Integer}
        d1::I = _rm_cluster!(s_pp, usage_s, unused_s, s_s)
        d2::I = _rm_cluster!(s_pn, usage_s, unused_s, s_s)
        c1::I = _add_cluster!(s_np, usage_s, unused_s, s_s)
        c2::I = _add_cluster!(s_nn, usage_s, unused_s, s_s)
        total::I = (I)(0)
        for k in keys(usage_s)
            total += length(usage_s[k])
        end
        # @assert minimum(s_s) >= 0
        # @assert total == length(s_s)
        return (filter(x -> x > -1, [d1, d2]), filter(x -> x > -1, [c1, c2]) )
    end
end

# using sampler
# @time SamplerType.__ping_sampler()
