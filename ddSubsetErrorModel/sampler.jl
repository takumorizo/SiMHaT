Include("config.jl")
Include("inputFileParser.jl")
Include("random.jl")
Include("buffPhyloMatrix.jl")
Include("tableGraph.jl")
Include("distanceParser.jl")

"""
ddSubsetErrorModel sampling script
"""
module __sampler
    using ..config
    using ..inputParser
    using Distributions
    using ..buffPhyloMatrix
    using ..random
    using ..tableGraph
    using SpecialFunctions
    function parseData(errScorePath::String,
                       patScorePath::String,
                       matScorePath::String)::Array{REAL, 3}
        errScore::Array{REAL, 2} = inputParser.parseInputSummary(errScorePath)
        matScore::Array{REAL, 2} = inputParser.parseInputSummary(matScorePath)
        patScore::Array{REAL, 2} = inputParser.parseInputSummary(patScorePath)
        # param::Parameters{REAL}  = inputParser.parseConfigFile(paramPath)
        S::INT, M::INT = size(errScore)
        lnP::Array{REAL, 3} = convert.(REAL, zeros(S,M,3))
        for (s,m) in Iterators.product((INT)(1):S,(INT)(1):M)
            lnP[s,m,1] = errScore[s,m]
            lnP[s,m,2] = matScore[s,m]
            lnP[s,m,3] = patScore[s,m]
        end
        return lnP
    end

    function exp_normalize!(ln_p::Array{R,1})::Nothing where {R <: Real}
        ln_p .= (ln_p .- maximum(ln_p))
        ln_p .= exp.(ln_p)
        ln_p .= ln_p ./ sum(ln_p)
        return nothing
    end

    function log_sum_exp(ln_p::Array{R,1})::R where {R <: Real}
        maxVal::R = maximum(ln_p)
        return log(ℯ, sum(exp.( (ln_p .- maxVal) ))) + maxVal
    end

    function ln_P_Link(L::TableGraph{I, R})::R where {I <: Integer, R <: Real}
        ans::R = (R)(0.0)
        ln_p = zeros(R, L.V)
        for from in (I)(1):L.V
            ln_p .= (R)(0.0)
            ln_p .= log.(ℯ, L.W[:, from])
            exp_normalize!(ln_p)
            ln_p .= log.(ℯ, ln_p)
            for to in (I)(1):L.V
                ans += L.E[to, from] * ln_p[to]
            end
        end
        return ans
    end

    function ln_P_CRP(usage::Dict{I,Array{I,1}}, α::R)::R where {R <: Real, I <: Integer}
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
    function ln_P_ber(x::I, p::R)::R where {R <: Real, I <: Integer}
        return x * log(ℯ, p) + (1-x) * log(ℯ, (1.0-p))
    end

    function ln_P_beta(p::R, α::R, β::R)::R where {R <: Real}
        return (α-(R)(1.0))*log(ℯ, p) + (β-(R)(1.0))*log(ℯ,((R)(1.0)-p)) - SpecialFunctions.lbeta(α, β)
    end

    function ln_P_er(er::Array{I,1}, erSet::Set{I}, p_err::R)::R where {R <: Real, I <: Integer}
        ln_p_err::R = log(ℯ, p_err)
        ln_1m_p_err::R = log(ℯ, (R)(1.0) - p_err)
        return (R)(length(erSet)) * ln_p_err + (R)(length(er) - length(erSet)) * ln_1m_p_err
    end

    function ln_P_B(B::Dict{Tuple{I,I},I},
                    usageS::Dict{I, Array{I,1}},
                    usageV::Dict{I, Array{I,1}},
                    δ::Array{R, 1})::R where {I <: Integer, R <: Real}
        ans::R = 0.0
        for c in keys(usageS)
            for m in keys(usageV)
                (c != 0 && m != 0) &&
                (value::I = B[c,m]; ans += log(ℯ, δ[value]);)
            end
        end
        return ans
    end

    function isErrorValid(B::Dict{Tuple{I,I},I},
                          usageS::Dict{I, Array{I,1}},
                          usageV::Dict{I, Array{I,1}},
                          erSet::Set{I},
                          treeCache::BuffPhyloMatrix{I};
                          blockType::I = (I)(1), debug = false)::Bool where {I <: Integer}
        buffPhyloMatrix.update!(treeCache, B, usageS, usageV, blockType = (I)(1))
        # (debug) && (print("erSet: "); println(erSet))
        # (debug) && (print("isTree: "); println(buffPhyloMatrix.isTree(treeCache)))
        return (length(erSet) == 0) || (!buffPhyloMatrix.isTree(treeCache))
    end

    function isMainValid(usageS::Dict{I, Array{I,1}},
                         usageV::Dict{I, Array{I,1}},
                         B::Dict{Tuple{I,I},I})::Bool where {I <: Integer}
        # return true
        ans::Bool = true
        for m in keys(usageV)
            numNonUnique::I = 0
            numUnique::I = 0
            for c in keys(usageS)
                numNonUnique += (I)(B[c,m] == 1) + (I)(B[c,m] == 2)
                numUnique += (I)(B[c,m] == 3)
            end
            ans = ans && (!(numNonUnique > 0 && numUnique > 0))
        end
        for m in keys(usageV)
            numShared::I = 0
            numMerged::I = 0
            for c in keys(usageS)
                numShared += (I)(B[c,m] == 1)
                numMerged += (I)(B[c,m] == 2)
            end
            ans = ans && ( (numMerged == 0) || (numShared > 0 &&  numMerged > 0) )
        end
        return ans
    end

    function ln_P_V(Z::Array{I, 2},
                    usageS::Dict{I, Array{I,1}},
                    usageV::Dict{I, Array{I,1}},
                    B::Dict{Tuple{I,I},I},
                    erSet::Set{I},
                    treeCache::BuffPhyloMatrix{I},
                    ln_p_penalties::Array{R, 1})::R where {I <: Integer, R <: Real}
        isValid::Bool = isMainValid(usageS, usageV, B) && isErrorValid(B, usageS, usageV, erSet, treeCache, blockType = (I)(1))
        return (R)(isValid) * ln_p_penalties[1] + (R)(!isValid) * ln_p_penalties[2]
    end

    function ln_P_a(a::Array{I, 1}, p_merge::R)::R where {I <: Integer, R <: Real}
        ans::R = (R)(0.0)
        ln_p_merge::R    = log(ℯ, p_merge)
        ln_1m_p_merge::R = log(ℯ, 1.0 - p_merge)
        for x in a
            ans += (x == 2) * ln_p_merge + (x == 1) * ln_1m_p_merge
        end
        return ans
    end

    function ln_P_f(f::Array{R, 1}, λ::Array{R, 2}, er::Array{I, 1})::R where {I <: Integer, R <: Real}
        ans::R = (R)(0.0)
        for j in 1:length(f)
            ans += ln_P_beta(f[j], λ[er[j],1], λ[er[j],2])
        end
        return ans
    end

    function ln_P_g(g::Array{I, 1}, β::Array{R, 1}, er::Array{I, 1})::R where {I <: Integer, R <: Real}
        ans::R = (R)(0.0)
        for j in 1:length(g)
            ans += ln_P_ber(g[j]-(I)(1), β[er[j]])
        end
        return ans
    end

    function ln_P_Y(Z::Array{I, 2},
                    H::Array{I, 2},
                    S_s::Array{I, 1},
                    S_v::Array{I, 1},
                    B::Dict{Tuple{I,I},I},
                    a::Array{I, 1},
                    f::Array{R, 1},
                    g::Array{I, 1},
                    u::Array{I, 1},
                    er::Array{I, 1},
                    param::Parameters{I, R};
                    rangeS::AbstractArray{I, 1} = (I)(1):(I)(0),
                    rangeV::AbstractArray{I, 1} = (I)(1):(I)(0))::R where {I <: Integer, R <: Real}
        ans::R = (R)(0.0)
        S::I, M::I = size(Z)
        if length(rangeS) == 0
            rangeS = (I)(1):S
        end
        if length(rangeV) == 0
            rangeV = (I)(1):M
        end
        for j in rangeV
            for i in rangeS
                c::I, m::I = (S_s[i], S_v[j])
                if     B[c,m] == 1 # shared
                    ans += (R)(g[j] == 2) * ln_P_ber(H[i,j]-(I)(1), param.p_hap)
                    ans += (R)(g[j] == 1) * ln_P_ber(H[i,j]-(I)(1), (R)(1.0) - param.p_hap)
                    ans += ln_P_ber(Z[i,j]-(I)(1), f[j])
                elseif B[c,m] == 2 # merged
                    ans += (R)(g[j] == 2) * ln_P_ber(H[i,j]-(I)(1), param.p_hap)
                    ans += (R)(g[j] == 1) * ln_P_ber(H[i,j]-(I)(1), (R)(1.0) - param.p_hap)
                    ans += (R)(a[i] == 1) * ln_P_ber(Z[i,j]-(I)(1), param.p_back)
                    ans += (R)(a[i] == 2) * ln_P_ber(Z[i,j]-(I)(1), f[j])
                elseif B[c,m] == 3 # unique
                    ans += (R)(g[j] == 2) * ln_P_ber(H[i,j]-(I)(1), param.p_hap)
                    ans += (R)(g[j] == 1) * ln_P_ber(H[i,j]-(I)(1), (R)(1.0) - param.p_hap)
                    ans += (R)(u[j] != i) * ln_P_ber(Z[i,j]-(I)(1), param.p_unique)
                    ans += (R)(u[j] == i) * ln_P_ber(Z[i,j]-(I)(1), f[j])
                elseif B[c,m] == 4 # error
                    ans += ln_P_ber(H[i,j]-(I)(1), (R)(0.50))
                    ans += ln_P_ber(Z[i,j]-(I)(1), f[j])
                else
                    throw("unexpected Block type @ ln_P_Y")
                end
            end
        end
        return ans
    end

    function ln_P_Data_Y(lnPData::Array{R, 3},
                         Z::Array{I, 2},
                         H::Array{I, 2})::R where {I <: Integer, R <: Real}
        S::I, M::I = size(Z)
        ans::R = (R)(0.0)
        for (i,j) in Iterators.product(1:S, 1:M)
            t::I = (I)(Z[i,j] == 2) * ( (I)(1) + H[i,j] ) + (I)(Z[i,j] == 1) * 1
            ans += lnPData[i,j,t]
        end
        return ans
    end

    function setNumsCandidates(usage::Dict{I, Array{I, 1}}, unUsed::Set{I},
                               nowCluster::I, α::R;
                               rmIrrerevant::Bool = true)::Tuple{Array{R, 1}, Array{I, 1}} where {I <: Integer, R <: Real}
        candidates::Array{I, 1} = collect(keys(usage))
        nums::Array{R, 1} = [ (R)(length(usage[k])) for k in keys(usage)]

        if rmIrrerevant && (0 ∈ keys(usage))
            # irrerevantAt::I = findfirst(candidates, (I)(0))
            irrerevantAt::I = something(findfirst(isequal((I)(0)),candidates), 0)
            deleteat!(candidates, irrerevantAt)
            deleteat!(nums, irrerevantAt)
        end

        if nowCluster != 0
            # nowClusterAt::I = findfirst(candidates, nowCluster)
            nowClusterAt::I = something(findfirst(isequal((I)(nowCluster)),candidates), 0)
            if nums[nowClusterAt] == (R)(1.0)
                deleteat!(candidates, nowClusterAt)
                deleteat!(nums, nowClusterAt)
            end
        end

        newCluster::I = 0
        (length(unUsed) >  0) && (newCluster = first(unUsed))
        (length(unUsed) == 0) && (newCluster = nowCluster)
        push!(candidates, newCluster)
        push!(nums, α)

        return (nums, candidates)
    end

    function sampleCRP(nowCluster::I, usage::Dict{I, Array{I, 1}}, unUsed::Set{I},
                       alpha::R)::I where {I <: Integer, R <: Real}
        @assert nowCluster ∈ keys(usage)
        nums::Array{R, 1}, candidates::Array{I, 1} = setNumsCandidates(usage, unUsed, nowCluster, alpha)
        ln_p::Array{R, 1} = convert.(R, log.(ℯ, nums))
        exp_normalize!(ln_p)

        sampledCluster::I = candidates[ argmax(random.sampleMultiNomial((I)(1), ln_p)) ]
        return sampledCluster
    end

    function updateCluster!(prevCluster::I, nextCluster::I, array::Array{I, 1}, usage::Dict{I, Array{I, 1}}, unUsed::Set{I},
                            dataIndex::I)::Nothing where {I <: Integer}
        @assert 0 ∉ unUsed
        array[dataIndex] = nextCluster
        # prevDataIndex::I = findfirst(usage[prevCluster], dataIndex)
        prevDataIndex::I = something(findfirst(isequal((I)(dataIndex)),usage[prevCluster]), 0)
        deleteat!(usage[prevCluster], prevDataIndex)
        if length(usage[prevCluster]) == 0
            pop!(usage, prevCluster)
            (prevCluster != 0) && push!(unUsed, prevCluster)
        end

        if nextCluster ∈ keys(usage)
            push!(usage[nextCluster], dataIndex)
        else
            usage[nextCluster] = [dataIndex]
            (nextCluster != 0) && pop!(unUsed, nextCluster)
        end

        return nothing
    end

    function updateSubset!(element::I, prevState::I, nextState::I, set::Set{I}, array::Array{I, 1};
                           positive::I = (I)(2), negative::I = (I)(1))::Nothing where {I <: Integer}
        (prevState == nextState) && (return nothing)

        if nextState == positive
            push!(set, element)
            array[element] = positive
        elseif nextState == negative
            pop!(set, element)
            array[element] = negative
        end
        return nothing
    end

    function updateCacheInError!(element::I, prevState::I, nextState::I,
                                 v::AbstractArray{I, 1},
                                 treeCache::BuffPhyloMatrix{I};
                                 positive::I = (I)(2), negative::I = (I)(1))::Nothing where {I <: Integer}
        (prevState == nextState) && (return nothing)

        if nextState == positive
            buffPhyloMatrix.add!(treeCache, element, v)
        elseif nextState == negative
            buffPhyloMatrix.rm!(treeCache, element)
        end
        return nothing
    end

    function updatePenalty!(iter::I,
                            ln_p_penalty::Array{R, 1},
                            ln_p_ladders::Array{R, 2},
                            period::I)::Nothing where {I <: Integer, R <: Real}
        ln_p_penalty .= ln_p_ladders[(div(iter, period) % length(ln_p_ladders[:,1]))+1, :]
        return nothing
    end

    function sampleLink(D::Array{R, 2}, i::I) where{I <: Integer, R <: Real}
        @assert size(D)[1] == size(D)[2]
        ln_p::Array{R, 1} = log.(ℯ, D[:, i])
        __sampler.exp_normalize!(ln_p)
        return argmax( random.sampleMultiNomial((I)(1), ln_p) )
    end

    function rmCluster!(S::Set{I}, usage::Dict{I, Array{I, 1}}, unused::Set{I}, S_s::Array{I, 1})::I where {I <: Integer}
        if length(S) > 0
            i::I = first(S)
            k::I = S_s[i]
            for n in usage[k]
                S_s[n] = -1
            end
            pop!(usage, k)
            push!(unused, k)
            return k
        else
            return -1
        end
    end

    function addCluster!(S::Set{I}, usage::Dict{I, Array{I, 1}}, unused::Set{I}, S_s::Array{I, 1})::I where {I <: Integer}
        if length(S) > 0
            i::I = first(S)
            k::I = pop!(unused)
            usage[k] = collect(S)
            for n in S
                S_s[n] = k
            end
            return k
        else
            return -1
        end
    end

    function updateSsByGroupDiff!(S_pp::Set{I}, S_pn::Set{I}, S_np::Set{I}, S_nn::Set{I},
                                  usageS::Dict{I, Array{I, 1}},
                                  unUsedS::Set{I}, S_s::Array{I, 1})::Tuple{Array{I, 1}, Array{I, 1}} where {I <: Integer}
        d1::I = rmCluster!(S_pp, usageS, unUsedS, S_s)
        d2::I = rmCluster!(S_pn, usageS, unUsedS, S_s)
        c1::I = addCluster!(S_np, usageS, unUsedS, S_s)
        c2::I = addCluster!(S_nn, usageS, unUsedS, S_s)
        total::I = (I)(0)
        for k in keys(usageS)
            total += length(usageS[k])
        end
        # @assert minimum(S_s) >= 0
        # @assert total == length(S_s)
        return (filter(x -> x > -1, [d1, d2]), filter(x -> x > -1, [c1, c2]) )
    end
end

module sampler
    using ..config
    using ..inputParser
    using ..__sampler
    using ..random
    using ..buffPhyloMatrix
    using ..tableGraph
    using ..distanceParser
    using Random

    mutable struct Sampler{I,R}
        Z::Array{I, 2} # Z == 1 err, Z == 2 mutation
        H::Array{I, 2} # H == 1 mat, H == 2 pat

        L::TableGraph{I} # Table Link graph

        S_s::Array{I, 1}
        S_v::Array{I, 1}
        usageS::Dict{I, Array{I, 1}}
        usageV::Dict{I, Array{I, 1}}
        unUsedS::Set{I}
        unUsedV::Set{I}

        a::Array{I, 1}  # a == 1 not merged, a == 2 merged
        f::Array{R, 1}  # freq of mutation for each col
        g::Array{I, 1}  # haplotype for each col
        u::Array{I, 1}  # u[j] : an unique sample having mutation j

        p_err::R        # A frequency in which mutation is false positive.
        er::Array{I, 1} # er[j] == 1 : not error, er[j] == 2 : error
        erSet::Set{I}   # {j | er[j] == 2}

        B::Dict{Tuple{I,I}, I}
        # B[c,m] == 1 : shared
        # B[c,m] == 2 : merged
        # B[c,m] == 3 : unique
        # B[c,m] == 4 : error

        treeCache::BuffPhyloMatrix{I}
        lnPData::Array{R, 3}
        param::Parameters{I, R}
        anneal::Annealer{I, R}
        lnProb::R
    end
    export Sampler

    function isValid(samp::Sampler{I, R}; debug::Bool = false)::Bool where {I <: Integer, R <: Real}
        return  __sampler.isMainValid(samp.usageS, samp.usageV, samp.B) &&
                __sampler.isErrorValid(samp.B, samp.usageS, samp.usageV, samp.erSet, samp.treeCache, blockType = (I)(1))
    end

    function ln_P_all(samp::Sampler{I,R}, debug::Bool = false)::R where {I<:Integer, R<:Real}
        ans::R = (R)(0.0)
        S::I, M::I = size(samp.Z)
        ans += __sampler.ln_P_Link(samp.L)
        ans += __sampler.ln_P_beta(samp.p_err, samp.param.γ_e[1], samp.param.γ_e[2])
        ans += __sampler.ln_P_er(samp.er, samp.erSet, samp.p_err)
        ans += __sampler.ln_P_CRP(samp.usageV, samp.param.α_v)
        ans += __sampler.ln_P_B(samp.B, samp.usageS, samp.usageV, samp.param.δ_s)
        ans += __sampler.ln_P_V(samp.Z, samp.usageS, samp.usageV, samp.B,
                                samp.erSet, samp.treeCache, samp.param.ln_p_v)
        ans += __sampler.ln_P_a(samp.a, samp.param.p_merge)
        ans += __sampler.ln_P_f(samp.f, samp.param.λ_s, samp.er)
        ans += __sampler.ln_P_g(samp.g, samp.param.β_s, samp.er)


        ans += __sampler.ln_P_Y(samp.Z, samp.H,
                                samp.S_s, samp.S_v, samp.B,
                                samp.a, samp.f, samp.g, samp.u, samp.er, samp.param)

        ans += __sampler.ln_P_Data_Y(samp.lnPData, samp.Z, samp.H)
        if debug
            print("__sampler.ln_P_Link(samp.L): ");println(__sampler.ln_P_Link(samp.L))
            print("__sampler.ln_P_er(samp.er, samp.erSet, samp.p_err): ");println(__sampler.ln_P_er(samp.er, samp.erSet, samp.p_err))
            print("__sampler.ln_P_CRP(samp.usageV, samp.param.α_v): ");println(__sampler.ln_P_CRP(samp.usageV, samp.param.α_v))
            print("__sampler.ln_P_B(samp.B, samp.usageS, samp.usageV, samp.param.δ_s): ");println(__sampler.ln_P_B(samp.B, samp.usageS, samp.usageV, samp.param.δ_s))
            print("__sampler.ln_P_V(samp.Z, samp.usageS, samp.usageV, samp.B, samp.erSet, samp.treeCache, samp.param.ln_p_v): ");println(__sampler.ln_P_V(samp.Z, samp.usageS, samp.usageV, samp.B, samp.erSet, samp.treeCache, samp.param.ln_p_v))
            print("__sampler.ln_P_a(samp.a, samp.param.p_merge): ");println(__sampler.ln_P_a(samp.a, samp.param.p_merge))
            print("__sampler.ln_P_f(samp.f, samp.param.λ_s, samp.er): ");println(__sampler.ln_P_f(samp.f, samp.param.λ_s, samp.er))
            print("__sampler.ln_P_g(samp.g, samp.param.β_s, samp.er): ");println(__sampler.ln_P_g(samp.g, samp.param.β_s, samp.er))
            print("__sampler.ln_P_Y(samp.Z, samp.H, samp.S_s, samp.S_v, samp.B, samp.a, samp.f, samp.g, samp.u, samp.er, samp.param): ");println(__sampler.ln_P_Y(samp.Z, samp.H, samp.S_s, samp.S_v, samp.B, samp.a, samp.f, samp.g, samp.u, samp.er, samp.param))
            print("__sampler.ln_P_Data_Y(samp.lnPData, samp.Z, samp.H): ");println(__sampler.ln_P_Data_Y(samp.lnPData, samp.Z, samp.H))
        end
        return ans
    end

    function sampleZ!(samp::Sampler{I, R})::Nothing where {I<:Integer, R <: Real}
        S::I, M::I = size(samp.Z)
        t::I = (I)(0)
        ln_p::Array{R,1} = [0.0, 0.0]
        ln_f::Array{R,1} = [0.0, 0.0]
        now_Z::I = (I)(0)
        for j in (I)(1):M
            for i in (I)(1):S
                now_Z = samp.Z[i,j]
                t = 1 + samp.H[i,j]
                ln_p[1] = samp.lnPData[i,j,1]; ln_p[2] = samp.lnPData[i,j,t];
                ln_f[1] = log(ℯ, 1.0 - samp.f[j]); ln_f[2] = log(ℯ, samp.f[j]);
                if     samp.er[j] == 1 && samp.B[(samp.S_s[i], samp.S_v[j])] == 1
                    ln_p[1] += ln_f[1]
                    ln_p[2] += ln_f[2]
                elseif samp.er[j] == 1 && samp.B[(samp.S_s[i], samp.S_v[j])] == 2
                    ln_p[1] += (samp.a[i] == 2) * ln_f[1] + (samp.a[i] == 1) * (R)(log(ℯ, 1.0-samp.param.p_back))
                    ln_p[2] += (samp.a[i] == 2) * ln_f[2] + (samp.a[i] == 1) * (R)(log(ℯ, samp.param.p_back))
                elseif samp.er[j] == 1 && samp.B[(samp.S_s[i], samp.S_v[j])] == 3
                    ln_p[1] += (samp.u[j] == i) * ln_f[1] + (samp.u[j] != i) * (R)(log(ℯ, 1.0-samp.param.p_unique))
                    ln_p[2] += (samp.u[j] == i) * ln_f[2] + (samp.u[j] != i) * (R)(log(ℯ, samp.param.p_unique))
                elseif samp.er[j] == 2
                    ln_p[1] += ln_f[1]
                    ln_p[2] += ln_f[2]
                    isV::Bool        = isValid(samp)
                    ln_p[now_Z]     += (R)(isV) * samp.param.ln_p_v[1] + (R)(!isV) * samp.param.ln_p_v[2]
                    buffPhyloMatrix.edit!(samp.treeCache, i, j, (I)((3 - now_Z) - 1))
                    isV              = isValid(samp)
                    ln_p[3 - now_Z] += (R)(isV) * samp.param.ln_p_v[1] + (R)(!isV) * samp.param.ln_p_v[2]
                end
                __sampler.exp_normalize!(ln_p)
                samp.Z[i,j] = (I)(argmax( random.sampleMultiNomial((I)(1), ln_p) ))
                if samp.er[j] == 2
                    buffPhyloMatrix.edit!(samp.treeCache, i, j, samp.Z[i,j]-(I)(1))
                end
            end
        end
        return nothing
    end

    function sampleH!(samp::Sampler{I, R})::Nothing where {I<:Integer, R <: Real}
        S::I, M::I = size(samp.Z)
        ln_p::Array{R, 1}   = [0.0, 0.0]
        ln_P_hap::Array{R, 1} = log.(ℯ, [ 1.0 - samp.param.p_hap, samp.param.p_hap ])
        ln_P_rn::Array{R, 1}  = convert.(R, log.(ℯ, [0.50, 0.50])) #convert.(R, log.(ℯ, [0.50, 0.50]))
        for j in (I)(1):M
            for i in (I)(1):S
                ln_p[1] = (samp.Z[i,j] == 2) * samp.lnPData[i,j, 1+1]
                ln_p[2] = (samp.Z[i,j] == 2) * samp.lnPData[i,j, 1+2]
                if     samp.er[j] == 1
                    ln_p[1] += (samp.g[j] == 1) * ln_P_hap[2] + (samp.g[j] == 2) * ln_P_hap[1]
                    ln_p[2] += (samp.g[j] == 2) * ln_P_hap[2] + (samp.g[j] == 1) * ln_P_hap[1]
                elseif samp.er[j] == 2
                    ln_p[1] += (samp.g[j] == 1) * ln_P_rn[2] + (samp.g[j] == 2) * ln_P_rn[1]
                    ln_p[2] += (samp.g[j] == 2) * ln_P_rn[2] + (samp.g[j] == 1) * ln_P_rn[1]
                end
                __sampler.exp_normalize!(ln_p)
                samp.H[i,j] = argmax( random.sampleMultiNomial((I)(1), ln_p) )
            end
        end
        return nothing
    end

    function sampleA!(samp::Sampler{I, R})::Nothing where {I <: Integer, R <: Real}
        ln_p::Array{R, 1}   = convert.(R, [0.0, 0.0])
        ln_p_merge::Array{R, 1} = convert.(R, log.(ℯ, [1.0 - samp.param.p_merge, samp.param.p_merge] ))
        ln_p_back::Array{R, 1}  = convert.(R, log.(ℯ, [1.0 - samp.param.p_back, samp.param.p_back] ))
        S::I, M::I = size(samp.Z)
        c::I = (I)(0)
        for i in (I)(1):S
            ln_p .= ln_p_merge
            c = samp.S_s[i]
            for m in keys(samp.usageV)
                if samp.B[c,m] == 2 # block merge
                    for j in samp.usageV[m]
                        ln_p[2] += (samp.Z[i,j] == 2) * log(ℯ, samp.f[j])
                        ln_p[2] += (samp.Z[i,j] == 1) * log(ℯ, 1.0 - samp.f[j])
                        ln_p[1] += (samp.Z[i,j] == 2) * ln_p_back[2]
                        ln_p[1] += (samp.Z[i,j] == 1) * ln_p_back[1]
                    end
                end
            end
            __sampler.exp_normalize!(ln_p)
            samp.a[i] = argmax( random.sampleMultiNomial((I)(1), ln_p) )
        end
        return nothing
    end

    function sampleF!(samp::Sampler{I, R})::Nothing where {I <: Integer, R <: Real}
        m::I = (I)(0)
        λ::Array{R ,1} = [0.0, 0.0]
        S::I, M::I = size(samp.Z)

        for j in (I)(1):M
            m = samp.S_v[j]
            λ[1] = (samp.er[j] == 1) * samp.param.λ_s[1, 1] + (samp.er[j] == 2) * samp.param.λ_s[2, 1]
            λ[2] = (samp.er[j] == 1) * samp.param.λ_s[1, 2] + (samp.er[j] == 2) * samp.param.λ_s[2, 2]
            for c in keys(samp.usageS)
                for i in samp.usageS[c]
                    isAdd::Bool = (samp.B[c,m] == 1 || samp.B[c,m] == 4) ||
                                  (samp.B[c,m] == 2 && samp.a[i] == 2)   ||
                                  (samp.B[c,m] == 3 && samp.u[j] == i)
                    λ[1] += (R)(isAdd) * (samp.Z[i,j] == 2)
                    λ[2] += (R)(isAdd) * (samp.Z[i,j] == 1)
                end
            end
            samp.f[j] = random.sampleBeta(λ[1], λ[2])
        end
        return nothing
    end

    function sampleG!(samp::Sampler{I, R})::Nothing where {I <: Integer, R <: Real}
        S::I, M::I = size(samp.Z)
        m::I = (I)(0)
        ln_p::Array{R, 1} = convert.(R,log.(ℯ, [0.5, 0.5]))
        ln_P_hap::Array{R, 1} = log.(ℯ, [ 1.0 - samp.param.p_hap, samp.param.p_hap ])
        for j in (I)(1):M
            m = samp.S_v[j]
            ln_p .= samp.param.β_s
            for c in keys(samp.usageS)
                for i in samp.usageS[c]
                    (samp.B[c,m] != 4) && (
                        ln_p[1] += ln_P_hap[ 1 + (samp.H[i,j] == 1) ];
                        ln_p[2] += ln_P_hap[ 1 + (samp.H[i,j] == 2) ];
                    )
                end
            end
            # ln_p_temp::Array{R, 1} = deepcopy(ln_p)
            __sampler.exp_normalize!(ln_p)
            samp.g[j] = argmax( random.sampleMultiNomial((I)(1), ln_p) )
        end
        return nothing
    end

    function sampleU!(samp::Sampler{I, R})::Nothing where {I <: Integer, R <: Real}
        S::I, M::I = size(samp.Z)
        m::I = (I)(0)
        ln_p::Array{R, 1} = zeros(R, S)
        ln_p_unique::Array{R, 1} = convert.(R, log.(ℯ, [1.0 - samp.param.p_unique, samp.param.p_unique]) )
        ln_p_f::Array{R, 1} = [0.0, 0.0]

        for j in (I)(1):M
            m = samp.S_v[j]
            ln_p .= (R)(0.0)
            ln_p_f[1] = log(ℯ, 1.0 - samp.f[j]); ln_p_f[2] = log(ℯ, samp.f[j]);

            ln_p_cons::R = (R)(0.0)
            for c in keys(samp.usageS)
                if samp.B[c,m] == 3
                    for i in samp.usageS[c]
                        consAt::R = (samp.Z[i,j] == 2) * ln_p_unique[2] + (samp.Z[i,j] == 1) * ln_p_unique[1]
                        prosAt::R = (samp.Z[i,j] == 2) * ln_p_f[2]      + (samp.Z[i,j] == 1) * ln_p_f[1]
                        ln_p[i] += (-consAt + prosAt);
                        ln_p_cons += consAt
                    end
                end
            end
            ln_p .+= ln_p_cons

            # ln_p_temp::Array{R, 1} = deepcopy(ln_p)
            __sampler.exp_normalize!(ln_p)
            samp.u[j] = argmax( random.sampleMultiNomial((I)(1), ln_p) )
        end
        return nothing
    end

    function sampleS_s!(samp::Sampler{I, R}, i::I)::Nothing where {I <: Integer, R <: Real}
        c::I = samp.S_s[i]
        prevCluster::I = samp.S_s[i]
        prevBs::Dict{Tuple{I,I}, I} = deepcopy(samp.B)
        nextCluster::I = __sampler.sampleCRP(samp.S_s[i], samp.usageS, samp.unUsedS, samp.param.α_s)
        nextBs::Dict{Tuple{I,I}, I} = deepcopy(samp.B)

        (   # previous cluster disappears
            if length(samp.usageS[prevCluster]) == 1;
                for m in keys(samp.usageV);
                    pop!(nextBs, (prevCluster, m));
                end;
            end;
            # novel cluster appears
            if ( nextCluster ∈ samp.unUsedS ) ||
               ( length(samp.unUsedS) == 0  && nextCluster == prevCluster);
                for m in keys(samp.usageV);
                    (m != 0) && (nextBs[nextCluster, m] = argmax( random.sampleMultiNomial((I)(1), samp.param.δ_s) ));
                    (m == 0) && (nextBs[nextCluster, m] = 4);
                end;
            end;
        )
        # acc from prev S_s, B related
        accPrev::R = (R)(0.0)
        (
            accPrev += __sampler.ln_P_V(samp.Z, samp.usageS, samp.usageV, samp.B, samp.erSet,
                                        samp.treeCache, samp.param.ln_p_v);
            accPrev += __sampler.ln_P_Y(samp.Z, samp.H, samp.S_s, samp.S_v, samp.B, samp.a, samp.f,
                                        samp.g, samp.u, samp.er, samp.param, rangeS = i:i);
        )

        accNext::R = (R)(0.0)
        (
            __sampler.updateCluster!(prevCluster, nextCluster, samp.S_s, samp.usageS, samp.unUsedS, i);
            samp.B = deepcopy(nextBs);
            buffPhyloMatrix.update!(samp.treeCache, samp.B, samp.usageS, samp.usageV);
            accNext += __sampler.ln_P_V(samp.Z, samp.usageS, samp.usageV, samp.B, samp.erSet,
                                        samp.treeCache, samp.param.ln_p_v);
            accNext += __sampler.ln_P_Y(samp.Z, samp.H, samp.S_s, samp.S_v, samp.B, samp.a, samp.f,
                                        samp.g, samp.u, samp.er, samp.param, rangeS = i:i);
        )
        # select prev/next
        accRate::R = min(1.0, exp(accNext - accPrev))
        if rand() > accRate # rejected
            # revert to a previous state
            __sampler.updateCluster!(nextCluster, prevCluster, samp.S_s, samp.usageS, samp.unUsedS, i)
            samp.B = deepcopy(prevBs)
            buffPhyloMatrix.update!(samp.treeCache, samp.B, samp.usageS, samp.usageV)
        end
        return nothing
    end

    function sampleP_err!(samp::Sampler{I, R})::Nothing where {I <: Integer, R <: Real}
        # prev_p_err::R = samp.p_err      # debug
        # prev_p_all::R = ln_P_all(samp)  # debug
        γ_e = [0.0, 0.0]
        γ_e .= samp.param.γ_e
        γ_e[1] += (R)(length(samp.erSet))
        γ_e[2] += (R)(length(samp.er) - length(samp.erSet))
        samp.p_err = random.sampleBeta(γ_e[1], γ_e[2])
        return nothing
    end

    function sampleS_v!(samp::Sampler{I, R}, j::I)::Nothing where {I <: Integer, R <: Real}
        prevCluster::I = samp.S_v[j]
        prevEr::I = samp.er[j]
        prevF::R  = samp.f[j]
        prevBs::Dict{Tuple{I,I}, I} = deepcopy(samp.B)

        nextEr::I = argmax( random.sampleMultiNomial((I)(1), [1.0 - samp.p_err, samp.p_err]) )
        nextCluster::I = (I)(0)
        (nextEr==1) && (nextCluster = __sampler.sampleCRP(samp.S_v[j], samp.usageV, samp.unUsedV, samp.param.α_v) )
        (nextEr==2) && (nextCluster = 0)
        nextF::R  = random.sampleBeta( samp.param.λ_s[nextEr, 1], samp.param.λ_s[nextEr, 2] )
        nextBs::Dict{Tuple{I,I}, I} = deepcopy(samp.B)
        (   # init nextBs
            # previous cluster disappears
            if length(samp.usageV[prevCluster]) == 1;
                for c in keys(samp.usageS);
                    pop!(nextBs, (c, prevCluster));
                end;
            end;
            if nextCluster == 0;
                for c in keys(samp.usageS);
                    nextBs[c, nextCluster] = 4;
                end;
            elseif nextCluster != 0; # novel cluster appears in non error cluster
                if ( nextCluster ∈ samp.unUsedV ) ||
                   ( length(samp.unUsedV) == 0  && nextCluster == prevCluster);
                    for c in keys(samp.usageS);
                        nextBs[c, nextCluster] = argmax( random.sampleMultiNomial((I)(1), samp.param.δ_s) );
                    end;
                end;
            end;
        )
        accPrev::R = (R)(0.0)  # acc from prev S_s, B related
        (
            accPrev += __sampler.ln_P_V(samp.Z, samp.usageS, samp.usageV, samp.B, samp.erSet,
                                        samp.treeCache, samp.param.ln_p_v);
            accPrev += __sampler.ln_P_Y(samp.Z, samp.H, samp.S_s, samp.S_v, samp.B, samp.a, samp.f,
                                        samp.g, samp.u, samp.er, samp.param, rangeV = j:j);
        )
        accNext::R = (R)(0.0) # update prev to next state
        (
            __sampler.updateSubset!(j, prevEr, nextEr, samp.erSet, samp.er);
            __sampler.updateCacheInError!(j, prevEr, nextEr, view(samp.Z, :, j).-(I)(1), samp.treeCache);
            __sampler.updateCluster!(prevCluster, nextCluster, samp.S_v, samp.usageV, samp.unUsedV, j);
            samp.B = deepcopy(nextBs);
            samp.f[j] = nextF;
            buffPhyloMatrix.update!(samp.treeCache, samp.B, samp.usageS, samp.usageV);
            accNext += __sampler.ln_P_V(samp.Z, samp.usageS, samp.usageV, samp.B, samp.erSet,
                                        samp.treeCache, samp.param.ln_p_v);
            accNext += __sampler.ln_P_Y(samp.Z, samp.H, samp.S_s, samp.S_v, samp.B, samp.a, samp.f,
                                        samp.g, samp.u, samp.er, samp.param, rangeV = j:j);
        )
        # select prev/next
        accRate::R = min((R)(1.0), exp(accNext - accPrev))
        if rand() > accRate # rejected and revert to a previous state
            __sampler.updateSubset!(j, nextEr, prevEr, samp.erSet, samp.er)
            __sampler.updateCacheInError!(j, nextEr, prevEr, view(samp.Z, :, j).-(I)(1), samp.treeCache)
            __sampler.updateCluster!(nextCluster, prevCluster, samp.S_v, samp.usageV, samp.unUsedV, j)
            samp.B = deepcopy(prevBs);
            samp.f[j]  = prevF
            buffPhyloMatrix.update!(samp.treeCache, samp.B, samp.usageS, samp.usageV)
        end
        return nothing
    end

    function ln_P_Y_marginal!(samp::Sampler{I, R}, c::I)::R where {I <: Integer, R <: Real}
        ans::R = (R)(0.0)
        ln_p::Array{R, 1} = [0.0, 0.0, 0.0]
        for m in keys(samp.usageV)
            ln_p .= (R)(0.0)
            now_B::I = samp.B[c,m]
            if now_B == 4
                continue
            end
            for b in (I)(1):(I)(3)
                ln_p[b] += log(ℯ, samp.param.δ_s[b])
                samp.B[c,m] = b
                ln_p[b] += __sampler.ln_P_Y(samp.Z, samp.H, samp.S_s, samp.S_v, samp.B, samp.a, samp.f,
                                            samp.g, samp.u, samp.er, samp.param,
                                            rangeS = samp.usageS[c], rangeV = samp.usageV[m])
            end
            samp.B[c,m] = now_B
            ans += __sampler.log_sum_exp(ln_p)
        end
        return ans
    end

    function sampleBPatiallyGibbs!(samp::Sampler{I,R}, c::I, m::I)::Nothing where {I <: Integer, R <: Real }
        (samp.B[c, m] == 4) && (return nothing)
        ln_p::Array{R, 1} = convert.(R, [0.0, 0.0, 0.0])
        for t in (I)(1):(I)(3)
            samp.B[c, m] = t
            ln_p[t] += log(ℯ, samp.param.δ_s[t])
            ln_p[t] += __sampler.ln_P_Y(samp.Z, samp.H, samp.S_s, samp.S_v, samp.B, samp.a, samp.f,
                                        samp.g, samp.u, samp.er, samp.param,
                                        rangeS = samp.usageS[c], rangeV = samp.usageV[m])
        end
        # ln_p_temp::Array{R, 1} = deepcopy(ln_p)
        __sampler.exp_normalize!(ln_p)
        samp.B[c, m] = argmax( random.sampleMultiNomial((I)(1), ln_p) )
        buffPhyloMatrix.update!(samp.treeCache, samp.B, samp.usageS, samp.usageV)
        return nothing
    end

    function sampleL!(samp::Sampler{I, R}, i::I, correct::Bool = true)::Nothing where {I <: Integer, R <: Real}
        prev_to::I = tableGraph.getLink(samp.L, i)
        next_to::I = __sampler.sampleLink(samp.L.W, i)
        S_pp::Set{I} = Set{I}([])
        S_pn::Set{I} = Set{I}([])
        S_np::Set{I} = Set{I}([])
        S_nn::Set{I} = Set{I}([])
        tableGraph.diffGroup!(samp.L, S_pp, S_pn, S_np, S_nn, i, prev_to, next_to, both = true)

        prevUsageS::Dict{I, Array{I, 1}} = deepcopy(samp.usageS)
        prevUnUsedS::Set{I}              = deepcopy(samp.unUsedS)
        prevS_s::Array{I}                = deepcopy(samp.S_s)
        prevBs::Dict{Tuple{I,I}, I}      = deepcopy(samp.B)
        accPrev::R = (R)(0.0)
        (
            if correct;
                accPrev += ln_P_Y_marginal!(samp, samp.S_s[prev_to]);
                if samp.S_s[prev_to] != samp.S_s[next_to];
                    accPrev += ln_P_Y_marginal!(samp, samp.S_s[next_to]);
                end
                accPrev += __sampler.ln_P_V(samp.Z, samp.usageS, samp.usageV, samp.B, samp.erSet,
                                            samp.treeCache, samp.param.ln_p_v);
            end;
        )
        accNext::R = (R)(0.0)
        (# update to proposal state
            tableGraph.rmEdge!(samp.L, i, prev_to);
            tableGraph.addEdge!(samp.L, i, next_to);
            (removed::Array{I, 1}, added::Array{I, 1}) = __sampler.updateSsByGroupDiff!(S_pp, S_pn, S_np, S_nn,
                                                                                        samp.usageS, samp.unUsedS, samp.S_s);
            errM::Set{I} = Set{I}([]);
            for c in removed; for m in keys(samp.usageV);
                (samp.B[c,m] == 4) && (push!(errM, m));
                pop!(samp.B, (c,m));
            end; end;
            for c in added; for m in keys(samp.usageV);
                samp.B[c,m] = 1;
                (m ∈ errM) && (samp.B[c,m] = 4);
            end; end;
            for c in added; for m in keys(samp.usageV);
                sampleBPatiallyGibbs!(samp, c, m);
            end; end;
            if correct;
                accNext += ln_P_Y_marginal!(samp, samp.S_s[prev_to]);
                if samp.S_s[prev_to] != samp.S_s[next_to];
                    accNext += ln_P_Y_marginal!(samp, samp.S_s[next_to]);
                end;
                accNext += __sampler.ln_P_V(samp.Z, samp.usageS, samp.usageV, samp.B, samp.erSet,
                                            samp.treeCache, samp.param.ln_p_v);
            end;
        )
        # select prev/next
        accRate::R = min((R)(1.0), exp(accNext - accPrev))

        if !correct
            return nothing
        elseif rand() > accRate # rejected
            tableGraph.rmEdge!(samp.L, i, next_to)
            tableGraph.addEdge!(samp.L, i, prev_to)
            samp.usageS  = deepcopy(prevUsageS)
            samp.unUsedS = deepcopy(prevUnUsedS)
            samp.S_s     = deepcopy(prevS_s)
            samp.B       = deepcopy(prevBs)
        end
        return nothing
    end

    function sampleB!(samp::Sampler{I, R}, c::I, m::I) where {I <: Integer, R <: Real}
        (samp.B[c, m] == 4) && (return nothing)
        # ln_p_correct::Array{R, 1} = [0.0, 0.0, 0.0]
        ln_p::Array{R, 1} = convert.(R, [0.0, 0.0, 0.0])
        rangeB::AbstractRange{I} = (I)(1):(I)(3)
        (samp.B[c,m] != 1) && (rangeB = (I)(3):(I)(-1):(I)(1))
        for t in (I)(1):(I)(3)
            samp.B[c, m] = t
            buffPhyloMatrix.update!(samp.treeCache, samp.B, samp.usageS, samp.usageV)
            ln_p[t] += log(ℯ, samp.param.δ_s[t])
            ln_p[t] += __sampler.ln_P_Y(samp.Z, samp.H, samp.S_s, samp.S_v, samp.B, samp.a, samp.f,
                                        samp.g, samp.u, samp.er, samp.param,
                                        rangeS = samp.usageS[c], rangeV = samp.usageV[m])
            ln_p[t] += __sampler.ln_P_V(samp.Z, samp.usageS, samp.usageV, samp.B, samp.erSet,
                                        samp.treeCache, samp.param.ln_p_v)
        end
        # ln_p_temp::Array{R, 1} = deepcopy(ln_p)
        __sampler.exp_normalize!(ln_p)
        samp.B[c, m] = argmax( random.sampleMultiNomial((I)(1), ln_p) )
        buffPhyloMatrix.update!(samp.treeCache, samp.B, samp.usageS, samp.usageV)
        # if abs( (ln_p_temp[now_B] - ln_p_temp[new_B]) - (ln_p_correct[now_B] - ln_p_correct[new_B]) ) > 0.001
        #     print("(ln_p_temp[now_B] - ln_p_temp[new_B])"); println((ln_p_temp[now_B] - ln_p_temp[new_B]))
        #     print("ln_p_temp[now_B]"); println(ln_p_temp[now_B])
        #     print("ln_p_temp[new_B]"); println(ln_p_temp[new_B])
        #     print("(ln_p_correct[now_B] - ln_p_correct[new_B])"); println((ln_p_correct[now_B] - ln_p_correct[new_B]))
        #     print("ln_p_correct[now_B]"); println(ln_p_correct[now_B])
        #     print("ln_p_correct[new_B]"); println(ln_p_correct[new_B])
        #     @assert abs( (ln_p_temp[now_B] - ln_p_temp[new_B]) - (ln_p_correct[now_B] - ln_p_correct[new_B]) ) < 0.001
        # end
        return nothing
    end


    function sampleMAP!(samp::Sampler{I, R};
                        seed::I = (I)(0),
                        iter::I = (I)(100000),
                        thin::I = (I)(1),
                        burnin::I = (I)(0),
                        progressCount::I = (I)(100))::Tuple{Sampler{I, R}, Array{R, 1}} where {I <:Integer, R <: Real}
        Random.seed!(seed)
        ln_p_v_true::Array{R, 1} = samp.param.ln_p_v
        mapState::Sampler{I,R} = deepcopy(samp)
        lnProbs::Array{R, 1} = []
        maxLn::R = (R)(-Inf)
        for count in (I)(1):(iter+burnin)
            S::I, M::I = size(samp.Z)
            __sampler.updatePenalty!((I)(count),
                                    samp.param.ln_p_v,
                                    samp.anneal.ln_p_ladders,
                                    samp.anneal.period)

            sampleZ!(samp)
            sampleH!(samp)

            for i in (I)(1):S
                sampleL!(samp, (I)(i), true)
            end

            sampleP_err!(samp)
            for j in (I)(1):M
                sampleS_v!(samp,  (I)(j))
            end

            sampleA!(samp)
            sampleF!(samp)
            sampleG!(samp)
            sampleU!(samp)

            for (c,m) in Iterators.product(keys(samp.usageS), keys(samp.usageV))
                sampleB!(samp,  (I)(c), (I)(m))
            end

            if count > burnin && count % thin == 0
                samp.param.ln_p_v .= ln_p_v_true
                samp.lnProb = ln_P_all(samp)
                push!(lnProbs, samp.lnProb)
                if maxLn < samp.lnProb
                    maxLn = samp.lnProb
                    mapState = deepcopy(samp)
                end
            end
            if count % progressCount == 0
                println(count)
            end
        end
        return (mapState, lnProbs)
    end

    # return the (state of MAP, iterCount, lnProb) with in this iterations
    function sampleAll!(samp::Sampler{I, R};
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
            __sampler.updatePenalty!((I)(count),
                                    samp.param.ln_p_v,
                                    samp.anneal.ln_p_ladders,
                                    samp.anneal.period)
            sampleZ!(samp)
            sampleH!(samp)

            for i in (I)(1):S
                sampleL!(samp, i)
            end

            for j in (I)(1):M
                sampleS_v!(samp, j)
            end

            sampleA!(samp)
            sampleF!(samp)
            sampleG!(samp)
            sampleU!(samp)

            for (c,m) in Iterators.product(keys(samp.usageS), keys(samp.usageV))
                sampleB!(samp, c, m)
            end

            if count > burnin && count % thin == 0
                samp.param.ln_p_v .= ln_p_v_true
                samp.lnProb = ln_P_all(samp)
                now::Sampler{I,R} = deepcopy(samp)
                push!(sampled, (now, count))
            end

            if count % 100 == 0
                println(count)
            end
        end
        # S, M = size(samp.Z)
        # println("time sampleZ")
        # @time for i in 1:1000;  for (i,j) in Iterators.product(1:S,1:M); sampleZ!(samp, i, j); end; end;
        # println("time sampleZ All")
        # @time for i in 1:1000;  sampleZ!(samp); end;
        # println("time sampleH")
        # @time for i in 1:1000;  for (i,j) in Iterators.product(1:S,1:M); sampleH!(samp, i, j); end; end;
        # println("time sampleH All")
        # @time for i in 1:1000;  sampleH!(samp); end;
        #
        # println("time sampleS_s")
        # @time for i in 1:1000; i::I = rand(1:S); sampleS_s!(samp, i); end;
        # println("time sampleS_v")
        # @time for i in 1:1000; j::I = rand(1:M); sampleS_v!(samp, j); end;
        #
        # println("time sampleB")
        # @time for i in 1:1000; for (c, m) in Iterators.product(keys(samp.usageS), keys(samp.usageV) ); sampleB!(samp, c, m); end; end;
        #
        # println("time sampleA")
        # @time for i in 1:1000; for i in 1:S; sampleA!(samp, i); end; end;
        # println("time sampleA All")
        # @time for i in 1:1000; sampleA!(samp); end;
        #
        # println("time sampleF")
        # @time for i in 1:1000; for j in 1:M; sampleF!(samp, j); end; end;
        # println("time sampleF All")
        # @time for i in 1:1000; sampleF!(samp); end;
        # println("time sampleG")
        # @time for i in 1:1000; for j in 1:M; sampleG!(samp, j); end; end;
        # println("time sampleG All")
        # @time for i in 1:1000; sampleG!(samp); end;
        # println("time sampleU")
        # @time for i in 1:1000; for j in 1:M; sampleU!(samp, j); end; end;
        # println("time sampleU All")
        # @time for i in 1:1000; sampleU!(samp); end;
        return sampled
    end

    function init(errScorePath::String, patScorePath::String, matScorePath::String, paramPath::String)
        lnP_D::Array{REAL, 3}   = __sampler.parseData(errScorePath, patScorePath, matScorePath)
        param::Parameters{INT,REAL}, anneal::Annealer{INT, REAL} = inputParser.parseConfigFile(paramPath::String)
        println("=========== data matrix ==========")
        println(lnP_D)
        println("==================================")

        S::INT = size(lnP_D, 1)
        M::INT = size(lnP_D, 2)

        Z::Array{INT, 2}   = convert.(INT, fill(1,S,M)) # init ℤ, Z[i,j] ∈ {1,2}, 1: error, 2: tumor
        H::Array{INT, 2}   = convert.(INT, fill(1,S,M)) # init H, H[i,j] ∈ {1,2}, 1: mat,   2: pat

        D::Array{REAL, 2}  = distanceParser.parseBFHammingDistance(errScorePath, patScorePath, matScorePath, paramPath, "alpha_s")
        decayFunc = distanceParser.parseDecayFunction(paramPath)
        for (i,j) in Iterators.product((INT)(1):S,(INT)(1):S)
            (i != j) && ( D[i,j] = decayFunc(D[i,j]))
        end
        println("f(distance)")
        println(D)
        L::TableGraph{INT, REAL} = tableGraph.init(S, D)

        s_s::Array{INT, 1} = convert.(INT, collect(1:S)) # init sample wise cluster
        s_v::Array{INT, 1} = convert.(INT, collect(1:M)) # init mutation wise cluster
        usageS::Dict{INT,Array{INT,1}} = Dict{INT,Array{INT,1}}()
        usageV::Dict{INT,Array{INT,1}} = Dict{INT,Array{INT,1}}()
        for s in (INT)(1):S; usageS[s] = [s]; end
        for m in (INT)(1):M; usageV[m] = [m]; end

        a::Array{INT, 1} = convert.(INT, fill(1, S))    # init merge sample indicator
        f::Array{REAL,1} = convert.(REAL,fill(0.90, M)) # init mutation freq for each mutation
        g::Array{INT, 1} = convert.(INT, fill(1, M))    # init haplotype for each mutation
        u::Array{INT, 1} = convert.(INT, fill(1, M))    # init mutation to unique sample indicator
        p_err::REAL      = (REAL)(0.05)                      # init false positive mutation rate
        er::Array{INT,1} = convert.(INT, fill(1, M))    # init error indicator for each mutation

        # init block cluster ∈ {1,2,3,4}, 1:shared, 2:meged, 3:unique, 4:error
        B::Dict{Tuple{INT,INT},INT} = Dict{Tuple{INT,INT},INT}()
        for (s,m) in Iterators.product((INT)(1):S,(INT)(1):M); B[(s,m)] = 1; end

        treeCache::BuffPhyloMatrix{INT} = buffPhyloMatrix.init(S, M, bufferSize = M)

        samp::Sampler{INT, REAL} = Sampler{INT, REAL}(Z, H, L, s_s, s_v, usageS, usageV, Set{INT}(), Set{INT}(),
                                                      a, f, g, u, p_err, er, Set{INT}(),  B, treeCache, lnP_D, param, anneal, (REAL)(0.0))
        samp.lnProb = ln_P_all(samp)
        return samp
    end

    function execMAP(errScores::String, matScores::String, patScores::String, iniFile::String;
                     seed::I = (I)(0), iter::I = (I)(100000), thin::I = (I)(10), burnin::I = (I)(10)) where {I <: Integer}
        samp = sampler.init(errScores, matScores, patScores, iniFile)
        map, lnProbs = sampler.sampleMAP!(samp, seed = seed, iter = iter, thin = thin, burnin = burnin)
        return (map, lnProbs)
    end


    function pingSampler(errScores::String, matScores::String, patScores::String, iniFile::String)
        samp = sampler.init(errScores, matScores, patScores, iniFile)
        sampled = sampler.sampleAll!(samp)
        return sampled
    end

    function __pingSampler()
        samp = sampler.init("../../simulationTree/err.score.txt", "../../simulationTree/mat.score.txt",
                           "../../simulationTree/pat.score.txt", "./simpleModel.ini")
        sampled = sampler.sampleAll!(samp)
        return sampled
    end
end

# using sampler
# @time sampler.__pingSampler()
