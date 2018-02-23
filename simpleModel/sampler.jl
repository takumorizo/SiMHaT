include("config64.jl")
include("inputFileParser.jl")
include("phylogenyTree.jl")
include("random.jl")


module __sampler
    using config64
    using inputParser
    using Distributions
    using phylogenyTree

    function parseData(errScorePath::String,
                       patScorePath::String,
                       matScorePath::String)::Array{REAL, 3}
        errScore::Array{REAL, 2} = inputParser.parseInputSummary(errScorePath)
        patScore::Array{REAL, 2} = inputParser.parseInputSummary(patScorePath)
        matScore::Array{REAL, 2} = inputParser.parseInputSummary(matScorePath)
        # param::Parameters{REAL}  = inputParser.parseConfigFile(paramPath)
        S = size(errScore,1)
        M = size(errScore,2)
        lnP::Array{REAL, 3} = convert.(REAL, zeros(S,M,3))
        for (s,m) in Iterators.product(1:S,1:M)
            lnP[s,m,1] = errScore[s,m]
            lnP[s,m,2] = patScore[s,m]
            lnP[s,m,3] = matScore[s,m]
        end
        return lnP
    end

    function areBlocksTree(B::Dict{Tuple{I,I}, I},
                          c::I, m::I,
                          usageS::Dict{I,Array{I,1}},
                          usageV::Dict{I,Array{I,1}},
                          positiveSet::Set{I};
                          seeDriver::Bool = false,
                          minX::I = 1,
                          minY::I = 1)::Array{Bool, 1} where {I <: Integer}
       ans::Array{Bool,1} = [true, true]
       C::I                = length(usageS)
       M::I                = length(usageV)
       blocks::Array{I, 2} = zeros(I, C, M)
       Cs::Array{I,1}      = collect(keys(usageS))
       Ms::Array{I,1}      = collect(keys(usageV))
       for (cc,mm) in Iterators.product(1:C,1:M)
           cluster_c::I = Cs[cc];
           cluster_m::I = Ms[mm];
           blocks[cc,mm] = (I)(B[cluster_c, cluster_m] ∈ positiveSet)
       end
       cIdx::I = findfirst(Cs, c)
       mIdx::I = findfirst(Ms, m)
       blocks[cIdx,mIdx] = (I)(0)
       ans[1] = phylogenyTree.isPhylogenic(blocks, checkDriver = seeDriver, minX = minX, minY = minY)
       blocks[cIdx,mIdx] = (I)(1)
       ans[2] = phylogenyTree.isPhylogenic(blocks, checkDriver = seeDriver, minX = minX, minY = minY)
       return ans
    end

    function isBlocksTree(B::Dict{Tuple{I,I}, I},
                         usageS::Dict{I,Array{I,1}},
                         usageV::Dict{I,Array{I,1}},
                         positiveSet::Set{I};
                         seeDriver::Bool = false,
                         minX::I = 1,
                         minY::I = 1)::Bool where {I <: Integer}
       C::I = length(usageS)
       M::I = length(usageV)
       blocks::Array{I, 2} = zeros(I, C, M)
       Cs::Array{I,1}      = collect(keys(usageS))
       Ms::Array{I,1}      = collect(keys(usageV))
       for (cc,mm) in Iterators.product(1:C,1:M)
           cluster_c::I = Cs[cc];
           cluster_m::I = Ms[mm];
           blocks[cc,mm] = (I)(B[cluster_c, cluster_m] ∈ positiveSet)
       end
       return phylogenyTree.isPhylogenic(blocks, checkDriver = seeDriver, minX = minX, minY = minY)
    end

    function isBlockValid(needsTree::Bool,
                          Z::Array{I, 2},
                          usageS::Dict{I,Array{I,1}},
                          usageV::Dict{I,Array{I,1}},
                          c::I, m::I;
                          offset::I = 1,
                          seeDriver::Bool = true,
                          minX::I = 1,
                          minY::I = 1)::Bool where {I <: Integer}
      return (!needsTree) || (phylogenyTree.isPhylogenic(view(Z, usageS[c], usageV[m]),
                                                         offset = offset, checkDriver = seeDriver,
                                                         minX = minX, minY = minY))
    end

    function lnAF(α::R, n::I)::R where {R <: Real, I <: Integer}
        ans = (R)(0.0)
        for i in 0:(n-1); ans += log(e, α + i); end
        return ans
    end

    function ln_P_CRP(usage::Dict{I,Array{I,1}}, α::R)::R where {R <: Real, I <: Integer}
        ans::R = (R)(0.0)
        ans += (length(keys(usage))) * log(e, α)
        totalNum::I = 0
        for (cluster, indexes) in usage
            number::I = (I)(length(indexes))
            ans += lgamma(number)
            totalNum += number
        end
        ans -= lnAF(α, totalNum)
        return ans
    end

    # x ∈ {0,1}, x = 1 w.p. p
    function ln_P_ber(x::I, p::R)::R where {R <: Real, I <: Integer}
        return x * log(e, p) + (1-x) * log(e, (1-p))
    end

    function ln_P_beta(p::R, α::R, β::R)::R where {R <: Real}
        return (α-(R)(1.0))*log(e, p) + (β-(R)(1.0))*log(e,((R)(1.0)-p)) - lbeta(α, β)
    end

    function ln_P_Z(Z::Array{I, 2},
                   S_s::Array{I, 1},
                   S_v::Array{I, 1},
                   f::Array{R, 3})::R where {R <: Real, I <: Integer}
        II::I = size(Z,1); J::I = size(Z,2)
        ans::R = (R)(0.0)
        for (i,j) in Iterators.product(1:II,1:J)
            ans += ln_P_ber(Z[i,j]-1, f[ S_s[i], S_v[j], j ])
        end
        return ans
    end

    function ln_P_Data_Z(lnPData::Array{R, 3},
                          Z::Array{I, 2})::R where {I <: Integer, R <: Real}
        II::I = (I)(size(Z,1))
        J::I  = (I)(size(Z,2))
        ans::R = (R)(0.0)
        for (i,j) in Iterators.product(1:II,1:J)
            # ans += lnPData[i,j, 1 + (Z[i,j]-1) * Y[i,j] ]
            ans += lnPData[i,j,Z[i,j]]
        end
        return ans
    end

    function ln_P_B(B::Dict{Tuple{I,I},I},
                   λ::Array{R, 1})::R where {I <: Integer, R <: Real}
        ans::R = 0.0
        for (key, val) in B; ans += log(e, λ[val]); end
        return ans
    end

    function ln_P_f_B(f::Array{R, 3},
                     B::Dict{Tuple{I,I},I},
                     βs::Array{R,2})::R where {I <: Integer, R <: Real}
        J::I = size(f,3)
        ans::R = (R)(0.0)
        for (cluster, value) in B
            c = cluster[1]
            m = cluster[2]
            for j in 1:J
                ans += ln_P_beta(f[c,m,j], βs[value,1], βs[value,2])
            end
        end
        return ans
    end

    function ln_P_t(Z::Array{I, 2},
                    usageS::Dict{I,Array{I,1}},
                    usageV::Dict{I,Array{I,1}},
                    B::Dict{Tuple{I,I},I},
                    ln_g::R, ln_1_m_g::R,
                    minX = 1,
                    minY = 1)::R where{I <: Integer, R <: Real}
        ans::R = (R)(0.0)
        for (cluster2D, value) in B
            (c::I, m::I) = cluster2D
            isTree::Bool = __sampler.isBlockValid(value == 3, Z, usageS, usageV, c, m,
                                        offset = 1, seeDriver = true, minX = minX, minY = minY)
            ans += (R)(isTree) * ln_g + (1.0 - (R)(isTree)) * ln_1_m_g
        end
        return ans
    end

    function ln_P_w_B(B::Dict{Tuple{I,I},I},
                     usageS::Dict{I,Array{I,1}},
                     usageV::Dict{I,Array{I,1}},
                     ln_g::R, ln_1_m_g::R,
                     minX = 1,
                     minY = 1)::R where {I <: Integer, R <: Real}
        ans::R = (R)(0.0)
        isTree::Bool = isBlocksTree(B,usageS, usageV, Set([3]), seeDriver = false, minX = minX, minY = minY)
        return (R)(isTree) * ln_g + (1.0 - (R)(isTree)) * ln_1_m_g
    end

    function exp_normalize!(ln_p::Array{R,1})::Void where {R <: Real}
        ln_p .= (ln_p .- maximum(ln_p))
        ln_p .= exp.(ln_p)
        ln_p .= ln_p ./ sum(ln_p)
        return nothing
    end

    function log_sum_exp(ln_p::Array{R,1})::R where {R <: Real}
        maxVal::R = maximum(ln_p)
        return log(e, sum(exp.( (ln_p .- maxVal) ))) + maxVal
    end


end

module sampler
    using config64
    using inputParser
    using __sampler
    using random
    using phylogenyTree

    export Sampler
    type Sampler{I,R}
        Z::Array{I, 2}
        S_s::Array{I, 1}
        S_v::Array{I, 1}
        f::Array{R, 3}
        B::Dict{Tuple{I,I}, I}
        usageS::Dict{I,Array{I,1}}
        unUsedS::Array{I, 1}
        usageV::Dict{I,Array{I,1}}
        unUsedV::Array{I, 1}
        lnPData::Array{R, 3}
        param::Parameters{I,R}
        lnProb::R
        # validity data for debug
        # T::Dict{Tuple{I,I},Bool}
        # W::Bool
    end

    function deepcopySampler!(from::Sampler{I, R}, to::Sampler{I, R})::Void where {I<:Integer, R <: Real}
        to.Z = deepcopy(from.Z)
        to.S_s = deepcopy(from.S_s)
        to.S_v = deepcopy(from.S_v)
        to.f = deepcopy(from.f)
        to.B = deepcopy(from.B)
        to.usageS = deepcopy(from.usageS)
        to.unUsedS = deepcopy(from.unUsedS)
        to.usageV = deepcopy(from.usageV)
        to.unUsedV = deepcopy(from.unUsedV)
        to.lnPData = deepcopy(from.lnPData)
        to.param = deepcopy(from.param)
        to.lnProb = deepcopy(from.lnProb)
        return nothing
    end

    # # # return lnProb after the sampling state
    function sampleZ!(samp::Sampler{I, R},
                      i::I, j::I)::Void where {I<:Integer, R <: Real}
        now_Z::I = samp.Z[i,j]
        ln_p::Array{R,1} = convert.(R, [0.0, 0.0])
        c::I, m::I = samp.S_s[i], samp.S_v[j]

        # data generation probability part
        ln_p[1] += samp.lnPData[i,j,1]
        ln_p[2] += samp.lnPData[i,j,2]

        # tree probability penalty part
        isValid::Array{Bool, 1} = [true, true]
        samp.Z[i,j] = (I)(1)
        isValid[1] = __sampler.isBlockValid(samp.B[c,m] == 3,
                     samp.Z, samp.usageS, samp.usageV, c, m,
                     seeDriver = true,
                     minX = samp.param.minXIn,
                     minY = samp.param.minYIn)
        samp.Z[i,j] = (I)(2)
        isValid[2] = __sampler.isBlockValid(samp.B[c,m] == 3,
                     samp.Z, samp.usageS, samp.usageV, c, m,
                     seeDriver = true,
                     minX = samp.param.minXIn,
                     minY = samp.param.minYIn)

        ln_p[1] += (R)(isValid[1]) * samp.param.ln_G_t[1] + (1.0-(R)(isValid[1])) * samp.param.ln_G_t[2]
        ln_p[2] += (R)(isValid[2]) * samp.param.ln_G_t[1] + (1.0-(R)(isValid[2])) * samp.param.ln_G_t[2]

        # prior prob part
        ln_p[1] += log(e, (R)(1.0)-samp.f[c,m,j])
        ln_p[2] += log(e, samp.f[c,m,j])

        __sampler.exp_normalize!(ln_p)
        new_Z::I = indmax( random.sampleMultiNomial(1, ln_p) )
        samp.Z[i,j] = new_Z

        return nothing
    end

    function sampleF!(samp::Sampler{I, R}, c::I, m::I, j::I)::Void where {I <: Integer, R <: Real}
        prev::R = samp.lnProb
        now_F::R = samp.f[c,m,j]
        β::Array{R,1} = samp.param.βs[ samp.B[c,m], : ]
        β[1] += sum(samp.Z[samp.usageS[c],j] .== 2)
        β[2] += sum(samp.Z[samp.usageS[c],j] .== 1)
        new_F::R = random.sampleBeta( β[1], β[2] )
        samp.f[c,m,j] = new_F
        return nothing
    end


    function sampleB!(samp::Sampler{I, R}, c::I, m::I)::Void where {I <: Integer, R <: Real}
        now_B::I = samp.B[c,m]
        ln_p::Array{R, 1} = convert.(R, log.(e, samp.param.λs))
        J::I = size(samp.Z, 2)
        # f_[c,m,j] lnprob part
        for b in 1:3
            for j in 1:J
                ln_p[b] -= lbeta(samp.param.βs[b,1], samp.param.βs[b,2])
                ln_p[b] += (R)(samp.param.βs[b,1])*(R)(log(e, (R)(samp.f[c,m,j])))
                ln_p[b] += (R)(samp.param.βs[b,2])*(R)(log(e, (R)(1.0-samp.f[c,m,j])))
            end
        end
        # Block tree lnProb part
        isBlocksTree::Array{Bool,1} = [true, true, true]
        isBlocksTree[[1,3]] = __sampler.areBlocksTree(samp.B, c, m, samp.usageS, samp.usageV, Set([3]),
                                                      minX = samp.param.minXOut,
                                                      minY = samp.param.minYOut)
        isBlocksTree[2] = isBlocksTree[1]

        for b in 1:3
            ln_p[b] += (R)(isBlocksTree[b]) * (R)(samp.param.ln_G_B[1])
            ln_p[b] += ((R)(1.0)-(R)(isBlocksTree[b])) * (R)(samp.param.ln_G_B[2])
        end
        isInnerBlockValid::Bool = __sampler.isBlockValid(samp.B[c,m] == 3,
                                               samp.Z, samp.usageS, samp.usageV, c, m,
                                               seeDriver = true,
                                               minX = samp.param.minXIn,
                                               minY = samp.param.minYIn)
        ln_p[3] += (R)(isInnerBlockValid) * (R)(samp.param.ln_G_t[1])
        ln_p[3] += ((R)(1.0)-(R)(isInnerBlockValid)) * (R)(samp.param.ln_G_t[2])

        __sampler.exp_normalize!(ln_p)
        new_B::I = indmax(random.sampleMultiNomial(1,ln_p))
        samp.B[c,m] = new_B

        return nothing
    end

    function getAcceptRatioRow(prev::Sampler{I, R}, post::Sampler{I,R}, i::I) where {I <: Integer, R <: Real}
        a1::R = 0.0
        new_S::I = post.S_s[i]
        # Z part
        for m in keys(post.usageV)
            for j in post.usageV[m]
                a1 += __sampler.ln_P_ber(post.Z[i,j]-1, post.f[ post.S_s[i], post.S_v[j], j ])
            end
        end
        # internal tree part
        for m in keys(post.usageV)
            isValidPost::Bool = __sampler.isBlockValid(post.B[new_S,m] == 3,
                                             post.Z, post.usageS, post.usageV, new_S, m,
                                             seeDriver = true,
                                             minX = post.param.minXIn,
                                             minY = post.param.minYIn)
            a1 += (R)(isValidPost) * post.param.ln_G_t[1] + (1.0-(R)(isValidPost)) * post.param.ln_G_t[2]
        end
        # block tree part
        isBlocksTree::Bool = __sampler.isBlocksTree(post.B, post.usageS, post.usageV, Set([3]),
                                                   minX = post.param.minXOut,
                                                   minY = post.param.minYOut)
        a1 += (R)(isBlocksTree) * post.param.ln_G_B[1] + (1.0-(R)(isBlocksTree)) * post.param.ln_G_B[2]

        a2::R = 0.0
        now_S::I = prev.S_s[i]
        # Z part
        for m in keys(prev.usageV)
            for j in prev.usageV[m]
                a2 += __sampler.ln_P_ber(prev.Z[i,j]-1, prev.f[ prev.S_s[i], prev.S_v[j], j ])
            end
        end
        # internal tree part
        for m in keys(prev.usageV)
            isValidPrev::Bool = __sampler.isBlockValid(prev.B[now_S,m] == 3,
                                             prev.Z, prev.usageS, prev.usageV, now_S, m,
                                             seeDriver = true,
                                             minX = prev.param.minXIn,
                                             minY = prev.param.minYIn)
            a2 += (R)(isValidPrev) * prev.param.ln_G_t[1] + (1.0-(R)(isValidPrev)) * prev.param.ln_G_t[2]
        end
        # block tree part
        isBlocksTree = __sampler.isBlocksTree(prev.B, prev.usageS, prev.usageV, Set([3]),
                                              minX = prev.param.minXOut,
                                              minY = prev.param.minYOut)
        a2 += (R)(isBlocksTree) * prev.param.ln_G_B[1] + (1.0-(R)(isBlocksTree)) * prev.param.ln_G_B[2]

        return (R)(min(1.0, exp(a1-a2)))
    end

    function proposeNewClusterRow(samp::Sampler{I, R}, i::I)::Sampler{I, R} where{ I <: Integer, R <: Real }
        J::I     = size(samp.Z, 2)
        now_S::I = samp.S_s[i]
        propose::Sampler{I, R} = deepcopy(samp)

        deleteat!(propose.usageS[now_S], findin(propose.usageS[now_S], i))
        if length(propose.usageS[now_S]) == 0
            delete!(propose.usageS, now_S)
            push!(propose.unUsedS,  now_S)
            for m in keys(propose.usageV)
                delete!(propose.B, (now_S, m))
            end
        end
        candidates::Array{I, 1} = deepcopy(collect(keys(propose.usageS)))
        push!(candidates, pop!(propose.unUsedS))

        ln_p::Array{R, 1} = zeros(R, length(candidates))
        for c in 1:(length(candidates)-1)
            cluster::I = candidates[c]
            sizeAt::I  = length(propose.usageS[cluster])
            ln_p[c]    = (R)(log(e, (R)(sizeAt)))
        end
        ln_p[length(candidates)] = (R)(log(e, propose.param.α_s))

        __sampler.exp_normalize!(ln_p)
        new_S_Idx::I = indmax( random.sampleMultiNomial(1, ln_p) )
        new_S::I = candidates[new_S_Idx]

        # update data structure
        propose.S_s[i] = new_S
        if in(new_S, keys(propose.usageS))
            push!(propose.usageS[new_S], i)
            push!(propose.unUsedS, new_S)
        else
            propose.usageS[new_S] = [i]
            for m in keys(propose.usageV)
                b::I = indmax( random.sampleMultiNomial(1, propose.param.λs) )
                propose.B[new_S, m] = b
                for j in 1:J
                    propose.f[new_S,m,j] = random.sampleBeta(propose.param.βs[b,1], propose.param.βs[b,2])
                end
            end
        end
        return propose
    end
    # sample row cluster following the metropolis update rules
    function sampleS_s!(samp::Sampler{I, R}, i::I)::Void where{ I <: Integer, R <: Real}
        propose::Sampler{I, R} = proposeNewClusterRow(samp, i)
        accRatio::R            = getAcceptRatioRow(samp, propose, i)
        # print("now      (col, usageS): "); print(i); print(","); println(samp.usageS);
        # print("proposal (col, usageS): "); print(i); print(","); println(propose.usageS);
        # print("accRatio: "); println(accRatio)
        if rand(R) < accRatio
            # accept the proposal
            deepcopySampler!(propose, samp)
            # println("accepted")
        # else
            # println("rejected")
        end
        return nothing
    end

    function getAcceptRatioCol(prev::Sampler{I, R}, post::Sampler{I,R}, j::I) where {I <: Integer, R <: Real}
        b1::R = 0.0
        new_V::I = post.S_v[j]
        # Z part
        for c in keys(post.usageS)
            for i in post.usageS[c]
                b1 += __sampler.ln_P_ber(post.Z[i,j]-1, post.f[ post.S_s[i], post.S_v[j], j ])
            end
        end
        # internal tree part
        for c in keys(post.usageS)
            isValidPost::Bool = __sampler.isBlockValid(post.B[c, new_V] == 3,
                                             post.Z, post.usageS, post.usageV, c, new_V,
                                             seeDriver = true,
                                             minX = post.param.minXIn,
                                             minY = post.param.minYIn)
            b1 += (R)(isValidPost) * post.param.ln_G_t[1] + (1.0-(R)(isValidPost)) * post.param.ln_G_t[2]
        end
        # block tree part
        isBlocksTree::Bool = __sampler.isBlocksTree(post.B, post.usageS, post.usageV, Set([3]),
                                                    minX = post.param.minXOut,
                                                    minY = post.param.minYOut)
        b1 += (R)(isBlocksTree) * post.param.ln_G_B[1] + (1.0-(R)(isBlocksTree)) * post.param.ln_G_B[2]

        b2::R = 0.0
        now_V::I = prev.S_v[j]
        # Z part
        for c in keys(prev.usageS)
            for i in prev.usageS[c]
                b2 += __sampler.ln_P_ber(prev.Z[i,j]-1, prev.f[ prev.S_s[i], prev.S_v[j], j ])
            end
        end
        # internal tree part
        for c in keys(prev.usageS)
            isValidPrev::Bool = __sampler.isBlockValid(prev.B[c, now_V] == 3,
                                             prev.Z, prev.usageS, prev.usageV, c, now_V,
                                             seeDriver = true,
                                             minX = prev.param.minXIn,
                                             minY = prev.param.minYIn)
            b2 += (R)(isValidPrev) * prev.param.ln_G_t[1] + (1.0-(R)(isValidPrev)) * prev.param.ln_G_t[2]
        end
        # block tree part
        isBlocksTree = __sampler.isBlocksTree(prev.B, prev.usageS, prev.usageV, Set([3]),
                                              minX = prev.param.minXOut,
                                              minY = prev.param.minYOut)
        b2 += (R)(isBlocksTree) * prev.param.ln_G_B[1] + (1.0-(R)(isBlocksTree)) * prev.param.ln_G_B[2]

        return (R)(min(1.0, exp(b1-b2)))
    end


    function proposeNewClusterCol(samp::Sampler{I, R}, j::I)::Sampler{I, R} where{ I <: Integer, R <: Real }
        J::I = size(samp.Z, 2)
        now_V::I = samp.S_v[j]
        propose::Sampler{I, R} = deepcopy(samp)

        deleteat!(propose.usageV[now_V], findin(propose.usageV[now_V], j))
        if length(propose.usageV[now_V]) == 0
            delete!(propose.usageV, now_V)
            push!(propose.unUsedV,  now_V)
            for c in keys(propose.usageS)
                delete!(propose.B, (c, now_V))
            end
        end
        candidates::Array{I, 1} = deepcopy(collect(keys(propose.usageV)))
        push!(candidates, pop!(propose.unUsedV))

        ln_p::Array{R, 1} = zeros(R, length(candidates))
        for m in 1:(length(candidates)-1)
            cluster::I = candidates[m]
            sizeAt::I  = length(propose.usageV[cluster])
            ln_p[m]    = (R)(log(e, (R)(sizeAt)))
        end
        ln_p[length(candidates)] = (R)(log(e, propose.param.α_v))

        __sampler.exp_normalize!(ln_p)
        new_V_Idx::I = indmax( random.sampleMultiNomial(1, ln_p) )
        new_V::I = candidates[new_V_Idx]

        # update data structure
        propose.S_v[j] = new_V
        if in(new_V, keys(propose.usageV))
            push!(propose.usageV[new_V], j)
            push!(propose.unUsedV, new_V)
        else
            propose.usageV[new_V] = [j]
            for c in keys(propose.usageS)
                b::I = indmax( random.sampleMultiNomial(1, propose.param.λs) )
                propose.B[c, new_V] = b
                for j in 1:J
                    propose.f[c,new_V,j] = random.sampleBeta(propose.param.βs[b,1], propose.param.βs[b,2])
                end
            end
        end
        return propose
    end

    function sampleV_s!(samp::Sampler{I, R}, j::I)::Void where{ I <: Integer, R <: Real}
        propose::Sampler{I, R} = proposeNewClusterCol(samp, j)
        accRatio::R            = getAcceptRatioCol(samp, propose, j)
        # print("now      (row, usageV): "); print(j); print(","); println(samp.usageV); println(samp.S_v)
        # print("proposal (row, usageV): "); print(j); print(","); println(propose.usageV); println(propose.S_v)
        # print("accRatio: "); println(accRatio)
        if rand(R) < accRatio
            # accept the proposal
            deepcopySampler!(propose, samp)
            # samp = deepcopy(propose)
            # println("accepted")
        # else
            # println("rejected")
        end
        return nothing
    end

    # # return the (state of MAP, iterCount, lnProb) with in this iterations
    function sampleAll!(samp::Sampler{I, R},
                        seed::I = 0,
                        iter::I = 100000,
                        thin::I = 100,
                        burnin::I = 0)::Array{Tuple{Sampler{I, R}, I}, 1} where {I <:Integer, R <: Real }
        # setting the given random seed
        srand(seed)
        sampled::Array{ Tuple{Sampler{I,R}, I}, 1} = []
        if samp.lnProb > 0.0
            error("invalid prob")
        end
        for count in 1:(iter+burnin)
            S::I, M::I = size(samp.Z)
            for (i,j) in Iterators.product(1:S,1:M)
                sampleZ!(samp, i, j)
            end
            for i in 1:S
                sampleS_s!(samp, i)
            end
            for j in 1:M
                sampleV_s!(samp, j)
            end
            for (c,m) in Iterators.product(keys(samp.usageS), keys(samp.usageV))
                sampleB!(samp, c, m)
            end

            for (c,m) in Iterators.product(keys(samp.usageS), keys(samp.usageV))
                for j in samp.usageV[m]
                    sampleF!(samp, c, m, j)
                end
            end

            if count > burnin && count % thin == 0
                samp.lnProb = ln_P_all(samp)
                now::Sampler{I,R} = deepcopy(samp)
                push!(sampled, (now, count))
            end
        end
        return sampled
    end

    function ln_P_all(samp::Sampler{I,R})::R where {I<:Integer, R<:Real}
        ans::R = 0.0
        ans += __sampler.ln_P_CRP(samp.usageS, samp.param.α_s)
        ans += __sampler.ln_P_CRP(samp.usageV, samp.param.α_v)
        ans += __sampler.ln_P_Z(samp.Z, samp.S_s, samp.S_v, samp.f)
        ans += __sampler.ln_P_Data_Z(samp.lnPData, samp.Z)
        ans += __sampler.ln_P_B(samp.B, samp.param.λs)
        ans += __sampler.ln_P_f_B(samp.f, samp.B, samp.param.βs)
        ans += __sampler.ln_P_t(samp.Z, samp.usageS, samp.usageV, samp.B,
                                samp.param.ln_G_t[1], samp.param.ln_G_t[2],
                                samp.param.minXIn, samp.param.minYIn)
        ans += __sampler.ln_P_w_B(samp.B, samp.usageS, samp.usageV,
                                  samp.param.ln_G_B[1], samp.param.ln_G_B[2],
                                  samp.param.minXOut, samp.param.minYOut)
        return ans
    end

    function init(errScorePath::String, patScorePath::String, matScorePath::String, paramPath::String)
        lnP_D::Array{REAL, 3}   = __sampler.parseData(errScorePath, patScorePath, matScorePath)
        param::Parameters{INT,REAL} = inputParser.parseConfigFile(paramPath::String)
        println("=========== data matrix ==========")
        println(lnP_D)
        println("==================================")

        S = size(lnP_D, 1)
        M = size(lnP_D, 2)
        # init ℤ, Z[i,j] ∈ {1,2}, 1: error, 2: tumor
        Z::Array{INT, 2}   = convert.(INT, fill(1,S,M))
        # init sample/mutation wise cluster
        s_s::Array{INT, 1} = convert.(INT, collect(1:S))
        s_v::Array{INT, 1} = convert.(INT, collect(1:M))
        # init mutation freq for each mutation @ eaach cluster
        f::Array{REAL, 3}  = convert.(REAL,fill(0.50, S, M, M))
        # init block cluster ∈ {1,2,3}, 1:error, 2:merged, 3:tree
        B::Dict{Tuple{INT,INT},INT} = Dict{Tuple{INT,INT},INT}()

        for (s,m) in Iterators.product(1:S,1:M); B[(s,m)] = 1; end
        # init
        usageS::Dict{INT,Array{INT,1}} = Dict{INT,Array{INT,1}}()
        usageV::Dict{INT,Array{INT,1}} = Dict{INT,Array{INT,1}}()

        for s in 1:S; usageS[s] = [s]; end
        for m in 1:M; usageV[m] = [m]; end

        samp::Sampler{INT,REAL} = Sampler{INT,REAL}(Z, s_s, s_v, f, B, usageS, [], usageV, [], lnP_D, param, 0.0)
        samp.lnProb = ln_P_all(samp)
        return samp
    end

    function sortBiClusteredMatrix(matrix, usageR, usageC)
        sortR = []
        sortC = []
        rBreaks = []
        cBreaks = []
        for (key, array) in usageR
            append!(sortR, array)
            push!(rBreaks, length(array))
        end
        for (key, array) in usageC
            append!(sortC, array)
            push!(cBreaks, length(array))
        end
        for i in 1:(length(rBreaks)-1)
            rBreaks[i+1] += rBreaks[i]
        end
        for i in 1:(length(cBreaks)-1)
            cBreaks[i+1] += cBreaks[i]
        end
        return (matrix[sortR, sortC], rBreaks, sortR, cBreaks, sortC)
    end

    function viewMatInHeatMap(matrix, rBreaks, sortR, cBreaks, sortC, title)
        heatmap(yflip=true, xticks = (1:size(matrix,2),sortC),
                yticks = (1:size(matrix,1),sortR), matrix, aspect_ratio = 1, title = title)
        hline!(rBreaks .+ 0.5)
        vline!(cBreaks .+ 0.5)
    end

    function getMAPState(sampled)
        maxLik     = -Inf
        maxCountAt = 0
        for (s, count) in sampled
            if s.lnProb > maxLik
                maxCountAt = count
                maxLik = s.lnProb
            end
        end
        sMax = nothing
        for (s, count) in sampled;
            if count == maxCountAt; sMax = s; end
        end
        return sMax
    end

    function __pingSampler()
        samp = sampler.init("../simulationTree/err.score.txt", "../simulationTree/pat.score.txt",
                    "../simulationTree/mat.score.txt", "./simpleModel.ini")
        sampled = sampler.sampleAll!(samp)
        return sampled
        ## averaging Z value
        """
        using Plots
        pyplot()

        lnProbs = []; iters   = []; for (s,c) in sampled; push!(lnProbs, s.lnProb); push!(iters,   c); end; plot(iters, lnProbs)

        sMax = getMAPState(sampled)
        mat, rbreak, sortR, cbreak, sortC = sortBiClusteredMatrix(sMax.Z, sMax.usageS, sMax.usageV)
        viewMatInHeatMap(mat, rbreak, sortR, cbreak, sortC, "MAP.Z")

        errData = samp.lnPData[:,:,1]
        patData = samp.lnPData[:,:,2]
        matData = samp.lnPData[:,:,3]
        heatmap(yflip=true, errData, aspect_ratio = 1, title = "errData")
        savefig("./debug/errData.png")
        heatmap(yflip=true, patData, aspect_ratio = 1, title = "patData")
        savefig("./debug/patData.png")
        heatmap(yflip=true, matData, aspect_ratio = 1, title = "matData")
        savefig("./debug/matData.png")

        z_avg = zeros(Float64, size(samp.Z) )
        for (s, count) in sampled
            z_avg .= (z_avg .+ convert.(Float64, s.Z) )
        end
        z_avg .= z_avg ./ length(sampled)
        z_avg .= z_avg .- 1.0

        heatmap(yflip=true, z_avg, aspect_ratio = 1, title = "z_avg")
        savefig("./debug/z_avg.png")

        avg_Cols = 0.0
        avg_Rows = 0.0
        for (s, count) in sampled
            avg_Cols += (Float64)(length(s.usageS))
            avg_Rows += (Float64)(length(s.usageV))
        end
        avg_Cols /= length(sampled)
        avg_Rows /= length(sampled)

        heatmap(yflip=true, sMax.Z, aspect_ratio = 1, title = "m.l.e Z")

        for (key,value) in sMax.B
            print("value: "); print(value); print(", usageS: "); print(sMax.usageS[key[1]]);
            print(", usageV: "); println(sMax.usageV[key[2]]);
        end

        """
        # return nothing

        # sampleAll!(samp::Sampler{I, R},
        #                     seed::I = 0,
        #                     iter::I = 2000,
        #                     thin::I = 4,
        #                     burnin::I = 1000)::Array{Tuple{Sampler{I, R}, I}, 1} where {I <:Integer, R <: Real }
    end
end
