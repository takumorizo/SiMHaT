Include("config.jl")
Include("buffPhyloMatrix.jl")
Include("sampler.jl")
Include("inputFileParser.jl")

module __result
    using Plots
    using ..config
    using ..sampler
    using ..inputParser

    export sortBiClusteredMatrix
    export viewMatInHeatMap
    export viewFInHeatMap
    export viewUInHeatMap
    export viewAInHeatMap
    export viewBlockTypesInHeatMap
    export getMAPState
    export ACF
    export evaluateFvalueRaw
    export evaluateFvalueFromSampleAns

    function sortBiClusteredMatrix(matrix, usageR, usageC, B)
        cKeys = []
        for c in keys(usageC)
            num = 0
            for r in keys(usageR)
                num += (B[r,c] == 1)
                num -= (B[r,c] == 4)
            end
            append!(cKeys, tuple([c, num]))
        end
        sort!(cKeys, by = x -> x[2], rev = true)
        sortR = []
        sortC = []
        rBreaks = []
        cBreaks = []
        for x in cKeys
            append!(sortC, usageC[x[1]])
            push!(cBreaks, length(usageC[x[1]]))
        end
        for (key, array) in usageR
            append!(sortR, array)
            push!(rBreaks, length(array))
        end
        for i in 1:(length(rBreaks)-1)
            rBreaks[i+1] += rBreaks[i]
        end
        for i in 1:(length(cBreaks)-1)
            cBreaks[i+1] += cBreaks[i]
        end
        return (matrix[sortR, sortC], rBreaks, sortR, cBreaks, sortC)
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

    function viewMatInHeatMap(matrix, rBreaks, sortR, cBreaks, sortC, title;
                              filePath = "", clim = nothing)
        pyplot()
        (clim == nothing) && heatmap(yflip=true, xticks = (1:size(matrix,2),sortC),
                             yticks = (1:size(matrix,1),sortR), matrix, title = title, size=(2000,2000), c=:Blues)
        (clim != nothing) && heatmap(yflip=true, xticks = (1:size(matrix,2),sortC),
                             yticks = (1:size(matrix,1),sortR), matrix, title = title, clim = clim, size=(2000,2000), c=:Blues)
        hline!(rBreaks .+ 0.5)
        vline!(cBreaks .+ 0.5)
        if filePath != ""
            savefig(filePath)
        end
    end

    function viewFInHeatMap(f, usageR, rBreaks, sortR, usageC, cBreaks, sortC, title;
                            filePath = "")
        pyplot()
        matrix = zeros(Float64, length(sortR), length(sortC))
        X, Y = size(matrix)
        R = length(sortR)
        C = length(sortC)
        for y in 1:Y
            for x in 1:X
                matrix[x,y] = f[y]
            end
        end
        matrix = matrix[sortR, sortC]
        heatmap(yflip=true, xticks = (1:size(matrix,2),sortC),
                yticks = (1:size(matrix,1),sortR), clim = (0.7, 1.0), matrix, title = title, size=(2000,2000), c=:Blues)
        hline!(rBreaks .+ 0.5)
        vline!(cBreaks .+ 0.5)
        if filePath != ""
            savefig(filePath)
        end
    end

    function viewUInHeatMap(u, usageR, rBreaks, sortR, usageC, cBreaks, sortC, title;
                            filePath = "")
        pyplot()
        matrix = zeros(Float64, length(sortR), length(sortC))
        X, Y = size(matrix)
        R = length(sortR)
        C = length(sortC)
        for y in 1:Y
            for x in 1:X
                matrix[x,y] = (u[y] == x)
            end
        end
        matrix = matrix[sortR, sortC]
        heatmap(yflip=true, xticks = (1:size(matrix,2),sortC),
                yticks = (1:size(matrix,1),sortR), clim = (0,1), matrix, title = title, size=(2000,2000), c=:Blues)
        hline!(rBreaks .+ 0.5)
        vline!(cBreaks .+ 0.5)
        if filePath != ""
            savefig(filePath)
        end
    end

    function viewAInHeatMap(a, usageR, rBreaks, sortR, usageC, cBreaks, sortC, title;
                            filePath = "")
        pyplot()
        matrix = zeros(Float64, length(sortR), length(sortC))
        X, Y = size(matrix)
        R = length(sortR)
        C = length(sortC)
        for y in 1:Y
            for x in 1:X
                matrix[x,y] = a[x]
            end
        end
        matrix = matrix[sortR, sortC]
        heatmap(yflip=true, xticks = (1:size(matrix,2),sortC),
                yticks = (1:size(matrix,1),sortR), clim = (1,2), matrix, title = title, size=(2000,2000), c=:Blues)
        hline!(rBreaks .+ 0.5)
        vline!(cBreaks .+ 0.5)
        if filePath != ""
            savefig(filePath)
        end
    end

    function viewBlockTypesInHeatMap(B, usageR, rBreaks, sortR, usageC, cBreaks, sortC, title;
                                     filePath = "")
        pyplot()
        matrix = zeros(Float64, length(sortR), length(sortC))
        R = length(sortR)
        C = length(sortC)
        for (r,c) in Iterators.product(keys(usageR),keys(usageC))
            for (i,j) in Iterators.product(usageR[r],usageC[c])
                matrix[i,j] = B[r,c]
            end
        end
        matrix = matrix[sortR, sortC]
        print(matrix)
        heatmap(yflip=true, xticks = (1:size(matrix,2),sortC),
                yticks = (1:size(matrix,1),sortR), clim = (1,4), matrix, title = title, size=(2000,2000), c=:pu_or)
        hline!(rBreaks .+ 0.5)
        vline!(cBreaks .+ 0.5)
        if filePath != ""
            savefig(filePath)
        end
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

    function ACF(s::Array{R, 1})::Array{R, 1} where { R <: Real }
        S = length(s)
        ans::Array{R, 1} = zeros(R, S)
        s_avg::R = sum(s) / S
        s_se::R  = sum( (s .- s_avg) .* (s .- s_avg) ) / (S - 1.0)

        for t in 1:(S-1)
            corr::R = 0.0
            for i in 1:(S-t)
                corr += (s[i] - s_avg) * (s[i+t] - s_avg)
            end
            corr = corr / (S-t)
            ans[t] = corr / s_se
        end
        return ans
    end

    function evaluateFvalueRaw(summaryPath::String, score::Array{REAL, 2}, thres::REAL = (REAL)(0.0))::Tuple{REAL,REAL,REAL}
        summaryTemp::Array{REAL, 2} = inputParser.parseInputSummary(summaryPath)
        summary::Array{REAL, 2}     = zeros(size(score))
        for (x,y) in Iterators.product(1:size(summaryTemp)[1],1:size(summaryTemp)[2])
            summary[x,y] = summaryTemp[x,y]
        end
        TP::INT = sum( (score .>  thres) .* (summary .>  0.0) )
        TN::INT = sum( (score .<= thres) .* (summary .<= 0.0) )
        FP::INT = sum( (score .>  thres) .* (summary .<= 0.0) )
        FN::INT = sum( (score .<= thres) .* (summary .>  0.0) )
        precision::REAL = TP / (TP + FP)
        recall::REAL = TP / (TP + FN)
        print("Precision: "); println(precision);
        print("Recall: ");    println(recall);
        fvalue::REAL = 2.0 * precision * recall / (precision + recall)
        print("Fvalue: ");    println(fvalue);
        return (precision, recall, fvalue)
    end

    function evaluateFvalueFromSampleAns(summaryPath::String, sampled::Sampler{INT, REAL},
                                         thres::REAL = (REAL)(0.0),
                                         rmError::Bool = true)::Tuple{REAL,REAL,REAL}
        score::Array{REAL, 2}       = convert.(REAL, sampled.Z) .- 1.0
        S::INT, M::INT = size(score)
        if rmError
            for (i,j) in Iterators.product(1:S, 1:M)
                (sampled.er[j] == 2) && (score[i,j] = 0.0)
            end
        end
        return evaluateFvalueRaw(summaryPath, score, thres)
    end

end

module result
    using ..__result
    using Plots
    using ..sampler
    # using JLD
    using ..config
    using ..inputParser

    # function writeToJLD(sMax, lnProbs, JLDFilePath)
    #     jldopen(JLDFilePath, "w") do file
    #         addrequire(file, sampler)
    #         write(file, "sMax", sMax)
    #         write(file, "lnProbs", lnProbs)
    #     end
    # end
    #
    # function readFromJLD(JLDFilePath)
    #     sMax = jldopen(JLDFilePath, "r") do file
    #         read(file, "sMax")
    #     end
    #     lnProbs = jldopen(JLDFilePath, "r") do file
    #         read(file, "lnProbs")
    #     end
    #     return sMax, lnProbs
    # end

    function viewMAP(sMax, lnProbs, outputDir::String)
        pyplot()
        iters = collect(1:length(lnProbs))
        plot(iters, lnProbs, size=(2000,2000))
        savefig(outputDir * "/lnProbs.png")

        acf = ACF(lnProbs)
        bar(acf)
        savefig(outputDir * "/acf.png")

        println("isValid state")
        println(sampler.isValid(sMax, debug = true))
        mat, rbreak, sortR, cbreak, sortC = sortBiClusteredMatrix(sMax.Z, sMax.usageS, sMax.usageV, sMax.B)

        # println("MAP, B")
        # println(sMax.B)
        # println("isError valid matrix")
        # println(BuffPhyloMatrix{Int64} <: BuffPhyloMatrix)
        # println(buffPhyloMatrix.BuffPhyloMatrix{Int64} <: BuffPhyloMatrix)
        # println(typeof(sMax.treeCache) <: BuffPhyloMatrix)
        # println(typeof(sMax.treeCache) <: BuffPhyloMatrix{Int64})
        # isTree(sMax.treeCache)

        viewMatInHeatMap(mat, rbreak, sortR, cbreak, sortC, "MAP.Z", filePath = outputDir * "/MAP_Z.png")
        viewMatInHeatMap(sMax.lnPData[sortR,sortC,1], rbreak, sortR, cbreak, sortC, "err", filePath = outputDir * "/err.png")
        viewMatInHeatMap(sMax.lnPData[sortR,sortC,2], rbreak, sortR, cbreak, sortC, "mat", filePath = outputDir * "/mat.png")
        viewMatInHeatMap(sMax.lnPData[sortR,sortC,3], rbreak, sortR, cbreak, sortC, "pat", filePath = outputDir * "/pat.png")

        viewMatInHeatMap(sMax.lnPData[:,:,1], rbreak, sortR, cbreak, sortC, "err", filePath = outputDir * "/err.nonSort.png")
        viewMatInHeatMap(sMax.lnPData[:,:,2], rbreak, sortR, cbreak, sortC, "mat", filePath = outputDir * "/mat.nonSort.png")
        viewMatInHeatMap(sMax.lnPData[:,:,3], rbreak, sortR, cbreak, sortC, "pat", filePath = outputDir * "/pat.nonSort.png")
        viewMatInHeatMap(sMax.H[sortR,sortC], rbreak, sortR, cbreak, sortC, "MAP.H", filePath = outputDir * "/MAP_H.png", clim = [1,2])
        viewFInHeatMap(sMax.f, sMax.usageS, rbreak, sortR, sMax.usageV, cbreak, sortC, "MAP.f", filePath = outputDir * "/MAP_f.png")
        viewUInHeatMap(sMax.u, sMax.usageS, rbreak, sortR, sMax.usageV, cbreak, sortC, "MAP.u", filePath = outputDir * "/MAP_u.png")
        viewAInHeatMap(sMax.a, sMax.usageS, rbreak, sortR, sMax.usageV, cbreak, sortC, "MAP.a", filePath = outputDir * "/MAP_a.png")
        viewBlockTypesInHeatMap(sMax.B, sMax.usageS, rbreak, sortR, sMax.usageV, cbreak, sortC, "MAP.B", filePath = outputDir * "/MAP_B.png")
    end

    function viewDataAll(sampled, outputDir::String)
        pyplot()

        cutoff = 200
        lnProbs = []; iters   = []; for (s,c) in sampled; push!(lnProbs, s.lnProb); push!(iters, c); end; plot(iters[cutoff:length(iters)], lnProbs[cutoff:length(iters)], size=(2000,2000))
        savefig(outputDir * "/lnProbs.png")

        sMax = getMAPState(sampled)
        println("isValid state")
        println(sampler.isValid(sMax, debug = true))
        mat, rbreak, sortR, cbreak, sortC = sortBiClusteredMatrix(sMax.Z, sMax.usageS, sMax.usageV)

        # println("MAP, B")
        # println(sMax.B)
        # println("isError valid matrix")
        # println(BuffPhyloMatrix{Int64} <: BuffPhyloMatrix)
        # println(buffPhyloMatrix.BuffPhyloMatrix{Int64} <: BuffPhyloMatrix)
        # println(typeof(sMax.treeCache) <: BuffPhyloMatrix)
        # println(typeof(sMax.treeCache) <: BuffPhyloMatrix{Int64})
        # isTree(sMax.treeCache)

        viewMatInHeatMap(mat, rbreak, sortR, cbreak, sortC, "MAP.Z", filePath = outputDir * "/MAP_Z.png")
        viewMatInHeatMap(sMax.lnPData[sortR,sortC,1], rbreak, sortR, cbreak, sortC, "err", filePath = outputDir * "/err.png")
        viewMatInHeatMap(sMax.lnPData[sortR,sortC,2], rbreak, sortR, cbreak, sortC, "mat", filePath = outputDir * "/mat.png")
        viewMatInHeatMap(sMax.lnPData[sortR,sortC,3], rbreak, sortR, cbreak, sortC, "pat", filePath = outputDir * "/pat.png")

        viewMatInHeatMap(sMax.H, rbreak, sortR, cbreak, sortC, "MAP.H", filePath = outputDir * "/MAP_H.png", clim = [1,2])
        viewFInHeatMap(sMax.f, sMax.usageS, rbreak, sortR, sMax.usageV, cbreak, sortC, "MAP.f", filePath = outputDir * "/MAP_f.png")
        viewUInHeatMap(sMax.u, sMax.usageS, rbreak, sortR, sMax.usageV, cbreak, sortC, "MAP.u", filePath = outputDir * "/MAP_u.png")
        viewAInHeatMap(sMax.a, sMax.usageS, rbreak, sortR, sMax.usageV, cbreak, sortC, "MAP.a", filePath = outputDir * "/MAP_a.png")
        viewBlockTypesInHeatMap(sMax.B, sMax.usageS, rbreak, sortR, sMax.usageV, cbreak, sortC, "MAP.B", filePath = outputDir * "/MAP_B.png")
    end

    # function evaluateFvalue(summaryPath::String, scorePath::String, thres::REAL = (REAL)(0.0), rmError::Bool = true)::Tuple{REAL,REAL,REAL}
    #     if endswith(scorePath, "jld")
    #         sampled::Sampler{INT, REAL}, lnProbs::Array{REAL, 1} = readFromJLD(scorePath)
    #         return evaluateFvalueFromSampleAns(summaryPath, sampled, thres, rmError)
    #     else
    #         score::Array{REAL, 2}       = inputParser.parseInputSummary(scorePath)
    #         return evaluateFvalueRaw(summaryPath, score, thres)
    #     end
    # end
    #
    # function storeFvalue(outputPath::String, prescision::REAL, recall::REAL, fvalue::REAL, tag::String)
    #     @assert endswith(outputPath, "jld")
    #     if isfile(outputPath)
    #         precisions::Dict{String, REAL} = jldopen(outputPath, "r") do file
    #             read(file, "precisions")
    #         end
    #         recalls::Dict{String, REAL} = jldopen(outputPath, "r") do file
    #             read(file, "recalls")
    #         end
    #         fvalues::Dict{String, REAL} = jldopen(outputPath, "r") do file
    #             read(file, "fvalues")
    #         end
    #         precisions[tag] = prescision
    #         recalls[tag] = recall
    #         fvalues[tag] = fvalue
    #     else
    #         precisions = Dict{String,REAL}(tag=>prescision)
    #         recalls = Dict{String,REAL}(tag=>recall)
    #         fvalues = Dict{String,REAL}(tag=>fvalue)
    #     end
    #     println(precisions)
    #     println(recalls)
    #     println(fvalues)
    #
    #     jldopen(outputPath, "w") do file
    #         write(file, "precisions", precisions)
    #         write(file, "recalls", recalls)
    #         write(file, "fvalues", fvalues)
    #     end
    # end
end
