include("buffPhyloMatrix.jl")
include("sampler.jl")

module __result
    using Plots
    export sortBiClusteredMatrix
    export viewMatInHeatMap
    export viewFInHeatMap
    export viewUInHeatMap
    export viewAInHeatMap
    export viewBlockTypesInHeatMap
    export getMAPState

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
                             yticks = (1:size(matrix,1),sortR), matrix, title = title, size=(2000,2000), c=:ice)
        (clim != nothing) && heatmap(yflip=true, xticks = (1:size(matrix,2),sortC),
                             yticks = (1:size(matrix,1),sortR), matrix, title = title, clim = clim, size=(2000,2000), c=:ice)
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
                yticks = (1:size(matrix,1),sortR), clim = (0.7, 1.0), matrix, title = title, size=(2000,2000), c=:ice)
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
                yticks = (1:size(matrix,1),sortR), clim = (0,1), matrix, title = title, size=(2000,2000), c=:ice)
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
                yticks = (1:size(matrix,1),sortR), clim = (1,2), matrix, title = title, size=(2000,2000), c=:ice)
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

end

module result
    using __result
    using Plots
    using sampler

    function viewMAP(sMax, lnProbs, outputDir::String)
        pyplot()
        iters = collect(1:length(lnProbs))
        plot(iters, lnProbs, size=(2000,2000))
        savefig(outputDir * "/lnProbs.png")

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

    function viewDataAll(sampled, outputDir::String)
        pyplot()

        cutoff = 100
        lnProbs = []; iters   = []; for (s,c) in sampled; push!(lnProbs, s.lnProb); push!(iters,   c); end; plot(iters[cutoff:length(iters)], lnProbs[cutoff:length(iters)], size=(2000,2000))
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
end

using sampler
using result
using DocOpt
using Plots
pyplot()

nowTime = (Dates.value(Dates.now())) #Int(Dates.now())
cwd = pwd()
doc = """SubsetErrorModel

Usage:
    result.jl MAP <errScore> <matScore> <patScore> <iniFile> [options]

Options:
  -o <DIR> --outDir=<DIR>  Output directory. We use current working directory if unspecified. [default: $cwd]
  -s <SEED> --seed=<SEED>  Seed of randomness. We use current time (msec) if unspecified. [default: $nowTime]
  -n <NUM> --number=<NUM>  Number of iteration after burnin. [default: 1000]
  -b <BURNIN> --burnin=<BURNIN>  Number of iteration during burnin. [default: 0]
  -t <THIN> --thin=<THIN>  Duration between sampling. [default: 1]
  -h --help  Show this screen.

"""
# println(doc)
args = docopt(doc)
println(args)
@time sMax, lnProbs  = sampler.execMAP(args["<errScore>"], args["<matScore>"], args["<patScore>"], args["<iniFile>"],
                                       seed = parse(args["--seed"]), iter = parse(args["--number"]), thin = parse(args["--thin"]), burnin = parse(args["--burnin"]))
result.viewMAP(sMax, lnProbs, args["--outDir"])
