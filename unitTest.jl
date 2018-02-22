include("phylogenyTree.jl")

module __phylogenyTreeTest
    function makeTreeMatrixRandomly!(mat::Array{Int, 2}, rangeX, rangeY, p = 0.5, divRange=1:3)::Void
        for (x,y) in Iterators.product(rangeX,rangeY)
             mat[x,y] = 0
        end
        if size(rangeX, 1) <= 1 || size(rangeY, 1) <= 1 || rand() <= p
            iter = 0
            for y in rangeY
                for x in rangeX
                    if x <= rangeX.stop - iter
                        mat[x,y] = 1
                    end
                end
                iter += 1
            end
            return
        else
            divNum  = rand(divRange)
            middleX = round(Int, ( (divNum*rangeX.start) + rangeX.stop) / (divNum+1))
            middleY = round(Int, ( (divNum*rangeY.start) + rangeY.stop) / (divNum+1))

            rangeX1 = rangeX.start:middleX
            rangeX2 = (middleX+1):rangeX.stop

            rangeY1 = rangeY.start:middleY
            rangeY2 = (middleY+1):rangeY.stop

            makeTreeMatrixRandomly!(mat, rangeX1, rangeY1, p)
            makeTreeMatrixRandomly!(mat, rangeX2, rangeY2, p)
        end
    end
end

module phylogenyTreeTest
    using Base.Test

    using phylogenyTree
    using __phylogenyTree
    using __phylogenyTreeTest

    function runAllTest()::Void
        """ radix sort test """
        radixResults      = []
        phyloCheckResults = []
        numTest = 100
        for i in 1:numTest
            binSize = 5
            range   = 1:((2^binSize)-1)
            sortNum = 15
            mat::Array{Int,2} = zeros(Int, binSize, sortNum)

            allNums = []
            for y in 1:sortNum
                randNum = rand(range)
                push!(allNums, randNum)
                for x in 1:binSize
                    mat[x,y] = parse(Int,bin(randNum, binSize)[x])
                end
            end

            order = __phylogenyTree.radixSort(mat, ascendDigit = false, base = 2, descend = false)
            ok = true
            for isSame in (allNums[order].== sort(allNums))
                ok = ok && isSame
            end
            push!(radixResults, ok)
        end
        """ phylogenyTree check Test """
        for i in 1:numTest
            X = 30
            Y = 50
            treeMat = zeros(Int, X, Y);
            __phylogenyTreeTest.makeTreeMatrixRandomly!(treeMat, 1:X, 1:Y);
            push!(phyloCheckResults, phylogenyTree.isPhylogenic( treeMat[ randperm(X), randperm(Y)] ))
        end

        @testset "radix sort"     begin
            for isValid in radixResults;      @test isValid==true; end
        end
        @testset "phyloTreeCheck" begin
            for isValid in phyloCheckResults; @test isValid==true; end
        end
        return nothing
    end
end
