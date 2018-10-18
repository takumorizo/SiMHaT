include("include.jl")
Include("phylogenyTree.jl")
Include("sortingCache2D.jl")
Include("phyloMatrix.jl")
Include("buffPhyloMatrix.jl")

module PhylogenyTreeTest
    using Test
    using Random

    using ..PhylogenyTree

    function _makeTreeMatrixRandomly!(mat::Array{Int, 2}, rangeX, rangeY, p = 0.5, divRange=1:3)::Nothing
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

            _makeTreeMatrixRandomly!(mat, rangeX1, rangeY1, p)
            _makeTreeMatrixRandomly!(mat, rangeX2, rangeY2, p)
        end
    end



    function runAllTest()::Nothing
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
                    # mat[x,y] = parse(Int,bin(randNum, binSize)[x])
                    mat[x,y] = parse(Int,string(randNum, base=2, pad=binSize)[x])
                end
            end

            order = PhylogenyTree._radixSort(mat, ascendDigit = false, base = 2, descend = false)
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
            _makeTreeMatrixRandomly!(treeMat, 1:X, 1:Y);
            push!(phyloCheckResults, PhylogenyTree.isPhylogenic( treeMat[ Random.randperm(X), Random.randperm(Y)] ))
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



module SortingCache2DTest
    using Test
    using ..PhylogenyTree
    using ..SortingCache2DType

    function _genRandom2DMatrices(D::I, N::I;
                                 descend::Bool = true,
                                 ascendDigit::Bool = false,
                                 base::I = 2)  where {I <: Integer}
        array = []
        for i in 1:N
            push!(array, rand(1:N))
        end
        matrix = zeros(I, D, N)
        for i in 1:N
            aDigits = digits(array[i], base=2, pad=D)

            aDigits = aDigits[length(aDigits):-1:1] # println(array[i]); println(aDigits)
            matrix[:,i] = aDigits
        end
        sortedAns = PhylogenyTree._radixSort(matrix, descend = true, ascendDigit = false, base = 2) # print("sorted: ");println(array[sortedAns])
        # sortedMat = deepcopy(matrix[:, sortedAns])
        matrixExt = zeros(Int, D, N+2)
        for r in 1:D;for c in 1:N;
            matrixExt[r,c+1] = matrix[r,c]
        end;end
        for i in 1:length(sortedAns)
            sortedAns[i] += 1
        end
        return matrixExt, sortedAns
    end

    function _makeLLLists(matrixExt::AbstractArray{I, 2}, sortedAns::AbstractArray{I, 1}; linkedValue = 1) where {I <: Integer }
        D, N = size(matrixExt)
        LLLeft  = zeros(I, D, N)
        LLRight = zeros(I, D, N)
        N -= 2
        start = 1
        fin = N+2
        for r in 1:D
            from = start
            for c in 1:N
                if matrixExt[r, sortedAns[c]] == 1
                    LLRight[r,from] = sortedAns[c]
                    from = sortedAns[c]
                end
            end
            LLRight[r,from] = fin
        end
        for r in 1:D
            from = fin
            for c in N:-1:1
                if matrixExt[r, sortedAns[c]] == 1
                    LLLeft[r,from] = sortedAns[c]
                    from = sortedAns[c]
                end
            end
            LLLeft[r,from] = start
        end
        return LLLeft, LLRight
    end

    function _makeTreeMatrixRandomly!(mat::Array{I, 2}, rangeX, rangeY, p = 0.5, divRange=1:3)::Nothing where {I <: Integer}
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

            _makeTreeMatrixRandomly!(mat, rangeX1, rangeY1, p)
            _makeTreeMatrixRandomly!(mat, rangeX2, rangeY2, p)
        end
    end

    function _genRandomTree(D::I, N::I) where {I <: Integer}
        matrix = zeros(I, D, N)
        _makeTreeMatrixRandomly!(matrix, 1:D, 1:N)
        sortedAns = PhylogenyTree._radixSort(matrix, descend = true, ascendDigit = false, base = 2)
        matrixExt = zeros(I, D, N+2)
        for r in 1:D;for c in 1:N;
            matrixExt[r,c+1] = matrix[r,c]
        end;end
        for i in 1:length(sortedAns)
            sortedAns[i] += 1
        end
        return matrixExt, sortedAns
    end



    function runAllTest()::Nothing
        # @time @testset "add!/rm! time check random naive"  begin
        #     for i in 1:1
        #         N = 1000
        #         D = 10
        #         matrixExt, sortedAns = _genRandom2DMatrices(D, N)
        #         cache = SortingCache2DType.init(D,N)
        #         for c in 1:N
        #             SortingCache2DType.add!(cache, c, matrixExt[:, c+1])
        #         end
        #         temp::Bool = true
        #         for i in 1:100
        #             for c in 1:N
        #                 rmAt = sortedAns[rand(1:N)] - 1
        #                 SortingCache2DType.rm!(cache, rmAt)
        #                 SortingCache2DType.add!(cache, rmAt, matrixExt[:, rmAt+1])
        #                 temp = xor(temp, PhylogenyTree.isPhylogenic(cache.matrix))
        #             end
        #         end
        #
        #         @test cache.matrix == matrixExt
        #         @test matrixExt[:,cache.sortedCols] == matrixExt[:,sortedAns]
        #
        #         LLLeft, LLRight = _makeLLLists(matrixExt, cache.sortedCols)
        #
        #         @test cache.LLLeft == LLLeft
        #         @test cache.LLRight == LLRight
        #     end
        # end

        @time @testset "add!/rm! time check random matrix"  begin
            for i in 1:1
                N = 1000
                D = 10
                matrixExt, sortedAns = _genRandom2DMatrices(D, N)
                cache = SortingCache2DType.init(D,N)
                for c in 1:N
                    SortingCache2DType.add!(cache, c, matrixExt[:, c+1])
                end

                for i in 1:100
                    for c in 1:N
                        rmAt = sortedAns[rand(1:N)] - 1
                        SortingCache2DType.rm!(cache, rmAt)
                        SortingCache2DType.add!(cache, rmAt, matrixExt[:, rmAt+1])
                    end
                end

                @test cache.matrix == matrixExt
                @test matrixExt[:,cache.sortedCols] == matrixExt[:,sortedAns]

                LLLeft, LLRight = _makeLLLists(matrixExt, cache.sortedCols)

                @test cache.LLLeft == LLLeft
                @test cache.LLRight == LLRight
            end
        end

        # @time @testset "add!/rm! time check cache phylo naive"  begin
        #     for i in 1:1
        #         N = 1000
        #         D = 10
        #         matrixExt, sortedAns = _genRandomTree(D, N)
        #         cache = SortingCache2DType.init(D,N)
        #         for c in 1:N
        #             SortingCache2DType.add!(cache, c, matrixExt[:, c+1])
        #         end
        #         temp::Bool = true
        #         for i in 1:100
        #             for c in 1:N
        #                 rmAt = sortedAns[rand(1:N)] - 1
        #                 SortingCache2DType.rm!(cache, rmAt)
        #                 SortingCache2DType.add!(cache, rmAt, matrixExt[:, rmAt+1])
        #                 temp = xor(temp, PhylogenyTree.isPhylogenic(cache.matrix))
        #             end
        #         end
        #
        #         @test cache.matrix == matrixExt
        #         @test matrixExt[:,cache.sortedCols] == matrixExt[:,sortedAns]
        #
        #         LLLeft, LLRight = _makeLLLists(matrixExt, cache.sortedCols)
        #
        #         @test cache.LLLeft == LLLeft
        #         @test cache.LLRight == LLRight
        #     end
        # end

        @time @testset "add!/rm! time check cache phylo"  begin
            for i in 1:1
                N = 1000
                D = 10
                matrixExt, sortedAns = _genRandomTree(D, N)
                cache = SortingCache2DType.init(D,N)
                for c in 1:N
                    SortingCache2DType.add!(cache, c, matrixExt[:, c+1])
                end

                for i in 1:100
                    for c in 1:N
                        rmAt = sortedAns[rand(1:N)] - 1
                        SortingCache2DType.rm!(cache, rmAt)
                        SortingCache2DType.add!(cache, rmAt, matrixExt[:, rmAt+1])
                    end
                end

                @test cache.matrix == matrixExt
                @test matrixExt[:,cache.sortedCols] == matrixExt[:,sortedAns]

                LLLeft, LLRight = _makeLLLists(matrixExt, cache.sortedCols)

                @test cache.LLLeft == LLLeft
                @test cache.LLRight == LLRight
            end
        end

        @time @testset "add! cache"  begin
            for i in 1:1000
                N = 10
                D = 10
                matrixExt, sortedAns = _genRandom2DMatrices(D, N)
                cache = SortingCache2DType.init(D,N)
                for c in 1:N
                    SortingCache2DType.add!(cache, c, matrixExt[:, c+1])
                end

                @test cache.matrix == matrixExt
                @test matrixExt[:,cache.sortedCols] == matrixExt[:,sortedAns]

                LLLeft, LLRight = _makeLLLists(matrixExt, cache.sortedCols)

                @test cache.LLLeft == LLLeft
                @test cache.LLRight == LLRight
            end
        end

        @time @testset "rm!/add! cache"  begin
            for i in 1:1000
                N = 10
                D = 10
                matrixExt, sort = _genRandom2DMatrices(D, N)
                cache = SortingCache2DType.init(D,N)
                for c in 1:N
                    SortingCache2DType.add!(cache, c, matrixExt[:, c+1])
                end

                matrixAns = deepcopy(cache.matrix)
                sortAns = deepcopy(cache.sortedCols)
                rmAt = sort[rand(1:N)] - 1
                SortingCache2DType.rm!(cache, rmAt)
                SortingCache2DType.add!(cache, rmAt, matrixExt[:, rmAt+1])

                @test cache.matrix == matrixAns == matrixExt
                @test matrixExt[:,cache.sortedCols] == matrixExt[:,sortAns]

                LLLeft, LLRight = _makeLLLists(matrixExt, cache.sortedCols)

                @test cache.LLLeft == LLLeft
                @test cache.LLRight == LLRight
            end
        end

        @time @testset "edit! cache"  begin
            for i in 1:1000
                N = 10
                D = 10
                matrixExt, sort = _genRandom2DMatrices(D, N)
                cache = SortingCache2DType.init(D,N)
                for r in 1:D
                    for c in 1:N
                        for i in 1:1000
                            SortingCache2DType.edit!(cache, r, c, matrixExt[r, c+1])
                        end
                    end
                end
                @test cache.matrix == matrixExt
                @test matrixExt[:,cache.sortedCols] == matrixExt[:,sort]

                LLLeft, LLRight = _makeLLLists(matrixExt, cache.sortedCols)

                @test cache.LLLeft == LLLeft
                @test cache.LLRight == LLRight
            end
        end


        @time @testset "__binSearchAscend! in 2Darray"  begin
            for i in 1:1000
                array = []
                N = 10
                D = 10
                for i in 1:N
                    push!(array, rand(1:N))
                end
                matrix = zeros(Int64,D, N)
                for i in 1:N
                    aDigits = digits(array[i], base=10, pad=D)
                    aDigits = aDigits[length(aDigits):-1:1] # println(array[i]); println(aDigits)
                    matrix[:,i] = aDigits
                end
                sortedAns = PhylogenyTree._radixSort(matrix, ascendDigit = false, base = 10) # print("sorted: ");println(array[sortedAns])
                rmAt = rand(1:N)
                col  = sortedAns[rmAt]
                v    = matrix[:,col]; val  = array[col]

                sorting = deepcopy(sortedAns)
                deleteat!(sorting, rmAt) # print("rmed sorted: ");println(array[sorting]) # print("value: "); println(val)
                insertAt = SortingCache2DType._binSearchAscend!(matrix, sorting, v, col)
                @test array[sorting] == array[sortedAns] # println(val) # println(array[sorting])
            end
        end
        @time @testset "__binSearchDescend! in 2Darray"  begin
            for i in 1:1000
                array = []
                N = 10
                D = 10
                for i in 1:N
                    push!(array, rand(1:N))
                end
                matrix = zeros(Int64,D, N)
                for i in 1:N
                    aDigits = digits(array[i], base=10, pad=D)
                    aDigits = aDigits[length(aDigits):-1:1] # println(array[i]); println(aDigits)
                    matrix[:,i] = aDigits
                end
                sortedAns = PhylogenyTree._radixSort(matrix, descend = true, ascendDigit = false, base = 10) # print("sorted: ");println(array[sortedAns])
                rmAt = rand(1:N)
                col  = sortedAns[rmAt]
                v    = matrix[:,col]; val  = array[col]

                sorting = deepcopy(sortedAns)
                deleteat!(sorting, rmAt) # print("rmed sorted: ");println(array[sorting]) # print("value: "); println(val)
                insertAt = SortingCache2DType._binSearchDescend!(matrix, sorting, v, col)
                @test array[sorting] == array[sortedAns] # println(val) # println(array[sorting])
            end
        end
        @time @testset "__gt in array"  begin
            for i in 1:100
                a = abs(rand(Int)) % 100000
                b = abs(rand(Int)) % 100000
                aDigits = digits(a, base=10, pad=10)
                bDigits = digits(b, base=10, pad=10)

                aDigits = aDigits[length(aDigits):-1:1]
                bDigits = bDigits[length(bDigits):-1:1]

                @test SortingCache2DType._gt(view(aDigits,:), view(bDigits,:)) == ( a > b )
            end
        end
        @time @testset "__lt in array"  begin
            for i in 1:100
                a = abs(rand(Int)) % 100000
                b = abs(rand(Int)) % 100000
                aDigits = digits(a, base=10, pad=10)
                bDigits = digits(b, base=10, pad=10)
                aDigits = aDigits[length(aDigits):-1:1]
                bDigits = bDigits[length(bDigits):-1:1]

                @test SortingCache2DType._lt(view(aDigits,:), view(bDigits,:)) == ( a < b )
            end
        end
        @time @testset "__eq in array"  begin
            for i in 1:100
                a = abs(rand(Int)) % 100000
                b = abs(rand(Int)) % 100000
                aDigits = digits(a, base=10, pad=10)
                bDigits = digits(b, base=10, pad=10)
                aDigits = aDigits[length(aDigits):-1:1]
                bDigits = bDigits[length(bDigits):-1:1]

                @test SortingCache2DType._eq(view(aDigits,:), view(bDigits,:)) == ( a == b )
            end
        end

        return nothing
   end
end


module PhyloMatrixTest
    using Test
    using ..SortingCache2DTest
    using ..PhyloMatrixType
    using ..PhylogenyTree
    using ..SortingCache2DType

    function _makeTreeMatrixRandomly!(mat::Array{I, 2}, rangeX, rangeY, p = 0.5, divRange=1:3)::Nothing where {I <: Integer}
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

            _makeTreeMatrixRandomly!(mat, rangeX1, rangeY1, p)
            _makeTreeMatrixRandomly!(mat, rangeX2, rangeY2, p)
        end
    end

    function _genRandomTree(D::I, N::I) where {I <: Integer}
        matrix = zeros(I, D, N)
        _makeTreeMatrixRandomly!(matrix, 1:D, 1:N)
        sortedAns = PhylogenyTree._radixSort(matrix, descend = true, ascendDigit = false, base = 2)
        matrixExt = zeros(I, D, N+2)
        for r in 1:D;for c in 1:N;
            matrixExt[r,c+1] = matrix[r,c]
        end;end
        for i in 1:length(sortedAns)
            sortedAns[i] += 1
        end
        return matrixExt, sortedAns
    end

    function runAllTest()::Nothing
        @time @testset "add!/rm! time check cache"  begin
            for i in 1:1
                N = 1000
                D = 10
                matrixExt, sortedAns = _genRandomTree(D, N)
                cache = SortingCache2DType.init(D,N)
                for c in 1:N
                    SortingCache2DType.add!(cache, c, matrixExt[:, c+1])
                end

                for i in 1:100
                    for c in 1:N
                        rmAt = sortedAns[rand(1:N)] - 1
                        SortingCache2DType.rm!(cache, rmAt)
                        SortingCache2DType.add!(cache, rmAt, matrixExt[:, rmAt+1])
                    end
                end

                @test cache.matrix == matrixExt
                @test matrixExt[:,cache.sortedCols] == matrixExt[:,sortedAns]

                LLLeft, LLRight = SortingCache2DTest._makeLLLists(matrixExt, cache.sortedCols)

                @test cache.LLLeft == LLLeft
                @test cache.LLRight == LLRight
            end
        end

        @time @testset "add! phyloMat random"  begin
            for i in 1:1000
                N = 10
                D = 5
                if (i % 2 == 0)
                    matrixExt, sortedAns =  SortingCache2DTest._genRandom2DMatrices(D, N)
                else
                    matrixExt, sortedAns =  _genRandomTree(D, N)
                end

                phylo = PhyloMatrixType.init(D,N)
                for c in 1:N
                    PhyloMatrixType.add!(phylo, c, matrixExt[:, c+1])
                end

                @test phylo.cache.matrix == matrixExt
                @test matrixExt[:,phylo.cache.sortedCols] == matrixExt[:,sortedAns]

                LLLeft, LLRight = SortingCache2DTest._makeLLLists(matrixExt, phylo.cache.sortedCols)

                @test phylo.cache.LLLeft == LLLeft
                @test phylo.cache.LLRight == LLRight

                isTree::Bool = PhyloMatrixType.isTree(phylo)
                if isTree != PhylogenyTree.isPhylogenic(phylo.cache.matrix[:, 2:(N+1)])
                    println("LLLeft")
                    println(LLLeft)

                    println("matrix")
                    println(phylo.cache.sortedCols)
                    println(matrixExt)
                    println(matrixExt[:, phylo.cache.sortedCols])
                    println(phylo.Ly)
                    for c in 1:N
                        println(matrixExt[:, c+1])
                    end
                end
                @test isTree == PhylogenyTree.isPhylogenic(phylo.cache.matrix[:, 2:(N+1)])
            end
        end

        @time @testset "rm! phyloMat random"  begin
            for i in 1:1000
                N = 10
                D = 5

                if (i % 2 == 0)
                    matrixExt, sortedAns =  SortingCache2DTest._genRandom2DMatrices(D, N)
                else
                    matrixExt, sortedAns =  _genRandomTree(D, N)
                end

                phylo = PhyloMatrixType.init(D,N)
                for c in 1:N
                    PhyloMatrixType.add!(phylo, c, matrixExt[:, c+1])
                end

                PhyloMatrixType.rm!(phylo, N)

                isTree::Bool = PhyloMatrixType.isTree(phylo)
                @test isTree == PhylogenyTree.isPhylogenic(matrixExt[:, 2:N])
            end
        end

        @time @testset "edit! phyloMat random"  begin
            for i in 1:1000
                N = 20
                D = 10
                if (i % 2 == 0)
                    matrixExt, sortedAns =  SortingCache2DTest._genRandom2DMatrices(D, N)
                else
                    matrixExt, sortedAns =  _genRandomTree(D, N)
                end

                phylo = PhyloMatrixType.init(D,N)
                for r in 1:D
                    for c in 1:N
                        PhyloMatrixType.edit!(phylo, r, c, matrixExt[r, c+1])
                    end
                end

                for r in 1:D
                    for c in 1:N
                        PhyloMatrixType.edit!(phylo, r, c, rand(0:1))
                    end
                end

                for r in 1:D
                    for c in 1:N
                        PhyloMatrixType.edit!(phylo, r, c, matrixExt[r, c+1])
                    end
                end

                @test phylo.cache.matrix == matrixExt
                @test matrixExt[:,phylo.cache.sortedCols] == matrixExt[:,sortedAns]

                LLLeft, LLRight = SortingCache2DTest._makeLLLists(matrixExt, phylo.cache.sortedCols)

                @test phylo.cache.LLLeft == LLLeft
                @test phylo.cache.LLRight == LLRight

                isTree::Bool = PhyloMatrixType.isTree(phylo)
                @test isTree == PhylogenyTree.isPhylogenic(phylo.cache.matrix[:, 2:(N+1)])
            end
        end
        return nothing
    end
end


module BuffPhyloMatrixTest
    using Test
    using ..SortingCache2DTest
    using ..PhyloMatrixType
    using ..PhyloMatrixTest
    using ..PhylogenyTree
    using ..BuffPhyloMatrixType
    function _convertToBiCluster(matrix::Array{I, 2},
                                ; width::I = 1, height::I = 1) where {I <: Integer}
        X, Y = size(matrix)
        usageS::Dict{I,Array{I,1}} = Dict{I,Array{I,1}}()
        usageV::Dict{I,Array{I,1}} = Dict{I,Array{I,1}}()
        B::Dict{Tuple{I,I}, I} = Dict{Tuple{I,I}, I}()

        rangeX = (width * 1):width:(width * X)
        rangeY = (height * 1):height:(height * Y)

        for x in rangeX
            v::Array{I,1} = []
            for w in (x-width+1):x
                push!(v, w)
            end
            usageS[x] = v
        end

        for y in rangeY
            v::Array{I,1} = []
            for h in (y-height+1):y
                push!(v, h)
            end
            usageV[y] = v
        end

        for c in rangeX
            for m in rangeY
                x = div(c, width)
                y = div(m, height)
                B[(c,m)] = matrix[x,y]
            end
        end

        return B, usageS, usageV
    end

    function runAllTest()::Nothing
        @time @testset "add! BuffPhyloMat random"  begin
            for i in 1:1000
                N = 200
                D = 10
                if (i % 2 == 0)
                    matrixExt, sortedAns =  SortingCache2DTest._genRandom2DMatrices(D, N)
                else
                    matrixExt, sortedAns =  PhyloMatrixTest._genRandomTree(D, N)
                end

                bPhylo = BuffPhyloMatrixType.init(D,N)
                for c in 1:N
                    BuffPhyloMatrixType.add!(bPhylo, c, matrixExt[:, c+1])
                end

                @test matrixExt[:,bPhylo.phylo.cache.sortedCols] == matrixExt[:,sortedAns]
                @test bPhylo.phylo.cache.matrix[:,bPhylo.phylo.cache.sortedCols] == bPhylo.phylo.cache.matrix[:,sortedAns]
                isTree::Bool = BuffPhyloMatrixType.isTree(bPhylo)
                @test isTree == PhylogenyTree.isPhylogenic(bPhylo.phylo.cache.matrix[:, 2:(N+1)])
            end
        end
        @time @testset "rm! BuffPhyloMat random"  begin
            for i in 1:1000
                N = 200
                D = 20

                if (i % 2 == 0)
                    matrixExt, sortedAns =  SortingCache2DTest._genRandom2DMatrices(D, N)
                else
                    matrixExt, sortedAns =  PhyloMatrixTest._genRandomTree(D, N)
                end

                bPhylo = BuffPhyloMatrixType.init(D,N)
                for c in 1:N
                    BuffPhyloMatrixType.add!(bPhylo, c, matrixExt[:, c+1])
                end

                BuffPhyloMatrixType.rm!(bPhylo, N)
                isTree::Bool = BuffPhyloMatrixType.isTree(bPhylo)
                @test isTree == PhylogenyTree.isPhylogenic(bPhylo.phylo.cache.matrix[:, 2:N])
            end
        end
        @time @testset "edit! BuffPhyloMat random"  begin
            for i in 1:1000
                N = 20
                D = 10
                if (i % 2 == 0)
                    matrixExt, sortedAns =  SortingCache2DTest._genRandom2DMatrices(D, N)
                else
                    matrixExt, sortedAns =  PhyloMatrixTest._genRandomTree(D, N)
                end

                bPhylo = BuffPhyloMatrixType.init(D,N)
                for r in 1:D
                    for c in 1:N
                        BuffPhyloMatrixType.edit!(bPhylo, r, c, matrixExt[r, c+1])
                    end
                end

                isTree::Bool = BuffPhyloMatrixType.isTree(bPhylo)
                @test isTree == PhylogenyTree.isPhylogenic(bPhylo.phylo.cache.matrix[:, 2:N])
            end
        end

        @time @testset "update! BuffPhyloMat"  begin
            for i in 1:10000
                N1, N2 = (5, 5)
                D = 5
                width = 1
                height = 1
                matrixExt, sortedAns = PhyloMatrixTest._genRandomTree(D, N1+N2)
                if (i % 2 != 0)
                    matrixExtR, sortedAnsR =   SortingCache2DTest._genRandom2DMatrices(D, N1+N2)
                    matrixExt[:, (1+N1+1):(1+N1+N2)] = deepcopy(matrixExtR[:, (1+1):(N2+1)])
                end

                BL, usageSL, usageVL = _convertToBiCluster(matrixExt[:,(1+1):(1+N1)], width = width, height = height)
                B,  usageS,  usageV  = _convertToBiCluster(matrixExt[:,(1+1):(1+N1+N2)], width = width, height = height)

                bPhylo = BuffPhyloMatrixType.init(D * width, 1, bufferSize = (N1+N2) * height)

                BuffPhyloMatrixType.update!(bPhylo, B, usageS, usageV)
                isTreeAns::Bool = BuffPhyloMatrixType.isTree(bPhylo)
                @test isTreeAns == PhylogenyTree.isPhylogenic( matrixExt[:, 2:(N1+N2+1)] )

                BuffPhyloMatrixType.update!(bPhylo, BL, usageSL, usageVL)

                @test PhylogenyTree.isPhylogenic(matrixExt[:,(1+1):(1+N1)])
                @test BuffPhyloMatrixType.isTree(bPhylo) == true

                BuffPhyloMatrixType.update!(bPhylo, B, usageS, usageV)
                isTree::Bool = BuffPhyloMatrixType.isTree(bPhylo)

                @test isTreeAns == PhylogenyTree.isPhylogenic(matrixExt[:, 2:(N1+N2+1)])
                @test isTreeAns == isTree
            end
        end

        @time @testset "update! BuffPhyloMat random"  begin
            for i in 1:1000
                N = 10
                D = 5
                width = 1
                height = 1
                matrixExt1, sortedAns1 = PhyloMatrixTest._genRandomTree(D, N)
                if (i % 2 != 0)
                    matrixExt2, sortedAns2 =   SortingCache2DTest._genRandom2DMatrices(D, N)
                else
                    matrixExt2, sortedAns2 =  PhyloMatrixTest._genRandomTree(D, N)
                end

                B1, usageS1, usageV1 = _convertToBiCluster(matrixExt1[:,(1+1):(1+N)], width = width, height = height)
                B2, usageS2, usageV2 = _convertToBiCluster(matrixExt2[:,(1+1):(1+N)], width = width, height = height)

                bPhylo = BuffPhyloMatrixType.init(D * width, 1, bufferSize = (N) * height)

                BuffPhyloMatrixType.update!(bPhylo, B1, usageS1, usageV1)
                isTree1::Bool = BuffPhyloMatrixType.isTree(bPhylo)
                @test isTree1 == PhylogenyTree.isPhylogenic( matrixExt1[:, 2:(N+1)] )

                BuffPhyloMatrixType.update!(bPhylo, B2, usageS2, usageV2)
                isTree2::Bool = BuffPhyloMatrixType.isTree(bPhylo)
                @test isTree2 == PhylogenyTree.isPhylogenic( matrixExt2[:, 2:(N+1)] )
            end
        end

        return nothing
    end
end
