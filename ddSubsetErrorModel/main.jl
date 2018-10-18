include("include.jl")
Include("result.jl")
Include("config.jl")
Include("sampler.jl")

using ..SamplerType
using ..Result
using DocOpt
using Plots
using Dates
using ..config
pyplot()


nowTime = (Dates.value(Dates.now())) #Int(Dates.now())
cwd = pwd()
doc = """ddSubsetErrorModel

Usage:
    main.jl MAP <errScore> <matScore> <patScore> <iniFile> [options]
    main.jl EVAL <score> <answer> <result> [options]

Options:
  -o <DIR> --outDir=<DIR>  Output directory. We use current working directory if unspecified. [default: $cwd]
  -s <SEED> --seed=<SEED>  Seed of randomness. We use current time (msec) if unspecified. [default: $nowTime]
  -n <NUM> --number=<NUM>  Number of iteration after burnin. [default: 1000]
  -b <BURNIN> --burnin=<BURNIN>  Number of iteration during burnin. [default: 0]
  -t <THIN> --thin=<THIN>  Duration between sampling. [default: 1]
  -m <THRES> --threshold=<THRES>  Threshold value for call mutation from file. [default: 0.0]
  --rmError  Remove error if sample answer is used.
  --tag=<TAG>  Tag information for storing Fvalues. [default: ]
  -h --help  Show this screen.

"""
# println(doc)
args = docopt(doc)
println(args)

if args["MAP"] == "MAP"
    @time smax, lnprobs  = SamplerType.exec_map(args["<errScore>"], args["<matScore>"], args["<patScore>"], args["<iniFile>"],
                                                seed = parse(INT, args["--seed"]), iter = parse(INT, args["--number"]), thin = parse(INT, args["--thin"]), burnin = parse(INT, args["--burnin"]))
    Result.writejld(smax, lnprobs, args["--outDir"] * "/sampleAns.jld2")
    Result.viewmap(smax, lnprobs, args["--outDir"])
elseif args["EVAL"] == "EVAL"
    precision, recall, fvalue = Result.evaluatefvalue(args["<answer>"], args["<score>"], parse(REAL, args["--threshold"]), args["--rmError"])
    println((precision, recall, fvalue))
    Result.savefvalue(args["<result>"], precision, recall, fvalue, String(args["--tag"]))
end
