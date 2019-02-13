using InfoSeries
using LinearAlgebra, DelimitedFiles, StatsBase, Plots, Clustering, Random, Statistics

subt = ["Healthy","Basal","Her2", "LumA","LumB"]
ch = ARGS[1]

for s in subt
    expr = readdlm("$(s)_chr_$ch.txt")

    writedlm("ListAllGenes-Ch$(ch)-$s.txt", [expr[:,1] expr[:,end]])
end
