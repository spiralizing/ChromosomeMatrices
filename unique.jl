using InfoSeries
using LinearAlgebra, DelimitedFiles, StatsBase, Plots, Clustering, Random, Statistics
################################################################################
#GLOBAL
ch = ARGS[1]
#ch = 10
subt = ["Basal","Her2", "LumA","LumB"]
#FILES
l = readdlm("IntersecAll/GenesInterAll-$ch.dat")
s = Array{Any,1}(undef,4)
for i =1:4
    s[i] = l[find(x-> x == subt[i], l[:,2]),1]
end

writedlm("GenesUnique-$ch.dat",intersect(s...))
