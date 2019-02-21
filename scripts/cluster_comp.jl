# Next files are the imports
using Distributed
#addprocs(4)
@everywhere include("/home/alfredo/Git/ChromosomeMatrices/scripts/funcs.jl")
@everywhere using SharedArrays, LinearAlgebra, DelimitedFiles, StatsBase, Plots, Clustering, Random, Statistics
#pyplot(guidefont=20, titlefont=20,xtickfont=10, ytickfont=10)
#fc1 = cgrad([:red,:white,:blue])
#fc2 = cgrad([:red,:white])
#fc3 = cgrad([:white,:black])
################################################################################
#parsing
#function parse_commandline()
#    s = ArgParseSettings()
#    @add_arg_table s begin
#        "-c"
#            help = "takes the number of the chromosome as input"
#            arg_type = Int
#            default = 1
#            required = true
#    end
#    return parse_args(s)
#end
#Variables
#subt = ["Healthy","Basal","Her2", "LumA","LumB"]
#s = ARGS[1]
#chs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"]
################################################################################

nmats = 1000

@everywhere function max_eigv(gene) #returns largest eigenvalue of a randomized matrix generated from data.
    gran = zeros(size(gene))
    for i = 1:size(gene,2)
        gran[:,i] = gene[shuffle(1:end),i]
    end
    mr = corr_mat(gran)
    return maximum(eigvals(mr))
end



################################################################################

#s = subt[5]
#expr = []
#for ch in chs
#    push!(expr, readdlm("/home/alfredo/Git/ChromosomeMatrices/data/$(s)_chr_$ch.txt"))
#end

#expr = vcat(expr...)
#writedlm("$s-All.txt",expr)

@everywhere expr = readdlm("/home/alfredo/Git/ChromosomeMatrices/Basal-All.txt")

@everywhere gene = convert(Array{Float64,2},expr[:, 2:end-1])


sh_evmax = SharedArray{Float64, 1}((nmats))
@sync @distributed for i in 1:nmats #doing the surrogate matrices and saving max eigv.
    sh_evmax[i] = max_eigv(gene)
end

l_evmax = fetch(sh_evmax)
m = corr_mat(gene)
m_ev = filter(x->x>0.0000001, eigvals(m))
nc = length(filter(x->x>maximum(l_evmax),m_ev))
dmat = map(x-> 1-abs(x),m)

ens_c = []
for i = 1:500
    push!(ens_c, kmedoids(dmat,nc))
end


inx_c = findmin(map(x -> x.totalcost, ens_c))[2]
c = ens_c[inx_c]
writedlm("ClusterData-ALL-$s.dat", [expr[:,1] c.assignments c.acosts])
