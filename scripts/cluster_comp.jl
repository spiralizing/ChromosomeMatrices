
# Next files are the imports, InfoSeries is my personal repository.
using InfoSeries
using LinearAlgebra, DelimitedFiles, StatsBase, Plots, Clustering, Random, Statistics
pyplot(guidefont=20, titlefont=20,xtickfont=10, ytickfont=10)
fc1 = cgrad([:red,:white,:blue])
#fc2 = cgrad([:red,:white])
#fc3 = cgrad([:white,:black])
################################################################################
#Variables
subt = ["Healthy","Basal","Her2", "LumA","LumB"]
s = ARGS[1]
chs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"]
#nmats = 100 #number of surrogate matrices
################################################################################
#for s in subt
function number_c(gene) #returns the number of clusters.
    nmats = 1000

    surr_ev = []

    for t = 1:nmats
        gran = zeros(size(gene))
        for i = 1:size(gene)[2]
            gran[:,i] = gene[shuffle(1:end),i]
        end
        mr = corr_mat(gran)
        push!(surr_ev,eigvals(mr))
    end
    surr_ev = filter(x->x>0.000001,vcat(surr_ev...))
    m = corr_mat(gene)
    m_ev = filter(x->x>0.0000001, eigvals(m))
    #histogram(surr_ev, normalized=true, alpha=0.5)
    #histogram!(m_ev,normalized=true,alpha=0.5)
    #plot!(size=(1200,900))

    nc = length(filter(x->x>maximum(surr_ev),m_ev))
    #println(nc)
    dmat = map(x-> 1-abs(x),m)

    return nc
end


################################################################################

#s = subt[sn]
expr = []
for ch in chs
    push!(expr, readdlm("$(s)_chr_$ch.txt"))
end

expr = vcat(expr...)
gene = convert(Array{Float64,2},expr[:, 2:end-1])


ens_c = []
nc = number_c(gene)
m = corr_mat(gene)
dmat = map(x-> 1-abs(x),m)


for i = 1:500
    push!(ens_c, kmedoids(dmat,nc))
end

inx_c = findmin(map(x -> x.totalcost, ens_c))[2]
c = ens_c[inx_c]
writedlm("ClusterData-ALL-$s.dat", [expr[:,1] c.assignments c.acosts])
