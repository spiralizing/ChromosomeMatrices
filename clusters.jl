using InfoSeries
using LinearAlgebra, DelimitedFiles, StatsBase, Plots, Clustering, Random, Statistics
pyplot(guidefont=20, titlefont=20,xtickfont=10, ytickfont=10)
fc1 = cgrad([:red,:white,:blue])
#fc2 = cgrad([:red,:white])
#fc3 = cgrad([:white,:black])
################################################################################
#Variables
subt = ["Healthy","Basal","Her2", "LumA","LumB"]
#ch = ARGS[1]
ch = 19
nmats = 100 #number of surrogate matrices
################################################################################
#for s in subt
s = subt[4]
expr = readdlm("$(s)_chr_$ch.txt")
gene = convert(Array{Float64,2},expr[:, 2:end-1])
surr_ev = []
cvals = []
N = size(gene)[1]

for t = 1:nmats
    gran = zeros(size(gene))
    for i = 1:size(gene)[2]
        gran[:,i] = gene[shuffle(1:end),i]
    end
    mr = corr_mat(gran)
    for i = 1:(N-1)
        for j=(i+1):N
            push!(cvals,mr[i,j])
        end
    end
    push!(surr_ev,eigvals(mr))
end

cvals = convert(Array{Float64,1}, cvals)
qs = quantile(cvals,[0.25,0.75])

iq_r = qs[2]-qs[1] #interquartile range
tb = qs[1] - 1.5*iq_r
tt = qs[2] + 1.5*iq_r
th = 1 - max(tt, abs(tb))

surr_ev = filter(x->x>0.000001,vcat(surr_ev...))
m = corr_mat(gene)
m_ev = filter(x->x>0.0000001, eigvals(m))
#histogram(surr_ev, normalized=true, alpha=0.5)
#histogram!(m_ev,normalized=true,alpha=0.5)
#plot!(size=(1200,900))

nc = length(filter(x->x>maximum(surr_ev),m_ev))
dmat = map(x-> 1-abs(x),m)

groups = kmedoids(dmat,nc)


i_eg = findall(x -> x == 5, groups.assignments)
pmat = zeros(N,N)
ng = length(i_eg)
for i = 1:ng-1
    for j = i+1:ng
        pmat[i_eg[i],i_eg[j]] = pmat[i_eg[j],i_eg[i]] = m[i_eg[i], i_eg[j]]
    end
end

heatmap(pmat,yflip=true, fillcolor=fc1, clim=(-1.0,1.0), title="$s Ch$ch")
plot!(size=(900,900))


m_p = m[i_eg,i_eg]
mean(m_p)
scatter(sort(groups.costs[i_eg]),label="cluster 1")



len = zeros(nc)
for i = 1:nc
    i_eg = findall(x -> x == i, groups.assignments)
    len[i] = length(i_eg)
    #d_c[i] = mean(groups.acosts[i_eg])
    #cm = filter(x-> x < 0.5, groups.acosts[i_eg])
    #scatter!(sort(groups.acosts[i_eg]), label="cluster $i")
end
len









gmat = zeros(size(m))
tam = size(m)[1]
for i = 1:(tam-1)
    for j = i:tam
        if groups.assignments[i] == groups.assignments[j]
            gmat[i,j] = gmat[j,i] = groups.assignments[i]
        end
    end
end
writedlm("assignments-Ch$(ch)-$s.dat", [expr[:,1] groups.assignments])

heatmap(gmat,yflip=true, fillcolor=:lighttest, title="$s Ch$ch")
plot!(size=(1200,1200))
savefig("clusters-Ch$(ch)-$s.png")
#end
################################################################################
#heatmap(mr, fillcolor=fc1, clim=(-1.0,1.0),title="Basal", yflip=true)
