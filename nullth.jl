using InfoSeries
using LinearAlgebra, DelimitedFiles, StatsBase, Plots, Random, Statistics, Clustering
pyplot(guidefont=20, titlefont=20,xtickfont=10, ytickfont=10)
fc1 = cgrad([:red,:white,:blue])
################################################################################
#Variables
subt = ["Healthy","Basal","Her2", "LumA","LumB"]
#ch = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"]
ch = ARGS[1]
ncalc = 99
#nmats = 100 #number of surrogate matrices
################################################################################
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

for s in subt
#s = subt[2]
    #expr = []

    expr = readdlm("$(s)_chr_$(ch).txt")

    gene = convert(Array{Float64,2},expr[:, 2:end-1])

    N = size(expr)[1]
    vals =[]
    for nc = 1:ncalc
        mran = zeros(size(gene))
        for k = 1:size(gene)[2]
            mran[:,k] = gene[shuffle(1:end),k]
        end
        for i = 1:(N-1)
            for j=(i+1):N
                push!(vals,cor(mran[i,:], mran[j,:]))
            end
        end
        #push!(vals, [ix c_d])
    end
            #now estimating quartiles 1 and 3
    vals = convert(Array{Float64,1}, vals)
    qs = quantile(vals,[0.25,0.75])
    iq_r = qs[2]-qs[1] #interquartile range
    tb = qs[1] - 1.5*iq_r; tt = qs[2] + 1.5*iq_r

    m_or = corr_mat(gene)


    m_tr = zeros(N,N)

    for i = 1:(N-1)
        for j = i:N
            if m_or[i,j] > tt
                m_tr[i,j] = m_tr[j,i] = m_or[i,j]
            elseif m_or[i,j] < tb
                m_tr[i,j] = m_tr[j,i] = m_or[i,j]
            end
        end
    end



    ens_c = []
    nc = number_c(gene)
    #m = corr_mat(gene)
    dmat = map(x-> 1-abs(x),m_or)


    for i = 1:500
        push!(ens_c, kmedoids(dmat,nc))
    end

    inx_c = findmin(map(x -> x.totalcost, ens_c))[2]
    c = ens_c[inx_c]
    writedlm("ClusterDataTh-$(ch)-$s.dat", [expr[:,1] c.assignments c.acosts])


    heatmap(m_tr,yflip=true, fillcolor=fc1,clim=(-1,1), title="Threshold $s Ch $(ch)")
    plot!(size=(1200,1200))

    savefig("PearsonThreshold-$ch-$s.png")
#writedlm("HMPearsonThresholdALL-$s.dat",m_tr)
#writedlm("HMPearsonMaxMin-$s.dat", m_mm)

end
