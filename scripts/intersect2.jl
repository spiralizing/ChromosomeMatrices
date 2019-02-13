using InfoSeries
using LinearAlgebra, DelimitedFiles, StatsBase, Plots, Clustering, Random, Statistics
pyplot(guidefont=20, titlefont=20,xtickfont=10, ytickfont=10)
fc1 = cgrad([:red,:white,:blue])
#fc2 = cgrad([:red,:white])
#fc3 = cgrad([:white,:black])
################################################################################
#Variables
subt = ["Healthy","Basal","Her2", "LumA","LumB"]
ch = ARGS[1]
#ch = 21
#nmats = 100 #number of surrogate matrices
################################################################################
for  ix1 = 1:4
    for  ix2 = ix1:5
        s1 = subt[ix1]
        c1 = readdlm("ClusterData-Ch$ch-$s1.dat")
        nc1 = maximum(c1[:,2])
        m1 = corr_mat(convert(Array{Float64,2},readdlm("$(s1)_chr_$ch.txt")[:, 2:end-1]))

        s2 = subt[ix2]
        c2 = readdlm("ClusterData-Ch$ch-$s2.dat")
        nc2 = maximum(c2[:,2])
        m2 = corr_mat(convert(Array{Float64,2},readdlm("$(s2)_chr_$ch.txt")[:, 2:end-1]))

        N = size(c2)[1]
        n_a = 0
        c = Array{Int64,1}(undef,2)

        for i = 1:nc1
            for j = 1:nc2
                ie1 = findall(x -> x == i, c1[:,2])
                ie2 = findall(x -> x == j, c2[:,2])
                ngc = length(intersect(ie1, ie2))
                tst = ngc/length(ie1) + ngc/length(ie2)
                if tst > n_a
                    n_a = tst
                    c[1] = i; c[2]=j
                end
            end
        end

        i_e1 = findall(x -> x == c[1], c1[:,2])
        i_e2 = findall(x -> x == c[2], c2[:,2])
        #pmat1 = zeros(N,N)
        ng1 = length(i_e1)
        #for i = 1:ng1-1
        #    for j = i+1:ng1
        #        pmat1[i_e1[i],i_e1[j]] = pmat1[i_e1[j],i_e1[i]] = m1[i_e1[i], i_e1[j]]
        #    end
        #end
        #pmat2 = zeros(N,N)
        ng2 = length(i_e2)
        #for i = 1:ng2-1
        #    for j = i+1:ng2
        #        pmat2[i_e2[i],i_e2[j]] = pmat2[i_e2[j],i_e2[i]] = m2[i_e2[i], i_e2[j]]
        #    end
        #end
        lbl1 = ["$s1" for i = 1:ng1]
        st1 = [c1[i_e1,1] lbl1]

        lbl2 = ["$s2" for i = 1:ng2]
        st2 = [c2[i_e2,1] lbl2]

        writedlm("GenesInterPerc-$s1-$s2-$ch.dat", vcat(st1,st2))

        #hm =[]

        #push!(hm,heatmap(pmat1,yflip=true, fillcolor=fc1,clim=(-1.0,1.0), title="$s1 Ch$ch"))
        #push!(hm,heatmap(pmat2,yflip=true, fillcolor=fc1,clim=(-1.0,1.0), title="$s2 Ch$ch"))

        #plot(hm..., size=(1200,600))
        #savefig("IntersectPerc-$ch-$s1-$s2.png")
    end
end
