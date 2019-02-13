using InfoSeries
using LinearAlgebra, DelimitedFiles, StatsBase, Plots, Random, Statistics
pyplot(guidefont=20, titlefont=20,xtickfont=10, ytickfont=10)
################################################################################
#Variables
subt = ["Healthy","Basal","Her2", "LumA","LumB"]
ch = ARGS[1]
#ch = 1
#ncalc = 800
#nmats = 100 #number of surrogate matrices
################################################################################

for s in subt
#s = subt[4]
    expr = readdlm("$(s)_chr_$ch.txt")
    gene = convert(Array{Float64,2},expr[:, 2:end-1])
    #mats_cor = []
    N = size(expr)[1]
    vals =[]
    for d = 1:(N-1)
        c_d = []
        for i = 1:(N-d)
            push!(c_d,cor(gene[i,:], gene[i+d,:]))
        end
        ix = [d for i = 1:length(c_d)]
        push!(vals, [ix c_d])
    end
    pcor_or = vcat(vals...)
#mran = zeros(size(gene))
#for i = 1:size(gene)[2]
#    mran[:,i] = gene[shuffle(1:end),i]
#end


#vals_r =[]
#for d = 418:(N-1)
#    c_d = []
#    for j = 1:ncalc
#        mran = zeros(size(gene))
#        for k = 1:size(gene)[2]
#            mran[:,k] = gene[shuffle(1:end),k]
#        end
#        coin = rand(1:2)
#        if coin==1
#            g1 = rand(1:(N-d-1))
#            g2 = g1+d
#        else
#            g1 = rand(d+1:N)
#            g2 = g1-d
#        end
#        push!(c_d,cor(mran[g1,:], mran[g2,:]))
#    end
#    ix = [d for i = 1:length(c_d)]
#    push!(vals_r, [ix c_d])
#end
#pcor_s = vcat(vals_r...)

    scatter(pcor_or[:,1], pcor_or[:,2], ms=3, label="Original", title="Chromosome $ch $s")
#scatter!(pcor_s[:,1], pcor_s[:,2], ms=3, label="Shuffle")
    plot!(size=(1200,800))
    savefig("PearsonDistance-Ch$(ch)-$s.png")
end
################################################################################
#heatmap(mr, fillcolor=fc1, clim=(-1.0,1.0),title="Basal", yflip=true)
