using InfoSeries
using LinearAlgebra, DelimitedFiles, StatsBase, Plots, Clustering, Random, Statistics
pyplot(guidefont=20, titlefont=20,xtickfont=10, ytickfont=10)
################################################################################
#Variables
subt = ["Healthy","Basal","Her2", "LumA","LumB"]
ch = ARGS[1]
ch = 8
#number of surrogate matrices
################################################################################
s = subt[4]
#for s in subt
cdata = readdlm("ClusterDataTh-$ch-$s.dat")
tam = size(cdata)[1]
inx = convert(Array{Int64,1},cdata[:,2])

gmat = zeros(tam,tam)

for i = 1:(tam-1)
    for j = i:tam
        if cdata[i,2] == cdata[j,2]
            gmat[i,j] = gmat[j,i] = inx[i]
        end
    end
end


heatmap(gmat,yflip=true, fillcolor=:lighttest, title="Chromosome $ch $s", xlabel="Gene", ylabel="Gene")
gui()

plot!(size=(1200,1200))

#savefig("clustersTh-Ch$(ch)-$s.png")
#end
