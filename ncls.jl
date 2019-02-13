#IMPORTS
using InfoSeries
using LinearAlgebra, DelimitedFiles, StatsBase, Plots, Clustering, Random, Statistics
pyplot(guidefont=20, titlefont=20,xtickfont=10, ytickfont=10)
################################################################
chrs = [i for i=1:22]
subt = ["Healthy","Basal","Her2", "LumA","LumB"]
##################################################
ncls = zeros(22,5)
for i = 1:22
    j = 1
    for s in subt
        ass = readdlm("Clusters/assignments-Ch$i-$s.dat")
        ncls[i,j] = maximum(ass[:,2])
        j+=1
    end
end

#ncls
plot(m=:o, ncls[:,1], label="$(subt[1])")
for i = 2:5
    plot!(m=:o, ncls[:,i], label="$(subt[i])")
end

plot!(size=(1200,800))
savefig("nclusters-todos.png")
