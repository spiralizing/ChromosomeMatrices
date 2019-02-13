using InfoSeries
using LinearAlgebra, DelimitedFiles, StatsBase, Plots, Clustering, Random, Statistics, ArgParse
pyplot(guidefont=20, titlefont=20,xtickfont=10, ytickfont=10)
#fc1 = cgrad([:red,:white,:blue])
fc2 = cgrad([:steelblue,:red,:black])
#fc3 = cgrad([:white,:black])
################################################################################
#Variables
subt = ["Healthy","Basal","Her2", "LumA","LumB"]

#ch = 5
#nmats = 100 #number of surrogate matrices
#argument section
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "-c"
            help = "takes the number of the chromosome as input"
            arg_type = Int
            default = 1
            required = true
    end
    return parse_args(s)
end
################################################################################
#LOADING DATA
ch = parse_commandline()["c"]
ch = 1
c1 = readdlm("ClusterDataTh-$ch-$(subt[2]).dat")
nc1 = maximum(c1[:,2])
#m1 = corr_mat(convert(Array{Float64,2},readdlm("$(subt[2])_chr_$ch.txt")[:, 2:end-1]))


c2 = readdlm("ClusterDataTh-$ch-$(subt[3]).dat")
nc2 = maximum(c2[:,2])
#m2 = corr_mat(convert(Array{Float64,2},readdlm("$(subt[3])_chr_$ch.txt")[:, 2:end-1]))


c3 = readdlm("ClusterDataTh-$ch-$(subt[4]).dat")
nc3 = maximum(c3[:,2])
#m3 = corr_mat(convert(Array{Float64,2},readdlm("$(subt[4])_chr_$ch.txt")[:, 2:end-1]))


c4 = readdlm("ClusterDataTh-$ch-$(subt[5]).dat")
nc4 = maximum(c4[:,2])
#m4 = corr_mat(convert(Array{Float64,2},readdlm("$(subt[5])_chr_$ch.txt")[:, 2:end-1]))

centro = convert(Int64,readdlm("centros.dat")[ch])

#########################################################################################

N = size(c2)[1]
n_a = 0
c = Array{Int64,1}(undef,4)
g_c = Array{Int64,1}

for i = 1:nc1
    for j = 1:nc2
        for k = 1:nc3
            for l = 1:nc4
                ie1 = findall(x -> x == i, c1[:,2])
                ie2 = findall(x -> x == j, c2[:,2])
                ie3 = findall(x -> x == k, c3[:,2])
                ie4 = findall(x -> x == l, c4[:,2])
                inx_int = intersect(ie1, ie2, ie3, ie4)
                ngc = length(inx_int)
                if ngc > n_a
                    g_c = inx_int
                    n_a = ngc
                    c[1] = i; c[2]=j; c[3] = k; c[4] = l
                end
            end
        end
    end
end



#i_e1 = findall(x -> x == c[1], c1[:,2])
#i_e2 = findall(x -> x == c[2], c2[:,2])
#i_e3 = findall(x -> x == c[3], c3[:,2])
#i_e4 = findall(x -> x == c[4], c4[:,2])

#lbl1 = ["$(subt[2])" for i = 1:length(i_e1)]
#st1 = [c1[i_e1,1] lbl1]
#lbl2 = ["$(subt[3])" for i = 1:length(i_e2)]
#st2 = [c2[i_e2,1] lbl2]
#lbl3 = ["$(subt[4])" for i = 1:length(i_e3)]
#st3 = [c2[i_e3,1] lbl3]
#lbl4 = ["$(subt[5])" for i = 1:length(i_e4)]
#st4 = [c2[i_e4,1] lbl4]

#writedlm("GenesInterAll-$ch.dat", vcat(st1,st2,st3,st4))

nmat = zeros(N,N)
for i = 1:(length(g_c)-1)
    for j = i:length(g_c)
        nmat[g_c[i],g_c[j]] = nmat[g_c[j],g_c[i]] = 1
    end
end

for ix = 1:N
    nmat[ix, centro] = nmat[centro,ix] = 2
end
#hm =[]
heatmap(nmat,yflip=true, fillcolor=fc2, title="Traslape Cromosoma $ch", leg=false)

#push!(hm,heatmap(pmat2,yflip=true, fillcolor=fc1,clim=(-1.0,1.0), title="$s2 Ch$ch"))

plot!(size=(1000,1000))
savefig("IntersectAllCenter-$ch.png")
