using InfoSeries, Plots, DelimitedFiles, Distributions, StatsBase

#file = ARGS[1]
#expr = readdlm(file)
labels = readcsv("mart_export.txt")
hl = labels[1,:]



for j = 1:22
    inx = find(x-> x == j, labels[:, 4]) #find the genes that correspond to chromosome j
    ch_mat = labels[inx,:]
    ch_mat = ch_mat[sortperm(ch_mat[:,3]), :] #order the chromosomes by place (order)
    writedlm("Ch$j-data.dat",ch_mat)
end
#vcat(Nmat_ch...)

for j = 1:22
    d1 = readdlm("Ch$j-data.dat")
    ix = findfirst(map(x->in('q', x), d1[:,5]))
    nb = d1[ix,3]

    genes = readdlm("Basal_chr_$j.txt")
    not_int = find(map(x->isa(x,SubString),genes[:,end]))
    for x in not_int
        genes[x,end] = 0
    end
    loc = convert(Array{Int64,1},genes[:,end])

    centros[j] = findfirst(x-> x>= nb, loc)
end
