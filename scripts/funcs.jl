#this are some functions for the scripts
################################################################################
function corr_mat(gene)
    ng = size(gene)[1]
    c_mat = zeros(ng,ng)
    for i = 1:(ng-1)
        c_mat[i,i] = 1.0
        for j = (i+1):ng
            c_mat[j,i] = c_mat[i,j] = cor(gene[i,:],gene[j,:])
            #sp_mat[i,j] = corspearman(gene[i,:],gene[j,:])
        end
    end
    return c_mat
end
################################################################################
