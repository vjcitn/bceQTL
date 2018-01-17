
runmodel = function (eset, snpset, ntokeep = 50) 
{
#
# let K be the length of snpset
# will run limma K times, once with each snp in snpset
# as the covariate after intercept, adjusting for age, dx year, and
# plate of assay
#
    require(limma)
    fmlas = lapply(paste0(paste0("~", snpset), "+agedx+dxyr+factor(PlateNum)"), 
        as.formula)
    ans = foreach(i=1:length(snpset)) %dopar%  topTable(eBayes(lmFit(eset, 
        model.matrix(fmlas[[i]], data = pData(eset)))), coef = 2, 
        n = ntokeep)
    names(ans) = snpset
    return(ans)
}

runmodelComb = function (eset, snpset, ntokeep = 100)
{
    require(limma)
    fmlas = lapply(paste0(paste0("~", snpset), "*isTumor+agedx+dxyr+factor(PlateNum)+factor(nhsid)"),
        as.formula)
    onemm = model.matrix(fmlas[[1]], data = pData(eset))
    index_interact = grep(":", colnames(onemm))
print(index_interact)
    ans = foreach(i = 1:length(snpset)) %dopar% topTable(eBayes(lmFit(eset,
        model.matrix(fmlas[[i]], data = pData(eset)))), coef = index_interact,
        n = ntokeep)
    names(ans) = snpset
    return(ans)
}

