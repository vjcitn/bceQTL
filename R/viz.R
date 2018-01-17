regressOut = function (sms, rhs, ...) 
{
    if (is(sms, "smlSet") | is(sms, "ExpressionSet")) {
        mm = model.matrix(rhs, data = pData(sms))
        ex = exprs(sms)
    }
    else if (is(sms, "SummarizedExperiment") || is(sms, "RangedSummarizedExperiment")) {
        mm = model.matrix(rhs, data = colData(sms))
        message("using assay() to extract 'expression' matrix from SummarizedExperiment")
        ex = assay(sms)
    }
    else stop("only works for smlSet, ExpressionSet or SummarizedExperiment")
    f = limma::lmFit(ex, mm, ...)
    r = ex - (f$coef %*% t(f$design))
    if (is(sms, "smlSet") | is(sms, "ExpressionSet")) 
        sms@assayData = assayDataNew("lockedEnvironment", exprs = r)
    else if (is(sms, "SummarizedExperiment") || is(sms, "RangedSummarizedExperiment")) 
        assay(sms) = r
    sms
}

.getinfo = function(...)
  gsub("(.*)([-,])(.*)", "\\1.\\3", getinfo(...))

getinfo = function(tcno) {
 rownames(tcann) = as.character(tcann[,1])
 tmp = tcann[tcno,]$public_transcript_id
 tmp = strsplit(tmp, " /// ")[[1]]
 select(EnsDb.Hsapiens.v75, keys=tmp, keytype="TXID", columns="GENENAME")$GENENAME[1]
}

prebox = function(es, snp, gene) {
 snpd = as.character(round(pData(es)[[snp]],0))
 ed = as.numeric(exprs(es)[gene,])
 ans = data.frame(ex=ed, snpd=snpd, stringsAsFactors=FALSE)
 gene = paste0(gene, "..",  .getinfo(gene))
 names(ans) = c(gene, snp)
 ans
}

triptych = function(snp, gene, negES, posES, normES) {
neg1 = cbind(prebox(negES, snp, gene), type="ER-") # "rs17529111", "TC0400338"),type="ER-")
pos1 = cbind(prebox(posES, snp, gene), type="ER+") #"rs17529111", "TC0400338"),type="ER+")
norm1 = cbind(prebox(normES, snp, gene), type="AdjNorm") #"rs17529111", "TC0400338"),type="ER+")
full1 = rbind(neg1, pos1, norm1)
gene = paste0(gene, "..",  .getinfo(gene))
p1 = ggplot(full1, aes_string(y=gene, x=snp)) + geom_boxplot() + facet_grid(.~type)
p1 
}

quadpanel = function(snp, gene, negES, posES, negnormES, posnormES, extraX="") {
labx = paste0(snp, extraX)
neg1 = cbind(prebox(negES, snp, gene), type="ER-") 
pos1 = cbind(prebox(posES, snp, gene), type="ER+") 
negnorm1 = cbind(prebox(negnormES, snp, gene), type="AdjNormER-") 
posnorm1 = cbind(prebox(posnormES, snp, gene), type="AdjNormER+") 
full1 = rbind(neg1, pos1, negnorm1, posnorm1)
gene = paste0(gene, "..",  .getinfo(gene))
p1 = ggplot(full1, aes_string(y=gene, x=snp)) + geom_boxplot() + 
   facet_grid(.~type) + xlab(labx)
p1 
}

geneLabelDF = function(syms) {
 require(Homo.sapiens)
 chrs = mapIds(Homo.sapiens, keys=syms, keytype="SYMBOL", column="TXCHROM")
 locs = mapIds(Homo.sapiens, keys=syms, keytype="SYMBOL", column="TXSTART")
 ans = data.frame(chrs, as.numeric(locs), as.numeric(locs), syms)
 ans = na.omit(ans)
 colnames(ans) = c("chr", "start", "start2", "sym")
 ans
}

 getglo = function(sym) {
   stopifnot(sym %in% keys(Homo.sapiens, keytype="SYMBOL"))
   mapIds(Homo.sapiens, keys=sym, keytype="SYMBOL", column="TXSTART")
   }
 getgchr = function(sym) {
   stopifnot(sym %in% keys(Homo.sapiens, keytype="SYMBOL"))
   mapIds(Homo.sapiens, keys=sym, keytype="SYMBOL", column="TXCHROM")
   }

l2df = function(x, coln=c("snpid", "geneloc")) {
  nrep = sapply(x,length)
  nn = rep(names(x), nrep)
  data=unlist(x)
  ans = data.frame(nn, data, stringsAsFactors=FALSE)
  names(ans) = coln
  ans
}

filterSyms = function(sgl, odb) {
 k = keys(odb, keytype="SYMBOL")
 sgl2 = lapply(sgl, function(x) intersect(x, k))
 names(sgl2) = names(sgl)
 le = sapply(sgl2,length)
 sgl2[which(le>0)]
}

sglistToDF = function(sgl, snpgr=NULL, slocpack=SNPlocs.Hsapiens.dbSNP144.GRCh37, filterSyms=TRUE) {
 if (filterSyms) sgl = filterSyms(sgl, Homo.sapiens)
 sn = names(sgl)
 if (is.null(snpgr)) snlocs = snpsById( slocpack, sn, ifnotfound="drop" )
 else {  # comply with snpsById protocol as much as possible
   snlocs = snpgr[sn]
   snlocs$RefSNP_id = names(snlocs)
   }
 avail = snlocs$RefSNP_id
 sgl = sgl[avail]  # reorder list if necessary
 require(GenomeInfoDb)
 seqlevelsStyle(snlocs) = "UCSC"
# snpchr = gsub("ch", "chr", names(snlocs))

 if (is(snlocs, "GPos")) sl=pos(snlocs)
 else sl = start(snlocs)
 sldf = data.frame(snpid=names(sgl), snpchr=as.character(seqnames(snlocs)), snploc=sl, stringsAsFactors=FALSE)
 genechr = l2df(lapply(sgl, getgchr), c("snpid", "genechr"))
 geneloc = l2df(lapply(sgl, getglo), c("snpid", "geneloc"))
 gdf = cbind(genechr, geneloc[,-1])
 names(gdf) = c("snpid", "gchr", "gloc")
 merge(sldf, gdf, by="snpid")
}

