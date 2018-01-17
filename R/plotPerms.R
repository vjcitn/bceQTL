
plotPerms = function() {
data(allpermFDR)
allpermFDR[which(allpermFDR$permfdr < 1e-6), "permfdr"] = 5e-3 # don't let 0 permutation FDR spoil log-log
with(allpermFDR, plot(limfdr, permfdr, log="xy", xlim=c(1e-4,.15),
   xlab="FDR by limma", ylab="FDR by 501 permutations of adjusted expr. against genotype", col=factor(type), pch=19, axes=FALSE))
axis(2)
axis(1, at=c(1e-4, 5e-4, 5e-3, 5e-2, 1e-1), labels=c(".0001", ".0005", ".005", ".05", ".1"))
abline(0,1)
abline(h=.1, lty=2)
bi= which(allpermFDR$permfdr>.1)
bad = allpermFDR[bi,]
with(bad[bad$permfdr>.2,], text(limfdr, permfdr+.04, as.character(snp), cex=.6))
with(bad[bad$permfdr < .2 & bad$permfdr > .1,], text(limfdr, permfdr+.01, as.character(snp), cex=.6))
legend(1e-4, .25, pch=19, legend=levels(factor(allpermFDR$type)), col=
  as.numeric(factor(levels(factor(allpermFDR$type)))))
abline(v=.1, lty=2)
}
