a <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(vegan))
library(permute)
L <- as.matrix(read.delim(a[1],header=FALSE,check.names=FALSE)); P <- as.matrix(read.delim(a[2],header=FALSE,check.names=FALSE)); m <- read.delim(a[3],check.names=FALSE); out <- a[4]; nperm <- as.integer(a[5]); set.seed(as.integer(a[6])); g <- factor(m$group,levels=unique(m$group)); if(any(table(g)<2)) stop('Need at least two patients per group')
perms <- shuffleSet(length(g), nset=nperm)
res <- lapply(list(fragment_length=L,midpoint=P),function(x){d<-dist(x); if(max(d)==0) stop('All patient profiles identical'); r<-adonis2(d~g,permutations=perms); bd<-betadisper(d,g,type='median',bias.adjust=TRUE); z<-permutest(bd,permutations=perms); data.frame(F=r$F[1],R2=r$R2[1],P=r$`Pr(>F)`[1],PERMDISP_P=z$tab$`Pr(>F)`[1])});result<-do.call(rbind,res);result$feature<-rownames(result);rownames(result)<-NULL;result$BH_q_BC<-p.adjust(result$P,'BH');write.table(result,file.path(out,'distribution_tests.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
if('central_coverage' %in% names(m)) {if(anyNA(m$central_coverage)) stop('central_coverage missing'); z<-adonis2(dist(m$central_coverage)~g,permutations=perms);write.table(data.frame(test='Griffin central ±30bp permutation ANOVA',N=nrow(m),F=z$F[1],R2=z$R2[1],P=z$`Pr(>F)`[1]),file.path(out,'central_coverage_test.tsv'),sep='\t',row.names=FALSE,quote=FALSE)}
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
