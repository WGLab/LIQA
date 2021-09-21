##########
## MAIN ##
##########
library("gcmr")
library("betareg")

input.file <- commandArgs(trailingOnly = TRUE)
#raw.data <- read.table("exonresult")
#raw.data <- read.table("isoest")

raw.data <- read.table(input.file[1])
gene.list <- unique(raw.data$V1)
generesultfile <- input.file[2]
exonresultfile <- input.file[3]

## set columns of data frame
test.result <- NA
gene.p.value.list <- c()
exon.p.value.list <- c()
result.gene.col <- c()
result.exon.col <- c()

##########################
## START TO DO ANALYSIS ##
##########################
### hide warning messages
options(warn = -1)

for(gene.name in gene.list) {
  
  ## get data for within the gene
  DS.tmp.list <- levels(factor(raw.data[raw.data$V1==gene.name,]$V2))
  data.gene.frame <- data.frame(factor(raw.data[raw.data$V1==gene.name,]$V2),factor(raw.data[raw.data$V1==gene.name,]$V3),raw.data[raw.data$V1==gene.name,]$V4,raw.data[raw.data$V1==gene.name,]$V6)
  colnames(data.gene.frame) <- c("DS","gp","est","spid")
  data.gene.frame <- data.gene.frame[order(data.gene.frame$spid),]
  numDS <- length(DS.tmp.list)
  
  #####################################################
  ## construct new design matrix for exon based test ##
  #####################################################
  allgp <- paste(data.gene.frame$DS,data.gene.frame$gp,sep=":")
  data.gene.frame <- data.frame(data.gene.frame,allgp)
  tmpgp <- allgp
  for(DS.tmp in DS.tmp.list){tmpgp[ tmpgp==paste(DS.tmp,1,sep=":") | tmpgp==paste(DS.tmp,2,sep=":") ] = paste(DS.tmp,"H0",sep="")}
  data.gene.frame <- data.frame(data.gene.frame,tmpgp)
  colnames(data.gene.frame)[6] <- "allgpH0"
  for(i in 1:length(DS.tmp.list)){
    DS.tmp <- DS.tmp.list[i]
    tmpgp <- allgp
    tmpgp[ tmpgp==paste(DS.tmp,1,sep=":") | tmpgp==paste(DS.tmp,2,sep=":") ] = "H0"
    data.gene.frame <- data.frame(data.gene.frame,tmpgp)
    tmpidx <- match(DS.tmp,DS.tmp.list)+6
    tmpidx <- i + 6
    colnames(data.gene.frame)[tmpidx] <- DS.tmp
  }
  
  ## reset value on boundary
  if(sum(data.gene.frame$est>.99999)>0) {data.gene.frame[data.gene.frame$est>.99999,]$est <- .99999}
  if(sum(data.gene.frame$est<.00001)>0) {data.gene.frame[data.gene.frame$est<.00001,]$est <- .00001}
  if(length(unique(data.gene.frame[data.gene.frame$gp==1,]$est)) <3 ) {next}
  if(length(unique(data.gene.frame[data.gene.frame$gp==2,]$est)) <3 ) {next}
  cat(paste(c("DAS testing: Gene", gene.name, "being processed...")))
  
  ###########################
  ## start to perform test ##
  ###########################
  if(numDS > 1) {
    ########################
    ## do gene based test ##
    ########################
    fit.H1 <- try(gcmr( est ~ allgp , marginal = beta.marg, cormat = cluster.cormat( id = spid, type = "exchangeable" ), data = data.gene.frame ))
    fit.H0 <- try(gcmr( est ~ allgpH0 , marginal = beta.marg, cormat = cluster.cormat( id = spid, type = "exchangeable" ), data = data.gene.frame ))
    if(grepl("Error",fit.H0)==T){next}
    if(grepl("Error",fit.H1)==T){next}
    LR.stat <- -2*(fit.H1$maximum-fit.H0$maximum)
    gene.p.value.list <- c(gene.p.value.list, 1-pchisq(LR.stat,numDS))
    result.gene.col <- c(result.gene.col,gene.name)
    ########################
    ## do exon based test ##
    ########################
    #for(i in 1:numDS) {
    #  tmpcol <- i+6
    #  fit.H1 <- try(gcmr( est ~ allgp , marginal = beta.marg, cormat = cluster.cormat( id = spid, type = "exchangeable" ), data = data.gene.frame ))
    #  fit.H0 <- try(gcmr( est ~ data.gene.frame[,tmpcol] , marginal = beta.marg, cormat = cluster.cormat( id = spid, type = "exchangeable" ), data = data.gene.frame ))
    #  if(grepl("Error",fit.H0)==T){next}
    #  if(grepl("Error",fit.H1)==T){next}
    #  LR.stat <- -2*(fit.H1$maximum-fit.H0$maximum)
    #  exon.p.value.list <- c(exon.p.value.list, 1-pchisq(LR.stat,1))
    #  result.exon.col <- c(result.exon.col,paste(gene.name,DS.tmp.list[i],sep=":"))
    #}
    
  }
  if(numDS == 1) {
    ###########################################################
    ## do gene based test and exon based test simultaneously ##
    ###########################################################
    fit.H1 <- try( betareg( est ~ allgp , data = data.gene.frame ) )
    fit.H0 <- try( betareg( est ~ 1 , data = data.gene.frame ) )
    if(grepl("Error",fit.H0)==T){next}
    if(grepl("Error",fit.H1)==T){next}
    LR.stat <- 2*(fit.H1$loglik-fit.H0$loglik)
    gene.p.value.list <- c(gene.p.value.list, 1-pchisq(LR.stat,1))
    result.gene.col <- c(result.gene.col,gene.name)
    exon.p.value.list <- c(exon.p.value.list, 1-pchisq(LR.stat,1))
    result.exon.col <- c(result.exon.col, paste(gene.name,DS.tmp.list[1],sep=":") )
  }

cat("Done!\n")
}

# summarize gene based results
test.result.gene <- data.frame(result.gene.col,gene.p.value.list)
if(is.null(test.result.gene) == F){
colnames(test.result.gene) <- c("gene","pvalue")
}
write.table(test.result.gene,file=generesultfile,row.names = FALSE,quote=F, sep="\t", col.names=F)

# summarize exon based results
#exon.loc <- c()
#for(i in result.exon.col) {
  #print(i)
#  tmp <- unlist(strsplit(i,":"))
#  tmploc <- factor(raw.data[ raw.data$V1==tmp[1] & raw.data$V2==tmp[2], ]$V5)
#  exon.loc <- c(exon.loc,levels(tmploc))
#}
#result.exon.col.tmp <- unlist(strsplit(result.exon.col,":"))
#result.exon.col.tmp <- matrix(result.exon.col.tmp,ncol=2,byrow=T)
#test.result.exon <- data.frame(result.exon.col.tmp[,1],result.exon.col.tmp[,2],exon.p.value.list,exon.loc)
#colnames(test.result.exon) <- c("gene","DS","pvalue","location")
#write.table(test.result.exon,file=exonresultfile,row.names = FALSE,quote=F, sep="\t", col.names=F)
