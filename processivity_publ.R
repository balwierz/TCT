library(purrr)
library(glmnet)
library(Hmisc)
library(openxlsx)
options(stringsAsFactors = FALSE)

# This script requires the following files:
#  HostGeneMatchingStagesSnoRNATPM.txt
#  di.genes.nAnTi.txt


rowNormalise <- function(mat)
{
    means <- rowMeans(mat)
    subtraction <- rep(means, dim(mat)[2]) %>% matrix(nrow = dim(mat)[1], ncol = dim(mat)[2])
    mat - subtraction
}

A <- read.table ("HostGeneMatchingStagesSnoRNATPM.txt", header=T, stringsAsFactors = F)

yr <- A[,c(2:8)]
yc <- A[,c(9:15)]
sno <- A[,c(16:22)]

# remove bad quality (snoRNA) 64C sample and CAGE along:
yr$YR_64cells <- NULL
yc$YC_64cells <- NULL
sno$snorna_64cells <- NULL

# differential usage
du = (yc / (yr + yc + 0.00000001)) # differential usage of initiators
apply(du, MARGIN = 1, FUN = median) %>% Ecdf; abline(v=0.75); abline(v=0.25)     # median instead of 24h, because in the next replicate nanog has zeros there
whichYcDom = apply(du, MARGIN = 1, FUN = median) >= 0.75
whichYrDom = apply(du, MARGIN = 1, FUN = median) <= 0.25

# sno data for the model.
# pseudocount, log-transform, normalise, form in one vector sno1,t1 sno1,t2, ... sno2,t1, ...
snoPc <- 5.0  # 5 tpm pseudocount for snoRNA
cagePc <- 5.0 # 5 tpm pseudocount for CAGE
Y <- sno %>% `+`(snoPc) %>% log() %>% rowNormalise() %>% t %>% as.list %>% do.call(what=rbind)

X <- data.frame(YR = yr %>% `+`(cagePc) %>% log2() %>% rowNormalise() %>% t %>% as.list %>% do.call(what=rbind),
                YC = yc %>% `+`(cagePc) %>% log2() %>% rowNormalise() %>% t %>% as.list %>% do.call(what=rbind)) %>%
    as.matrix()

whichYrDom.expanded <- rep(whichYrDom, each=dim(yr)[2])
whichYcDom.expanded <- rep(whichYcDom, each=dim(yr)[2])

# scatter plots:

{svg("yr_sno.scatter.svg", width = 5, height = 5)
plot(X[whichYrDom.expanded, "YR"], Y[whichYrDom.expanded], xlab="YR initiation log2 fold change",
     ylab="snoRNA log2 fold change",
     col="blue", cex=0.7,
     xlim=c(-4, 4), ylim=c(-3, 3))
points(X[whichYcDom.expanded, "YR"], Y[whichYcDom.expanded], col="red", cex=0.7)
points(X[!whichYcDom.expanded & !whichYrDom.expanded, "YR"], 
       Y[!whichYcDom.expanded & !whichYrDom.expanded], col="black", cex=0.7)
legend("topleft", bty = "n", paste0("r = ", sprintf(fmt = "%0.2f", cor(X[,"YR"], Y))))
dev.off()}


plot(X[, "YC"], Y, xlab="YC", ylab="sno")
{svg("yc_sno.scatter.svg", width = 5, height = 5)
plot(  X[whichYrDom.expanded, "YC"], Y[whichYrDom.expanded], xlab="YC initiation log2 fold change",
       ylab="snoRNA log2 fold change", 
       col="blue", cex=0.7,
       xlim=c(-4, 4), ylim=c(-3, 3))
points(X[whichYcDom.expanded, "YC"], Y[whichYcDom.expanded], col="red", cex=0.7)
points(X[!whichYcDom.expanded & !whichYrDom.expanded, "YC"], 
       Y[!whichYcDom.expanded & !whichYrDom.expanded], col="black", cex=0.7)
legend("topleft", bty = "n", paste0("r = ", sprintf(fmt = "%0.2f", cor(X[,"YC"], Y))))
dev.off()}


FIT.cv <- cv.glmnet(x=X, y=Y, intercept=F, alpha=1, standardize=F)
plot.cv.glmnet(FIT.cv)
coef.glmnet(FIT.cv, s=FIT.cv$lambda.min)
coef.glmnet(FIT.cv, s=FIT.cv$lambda.1se)

FIT = glmnet(x=X, y=Y, intercept=F, alpha=1, standardize=F)
plot(FIT, col=c("red", "blue"))

{svg("coefficients.svg", width = 5, height = 5)
plot.glmnet(FIT, xvar="lambda", col=c("blue", "red"), lwd=2)
grid()
legend("center", legend=c("YC coefficient", "YR coefficient"), col=c("red", "blue"), lwd=2, bty="n")
abline(v=log(FIT.cv$lambda.min), lty=2)
dev.off()}

{svg("coefficients_1se.svg", width = 5, height = 5)
    plot.glmnet(FIT, xvar="lambda", col=c("blue", "red"), lwd=2)
    grid()
    legend("center", legend=c(expression(""*A[YC]*" coefficient"),
                              expression(""*A[YR]*" coefficient")), col=c("red", "blue"), lwd=2, bty="n")
    abline(v=log(FIT.cv$lambda.1se), lty=2)
    dev.off()}



#### The second part, nAnTi CAGE:
nanti <- read.table("di.genes.nAnTi.txt", header = T)

rep1.IDs <- A$geneId %>% strsplit(":", T) %>% map(`[`, 1) %>% unlist()
rep2.IDs <- nanti$geneId %>% strsplit("_", F) %>% map(`[`, 1) %>% unlist()  # non-unique

yr.nanti <- nanti[match(rep1.IDs, rep2.IDs), 2:8]
yc.nanti <- nanti[match(rep1.IDs, rep2.IDs), 9:15]
yr.nanti$S02_128Cells <- NULL
yr.nanti$S03_512Cells <- NULL
yc.nanti$S02_128Cells.1 <- NULL
yc.nanti$S03_512Cells.1 <- NULL

# matching stages of nanti in the old snoRNA expression:
sno.nanti <- sno[, c("snorna_mergedEgg", 'snorna_30epiboly', 'snorna_12somite', 'snorna_prim5', 'snorna_prim16')]

# and remove the NAs:
sno.nanti <- sno.nanti[is.na(yr.nanti) %>% rowSums() == 0, ]
yr.nanti <- yr.nanti[is.na(yr.nanti) %>% rowSums() == 0, ]
yc.nanti <- yc.nanti[is.na(yc.nanti) %>% rowSums() == 0, ]

X.nanti <- data.frame(YR = yr.nanti %>% `+`(cagePc) %>% log2() %>% rowNormalise() %>% t %>% as.list %>% do.call(what=rbind),
                      YC = yc.nanti %>% `+`(cagePc) %>% log2() %>% rowNormalise() %>% t %>% as.list %>% do.call(what=rbind)) %>%
    as.matrix()
Y.nanti <- sno.nanti %>% `+`(snoPc) %>% log() %>% rowNormalise() %>% t %>% as.list %>% do.call(what=rbind)

du.nanti = (yc.nanti / (yr.nanti + yc.nanti+ 0.000000001)) # differential usage of initiators. warning: would be division by zero in nanog trascript!
apply(du.nanti, MARGIN = 1, FUN = median) %>% Ecdf; abline(v=0.75); abline(v=0.25)
whichYcDom.nanti = apply(du.nanti, MARGIN = 1, FUN = median) >= 0.75
whichYrDom.nanti = apply(du.nanti, MARGIN = 1, FUN = median) <= 0.25

whichYrDom.nanti.expanded <- rep(whichYrDom.nanti, each=dim(yr.nanti)[2])
whichYcDom.nanti.expanded <- rep(whichYcDom.nanti, each=dim(yr.nanti)[2])

svg("yr_sno.scatter.nanti.svg", width = 5, height = 5)
    plot(X.nanti[whichYrDom.nanti.expanded, "YR"], Y.nanti[whichYrDom.nanti.expanded], xlab="YR initiation log2 fold change",
         ylab="snoRNA log2 fold change",
         col="blue", cex=0.7,
         xlim=c(-4, 4), ylim=c(-3, 3)
         )
    points(X.nanti[whichYcDom.nanti.expanded, "YR"], Y.nanti[whichYcDom.nanti.expanded], col="red", cex=0.7)
    points(X.nanti[!whichYcDom.nanti.expanded & !whichYrDom.nanti.expanded, "YR"], 
           Y.nanti[!whichYcDom.nanti.expanded & !whichYrDom.nanti.expanded], col="black", cex=0.7)
    legend("topleft", bty = "n", paste0("r = ", sprintf(fmt = "%0.2f", cor(X.nanti[,"YR"], Y.nanti))))
    dev.off()


{svg("yc_sno.scatter.nanti.svg", width = 5, height = 5)
    plot(  X.nanti[whichYrDom.nanti.expanded, "YC"], Y.nanti[whichYrDom.nanti.expanded], xlab="YC initiation log2 fold change",
           ylab="snoRNA log2 fold change", 
           col="blue", cex=0.7,
           xlim=c(-4, 4), ylim=c(-3, 3))
    points(X.nanti[whichYcDom.nanti.expanded, "YC"], Y.nanti[whichYcDom.nanti.expanded], col="red", cex=0.7)
    points(X.nanti[!whichYcDom.nanti.expanded & !whichYrDom.nanti.expanded, "YC"], 
           Y.nanti[!whichYcDom.nanti.expanded & !whichYrDom.nanti.expanded], col="black", cex=0.7)
    legend("topleft", bty = "n", paste0("r = ", sprintf(fmt = "%0.2f", cor(X.nanti[,"YC"], Y.nanti))))
    dev.off()}

# fitting:

FIT.nanti.cv <- cv.glmnet(x=X.nanti, y=Y.nanti, intercept=F, alpha=1, standardize=F)
plot.cv.glmnet(FIT.nanti.cv)

coef.glmnet(FIT.nanti.cv, s=FIT.nanti.cv$lambda.min)
coef.glmnet(FIT.nanti.cv, s=FIT.nanti.cv$lambda.1se)

FIT.nanti = glmnet(x=X.nanti, y=Y.nanti, intercept=F, alpha=1, standardize=F)
plot(FIT.nanti, col=c("red", "blue"))

{svg("coefficients.nanti.svg", width = 5, height = 5)
    plot.glmnet(FIT.nanti, xvar="lambda", col=c("blue", "red"), lwd=2)
    grid()
    legend("center", legend=c(expression(""*A[YC]*" coefficient"),
                              expression(""*A[YR]*" coefficient")), col=c("red", "blue"), lwd=2, bty="n")
    abline(v=log(FIT.nanti.cv$lambda.1se), lty=2)
    dev.off()}

