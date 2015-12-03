library(gRapHD)
library(glasso)
# load the data
setwd("D://fall 2015/sta250/project")
load('LOI_data_414x1261.Rdata')
Y = data.frame(Y)

# Chapter 7 
# High Dimensional Modelling
dim(Y)
table(sapply(Y, class))

#Minimal BIC forest
bF <- minForest(Y)
plot(bF) 

# bad ideal to plot the whole HD graph. No structure visible
# interested in tamoxifen treatment and gene expression
nby <- neighbourhood(bF, orig=1260, rad=5)$v[,1]
plot(bF, vert=nby, numIter=1000, plotVert=TRUE, cex.vert.label=0.7,
     lwd.vert = 2,vert.hl=nby[1], col.hl = 'yellow',
     vert.radii = 0.02)

# The plot use the iterative layout algorithm of Fruchterman and Reingold.

# HHex inhibits the migration of breast and prostate epithelial cells through 
# direct transcriptional regulation of Endoglin.
HHEX = Y[,which(colnames(Y)=='HHEX')]
er = round(Y[,which(colnames(Y)=='er')])
mean(HHEX[er==0])
mean(HHEX[er==-2])
mean(HHEX[er==-4])
# the relation between er and HHEX(er-positive invasiv)

# plot components
# table(Degree(bF))
# col = c('red','blue','yellow','black','green')
# plot(bF, numIter=10000, vert.labels=1:bF@p, main='min BIC forest',
#      border=col[Degree(bF)], lwd.vert=2,vert.radii=0.005)


# decomposable stepwise search 
# This can lead to considerable efficientcy gains
bc.marg <- Y[,nby]
mbF <- minForest(bc.marg)
plot(mbF,numIter=1000,vert.hl=1,col.hl='yellow',
     cex.vert.label=0.7,lwd.vert = 2,vert.radii = 0.02)

# The gene adjacent to er is HHEX as we saw before. This suggests that the effect
# of the er status may show it's effect on gene expression of gen HHEX.

mbG <-stepw(model=mbF, data=bc.marg)
posn <- plot(mbF,numIter=1000,disp=F)
plot(mbG,numIter=0, coord=posn,vert.hl=1,col.hl='yellow',
     cex.vert.label=0.7,lwd.vert = 2,vert.radii = 0.02)

# Although the minimal BIC decomposable model is considerably less soarse, the 
# interpretation is unaltered.


# Selection by Approximation
# Using graphical lasso algorithm
S <- cor(Y[,nby[-1]])
res.lasso <- glasso(S, rho=0.45)   # rho can be changed
AM <- res.lasso$wi != 0
diag(AM) <- FALSE
rownames(AM) <- colnames(AM) <- names(Y)[nby[-1]]
g.lasso <- as(AM, 'graphNEL')
g.lasso

g.HD <- as(g.lasso, 'gRapHD')
plot(g.HD, numIt=1000,vert.hl=1,col.hl='yellow',
     cex.vert.label=0.7,lwd.vert = 2,vert.radii = 0.02)

# Above example selects a model to the variables in the neighbourhood of the er variable in
# the breast cancer dataset. (omitting the er variable itself since it is discrete)
# We choose parameter rho=0.45. This choice was made so as to obtain a grapg of comparable density
# to those obtained previously.