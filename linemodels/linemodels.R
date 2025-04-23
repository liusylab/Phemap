library(data.table)
library(openxlsx)
library(tidyverse)
library(linemodels)
rep<-read.xlsx('rep_results_all.xlsx') #comparison results with Taiwan Biobank GWAS data

for(pheno in unique(rep$Traits)){
  data<-rep[rep$Traits==pheno,]
  X = data[,c("BETA_meta","b_tw")]
  SE = data[,c("SE_meta","se_tw")]
  scales = c(0.2, 0.2)
  slopes = c(0,1)
  cors = c(0.999, 0.999)
  model.names=c('pregnancy','general')
  
  par.include=rbind(c(TRUE,FALSE,FALSE),
                    c(TRUE,TRUE,FALSE))
  set.seed(91)
  res<-line.models.optimize(X,SE,par.include = par.include,
                            init.scales = scales,
                            init.slopes = slopes,
                            init.cors = cors,
                            model.names = model.names,
                            model.priors = c(0.5,0.5),
                            r.lkhood = 0,print.steps = 2,tol.loglk = 1e-2)
  
  scales=c(res[1],res[2])
  slopes=c(0,res[4])
  cors=c(0.999,0.999)
  
  res.lm = line.models(X, SE,
                       scales = scales,
                       slopes = slopes,
                       cors = cors,
                       model.names = model.names,
                       model.priors = c(0.5,0.5),
                       r.lkhood = 0)
  colors = c("red","dodgerblue")
  cols = colors[apply(res.lm, 1, which.max)]
  ind = (apply(res.lm, 1, max) > 0.95)
  pdf(paste0(pheno,'.pdf'),width = 5.5, height = 5.5)
  visualize.line.models(scales, slopes, cors,
                        model.names = model.names, model.cols = colors,
                        legend.position = "topleft",
                        xlab = expression(beta[pregnancy]), ylab = expression(beta[TW]),
                        emphasize.axes = FALSE, cex.lab = 1.1, cex.axis = 1.1)
  points(X, pch = 1, col = "black", cex = 0.6,lwd=1)
  points(X[ind,], pch = 19, col = cols[ind], cex = 0.5)
  cols[!ind]='other'
  b1<-data.frame(res.lm,cols)
  data<-cbind(data,b1)
  dev.off()
  write.table(data,paste0(pheno,'_results.txt'),sep = '\t',quote = F,row.names = F)
  write.table(res,paste0(pheno,'_param.txt'),sep = '\t',quote = F,row.names = F)
}


