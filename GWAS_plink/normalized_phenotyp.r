#---
#to convert distribution of quantative traits into guassian distribution.
#---

rmoutliner<-function(x){
  qnt <- quantile(x, probs=c(.25, .75), na.rm = T)
  str(qnt)
  H <- 1.5 * IQR(x, na.rm = T)
  y <- x
  str(H)
  check1_t=qnt[1] - H
  check1<-y[x < (qnt[1] - H) & !is.na(x)]
  check2<-y[x > (qnt[2] + H) & !is.na(x)]
  check2_t=qnt[2] + H
  print("outliners_lower")
  str(check1_t)
  str(check1)
  print("outliners_upper")
  str(check2_t)
  str(check2)
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  x_na<-x[!is.na(x)]
  print("total x")
  str(x_na)
  y_na<-y[!is.na(y)]
  print("total y after remove NA")
  str(y_na)
  y
}

qtrans<-function(x){ #if you would like to remove outliner or not
   k<-which(!is.na(x))
   ran<-rank(as.numeric(x[k]))
   y<-qnorm((1:length(k)-0.5)/length(k))
   x[k]<-y[ran]
   str(y)
   str(x)
   x  
}

args=commandArgs(T)
file=args[1]
remove_ol=args[2] #it is the code for missing value
inmiss=args[3] #input missing code
outmiss=args[4] #output missing code
outdir=args[5]
print(paste("inputfile is",file))
#---
#reading data
#---
d=read.table(file,header=F,stringsAsFactors=F,comment='')
data_ori=d
#data_ori=d[,-ncol(d)]
data=as.numeric(data_ori[,ncol(data_ori)])
data[data==inmiss]="NA"
data=as.numeric(data)

if(remove_ol=="Y"){
   data_rout<-rmoutliner(data)
   outfile1=paste0(outdir,"/",basename(file),".rmout")
   outfile2=paste0(outdir,"/",basename(file),".rmout.qtrans")
   print("data_rout")
   str(data_rout)
   print("outfile1")
   print(outfile1)
   data_rout_print=cbind(data_ori[,-ncol(data_ori)],data_rout)
   write.table(data_rout_print,file=outfile1,quote=F,sep="\t",row.names=F,col.names=F)
}else{
   data_rout<-data
   str(data_rout)
   outfile1=paste0(outdir,"/",basename(file),".normout")
   outfile2=paste0(outdir,"/",basename(file),".normout.qtrans")
   write.table(data_rout,file=outfile1,quote=F,sep="\t",row.names=F,col.names=F)
}

data_normalized=qtrans(data_rout)
data_nout=data_normalized
data_nout[is.na(data_nout)]=outmiss
data_out=cbind(data_ori[,-ncol(data_ori)],data_nout)
#data_out=cbind(data_out,d[,ncol(d)])
print("data_out")
str(data_out)
print("outfile2")
print(outfile2)

write.table(data_out,file=outfile2,quote=F,sep="\t",row.names=F,col.names=F)

#---
#plot
#---
data_plot1=data[data!="-9" & data!="-999" & !is.na(data)]
data_plot2=data_rout[data_rout!="-9" & data_rout!="-999" & !is.na(data_rout)]
data_plot3=data_normalized[data_normalized!="-9" & data_normalized!="-999" & !is.na(data_normalized)]

print("data_plot1:")
str(data_plot1)
print("data_plot2:")
str(data_plot2)
print("data_plot3:")
str(data_plot3)


pdf(paste0(outdir,"/",basename(file),".checkroutnormalization.pdf"),hei=10,wid=15)
#png(paste0(outdir,"/",basename(file),".checkroutnormalization.png"),wid=960,hei=480)
par(mfrow=c(2,3))
hist(data_plot1,xlab="Before Rmout and Qtrans",main=paste(basename(file),paste0("N=",length(data_plot1)),sep="\n"))
hist(data_plot2,xlab="After Rmout",main=paste(basename(file),paste0("N=",length(data_plot2)),sep="\n"))
hist(data_plot3,xlab="After Rmout and Qtras",main=paste(basename(file),paste0("N=",length(data_plot3)),sep="\n"))

qqnorm(data_plot1,xlab="Before Rmout and Qtrans",main=basename(file))
qqline(data_plot1)
qqnorm(data_plot2,xlab="After Rmout",main=basename(outfile1))
qqline(data_plot2)
qqnorm(data_plot3,xlab="After Rmout and Qtrans",main=basename(outfile2))
qqline(data_plot3)
dev.off()

