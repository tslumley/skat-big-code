setwd("~/Documents/SKAT-BIG")
library(survey)
library(CompQuadForm)
wuweights<-function(maf) dbeta(maf,1,25)

source("srht.R")
source("skat-big.R")

makemacsdata<-function(N,length=15000,filter=0.05){
	f<-tempfile()
	system(paste("~/Documents/SKAT-BIG/macs-master/macs",2*N,length," -t .001 -r .001 2>/dev/null | ~/Documents/SKAT-BIG/macs-master/msformatter >", f))
	input<-readLines(f)[-(1:6)]
	unlink(f)
	haplo<-do.call(rbind,lapply(strsplit(input,""),as.integer))
	diplo<-haplo[1:N,]+haplo[(N+1):(2*N),]
	af<-colMeans(diplo)/2
	diplo[,af>0.5]<- 2-diplo[,af>0.5,drop=FALSE]
	maf<-colMeans(diplo)/2
	diplo[,af<=filter,drop=FALSE]
}
all.times<-list()
all.p<-list()

for(i in 1:20){
cat("Simulating\n")	
makemacsdata(5000,length=6e5)->G

cat("Calculating\n")

qchisqsum_satt<-function(x,M){
	tr<-sum(diag(M))
	tr2<-sum(M^2)
	scale<-tr2/tr
	nu<-(tr^2)/tr2
    qchisq(x, ceiling(nu),lower.tail=FALSE)*scale
}


scG<-scale(G)
wuweights<-function(maf) dbeta(maf,1,25)
tildeGt<-t(wuweights(attr(scG,"scaled:center")/2)*t(scG))
all.times$H<-rbind(all.times$H,system.time(H<- t(tildeGt)%*%tildeGt))

Ttail <-qchisqsum_satt(1e-6,H)


all.times$eigenH<-rbind(all.times$eigenH,system.time(e1<-eigen(H,only.values=TRUE,symmetric=TRUE)))
all.times$ssvd0H<-rbind(all.times$ssvd0H,system.time(all.p$ssvd0H<-c(all.p$ssvd0H,pchisqsum_ssvd(Ttail,H,n=100,q=0))))
all.times$ssvd1H<-rbind(all.times$ssvd1H,system.time(all.p$ssvd1H<-c(all.p$ssvd1H,pchisqsum_ssvd(Ttail,H,n=100,q=1))))

all.times$lanH<-rbind(all.times$lanH,system.time(all.p$lanH<-c(all.p$lanH,pchisqsum_partial(Ttail,H,n=100))))
 
all.times$ssvd0G<-rbind(all.times$ssvd0G, system.time(all.p$ssvd0G<-c(all.p$ssvd0G,(pchisqsum_rsvd(Ttail,tildeGt,n=100,q=0,tr2.sample.size=300)))))
all.times$ssvd2G<-rbind(all.times$ssvd2G, system.time(all.p$ssvd2G<-c(all.p$ssvd2G,pchisqsum_rsvd(Ttail,tildeGt,n=100,q=2,tr2.sample.size=300))))


all.times$lanG<-rbind(all.times$lanG,system.time(all.p$lanG<-c(all.p$lanG,pchisqsum_rpartial(Ttail,tildeGt,n=100,tr2.sample.size=300))))

all.times$liu<-rbind(all.times$liu,system.time(all.p$liu<-c(all.p$liu,pchisqsum_liu(Ttail,tildeGt,tr2.sample.size=200))))
all.p$liufull<-c(all.p$liufull,liu(Ttail,e1$values))

 
all.p$saddlefull<-c(all.p$saddlefull,pchisqsum(Ttail,rep(1,length(e1$values)),e1$values,lower.tail=FALSE, method="saddlepoint"))
all.p$intfull<-c(all.p$intfull,davies(Ttail,e1$values,acc=1e-8)$Qq)

all.p$satt<-c(all.p$satt,pchisqsum_partial(Ttail,H,n=0))


library(Matrix)
ww<-wuweights(attr(scG,"scaled:center")/2)
spG<-Matrix(G,sparse=TRUE)%*%Diagonal(x=ww/attr(scG,"scaled:scale"))
cntr<-ww*attr(scG,"scaled:center")/attr(scG,"scaled:scale")

sGmult<-function(X){
	 (t(t(spG%*%X)-colSums(cntr*as.matrix(X))))
	}
	
stGmult<-function(X){
	crossprod(spG,X)- outer(cntr,colSums(as.matrix(X)))
	}	
	
all.times$ssvd2Gs<-rbind(all.times$ssvd2Gs,system.time(all.p$ssvd2Gs<-c(all.p$ssvd2Gs,(pchisqsum_smf(Ttail,sGmult,stGmult,ncol(tildeGt),sum(tildeGt^2))	))))



all.times$lanGs<-rbind(all.times$lanGs,system.time(all.p$lanGs<-c(all.p$lanGs,(pchisqsum_lmf(Ttail,sGmult,stGmult,ncol(tildeGt),nrow(tildeGt),sum(tildeGt^2))))))

print(i)
}



	