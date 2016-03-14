##
##  x: the statistic value
##  M: the matrix whose eigenthingies we need
##  n: how many largest eigenvalues to use (time roughly proportion to this)
##  method: what method to use for computing the distribution based on the exact eigenvalues
##         (the 'integration' method has issues for very large matrices)


library(svd)

pchisqsum_partial<-function(x,M,n=100,method=c("saddlepoint","integration"),tol=1e-3){
	method<-match.arg(method)
	if (n>0){
		ee<-svd::trlan.eigen(M,n,opts=list(tol=tol))
	} else {
		ee<-list(values=numeric(0))
	}
	tr<-sum(diag(M))
	tr2<-sum(M^2)
	tr.small<-tr-sum(ee$d)
	tr2.small<-tr2-sum(ee$d^2)
	scale<-tr2.small/tr.small
	nu<-(tr.small^2)/tr2.small
    pchisqsum(x, c(rep(1,n),ceiling(nu)), c(ee$d, scale), method=method,lower.tail=FALSE)
}

pchisqsum_rpartial<-function(x,M,n=100, tr2.sample.size=300, method=c("saddlepoint","integration")){
	method<-match.arg(method)
	ee<-svd::trlan.svd(M,neig=n)$d[1:n]^2
	diags <- colSums(M^2)
	tr<-sum(diags)
	if (tr2.sample.size>0){
		tr2<-tracefht(M,k=tr2.sample.size,trace.full=tr)
	} else {
		Ms<-crossprod(M)
		tr2<- sum(Ms^2)
	}	
	tr.small<-tr-sum(ee)
	tr2.small<-tr2-sum(ee^2)
	scale<-tr2.small/tr.small
	nu<-(tr.small^2)/tr2.small
    pchisqsum(x, c(rep(1,n), ceiling(nu)), c(ee, scale), method=method,lower.tail=FALSE)
}	



pchisqsum_lmf<-function(x,Mmult,tMmult,ncolM,nrowM, traceM, n=100, tr2.sample.size=300, method=c("saddlepoint","integration")){
	method<-match.arg(method)
	extM<-extmat(function(x) as.numeric(Mmult(x)), function(x) as.numeric(tMmult(x)), nrow=nrowM, ncol= ncolM)
	ee<-svd::trlan.svd(extM, neig=n)$d^2
	tr <- traceM
	tr2<-trace_mf(Mmult,tMmult,k=tr2.sample.size,ncolM,trace.full=traceM)
	tr.small<-tr-sum(ee)
	tr2.small<-tr2-sum(ee^2)
	scale<-tr2.small/tr.small
	nu<-(tr.small^2)/tr2.small
    pchisqsum(x, c(rep(1,n), ceiling(nu)), c(ee, scale), method=method,lower.tail=FALSE)
}
