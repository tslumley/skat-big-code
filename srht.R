if (!is.loaded("mfwht")) dyn.load("fht.so")
##subsampled randomised hadamard transform

fht<-function(A){
	m<-NROW(A)
	mbig<-2^ceiling(log2(m))
	Abig<-matrix(0,nrow=mbig,ncol=NCOL(A))
	Abig[1:m,]<-A
	.C("mfwht", Abig, as.integer(mbig), as.integer(NCOL(Abig)))[[1]]
}


##
## square matrices: find an orthonormal basis for the range of a rank-k approximation to A
##
srfht<-function(A,k,q=0){
	m<-NROW(A)
	mbig<-2^ceiling(log2(m))
	R<-sample(c(-1,1),m,replace=TRUE)
	Astar<-A*R
	idx<-sample(mbig,k)
	AOmega<-fht(Astar)[idx,]/sqrt(k)
	Q<-qr.Q(qr(t(AOmega)))
	for(i in seq_len(q)){
		tildeQ<-qr.Q(qr(t(A)%*%Q))
		Q<-qr.Q(qr(A%*%tildeQ))
	}
	Q
	}



##
##  x: the statistic value
##  M: the matrix whose eigenthingies we need
##  n: how many largest eigenvalues to use (time roughly proportion to this)
##  p: oversampling factor (should be about log(NCOL(M)), maybe?)
##  method: what method to use for computing the distribution based on the exact eigenvalues

pchisqsum_ssvd<-function(x,M,n=100,p=10,q=0, method=c("saddlepoint","integration")){
	method<-match.arg(method)
	Q<-srfht(M,n+p,q=q)
	B<-t(Q)%*%M%*%Q
	ee<-eigen(B,symmetric=TRUE,only.values=TRUE)$values[1:n]
	tr<-sum(diag(M))
	tr2<-sum(M^2)
	tr.small<-tr-sum(ee)
	tr2.small<-tr2-sum(ee^2)
	scale<-tr2.small/tr.small
	nu<-(tr.small^2)/tr2.small
    pchisqsum(x, c(rep(1,n),ceiling(nu)), c(ee, scale), method=method,lower.tail=FALSE)
}	

##
## non-square matrices
##
srfht2<-function(A,k,q=0){
	m<-NROW(A)
	mbig<-2^ceiling(log2(m))
	R<-sample(c(-1,1),m,replace=TRUE)
	Astar<-A*R
	idx<-sample(mbig,k)
	AOmega<-fht(Astar)[idx,]/sqrt(k)
	Q<-qr.Q(qr(t(AOmega)))
	for(i in seq_len(q)){
		tildeQ<-qr.Q(qr(A%*%Q))
		Q<-qr.Q(qr(t(A)%*%tildeQ))
	}
	Q
	}


	
##
## Using a Hutchinson-style randomized trace estimator. 
##
pchisqsum_rsvd<-function(x,M,n=100,p=10,q=2, tr2.sample.size=100, method=c("saddlepoint","integration")){
	method<-match.arg(method)
	Q<-srfht2(M,n+p,q=q)
	B<-M%*%Q
	ee<-svd(B,nu=0,nv=0)$d[1:n]^2
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




## hutchinson-like trace estimator but with FHT. 
## improved with ratio estimator: about 50% std error reduction at k=100.
tracefht<-function(A,k,trace.full=NULL){
	m<-NROW(A)
	mbig<-2^ceiling(log2(m))
	R<-sample(c(-1,1),m,replace=TRUE)
	Astar<-A*R
	idx<-sample(mbig,k)
	AOmega<-fht(Astar)[idx,]
	AAOmega<-tcrossprod(AOmega,A)
	if (!is.null(trace.full))
		tr<-sum(rowSums(AOmega*AOmega))/k
	trsquared<-sum(rowSums(AAOmega*AAOmega))/k
	if (is.null(trace.full))
		trsquared
	else
		trsquared*(trace.full/tr)^2 
}

pchisqsum_Gsatt<-function (x, M, tr2.sample.size=100){
	  tr<-sum(M^2)
	  tr2<-tracefht(M,k=tr2.sample.size,trace.full=tr)
	  scale<-tr2/tr
	  nu<-(tr^2)/tr2
      pchisq(x/scale, nu,lower.tail=FALSE)
	}
##
## liu method 
##

pchisqsum_liu<-function (q, M, tr2.sample.size=300) 
{
	c1<-sum(M^2)
	traces<-tracefht_liu(M,k=tr2.sample.size,trace.full=c1)
	c2 <- traces[1]
    c3 <- traces[2]
    c4 <- traces[3]
    s1 <- c3/(c2^(3/2))
    s2 <- c4/c2^2
    muQ <- c1
    sigmaQ <- sqrt(2 * c2)
    tstar <- (q - muQ)/sigmaQ
    if (s1^2 > s2) {
        a <- 1/(s1 - sqrt(s1^2 - s2))
        delta <- s1 * a^3 - a^2
        l <- a^2 - 2 * delta
    }
    else {
        a <- 1/s1
        delta <- 0
        l <- c2^3/c3^2
    }
    muX <- l + delta
    sigmaX <- sqrt(2) * a
    Qq <- pchisq(tstar * sigmaX + muX, df = l, ncp = delta, lower.tail = FALSE)
    return(Qq)
}

## slower than using eigenthingies
pchisqsum_liufull<-function (q, H) 
{
	c1<-sum(diag(H))
	H2<-crossprod(H)
	c2 <- sum(diag(H2))
	H3<-crossprod(H,H2)
    c3 <- sum(diag(H3))
    c4 <- sum(H3^2)
    s1 <- c3/(c2^(3/2))
    s2 <- c4/c2^2
    muQ <- c1
    sigmaQ <- sqrt(2 * c2)
    tstar <- (q - muQ)/sigmaQ
    if (s1^2 > s2) {
        a <- 1/(s1 - sqrt(s1^2 - s2))
        delta <- s1 * a^3 - a^2
        l <- a^2 - 2 * delta
    }
    else {
        a <- 1/s1
        delta <- 0
        l <- c2^3/c3^2
    }
    muX <- l + delta
    sigmaX <- sqrt(2) * a
    Qq <- pchisq(tstar * sigmaX + muX, df = l, ncp = delta, lower.tail = FALSE)
    return(Qq)
}



tracefht_liu<-function(A,k,trace.full=sum(A^2)){
	m<-NROW(A)
	mbig<-2^ceiling(log2(m))
	R<-sample(c(-1,1),m,replace=TRUE)
	Astar<-A*R
	idx<-sample(mbig,k)
	AOmega<-fht(Astar)[idx,]
	AAOmega<-tcrossprod(AOmega,A)
	A3Omega<-AAOmega%*%A
	A4Omega<-tcrossprod(A3Omega,A)
	
	tr<-sum(rowSums(AOmega*AOmega))/k
	trAA2<-sum(rowSums(AAOmega*AAOmega))/k
	trAA3<-sum(rowSums(A3Omega*A3Omega))/k
	trAA4<-sum(rowSums(A4Omega*A4Omega))/k

	c(trAA2*(trace.full/tr)^2, trAA3*(trace.full/tr)^3,trAA4*(trace.full/tr)^4)
}


##
## Matrix-free for sparse G?
##

matfreeQ<-function(Mmult,tMmult,m,cols,q=0){
	Omega<-matrix(rnorm(m*cols), ncol=m)
	AOmega<-Mmult(Omega)
	Q<-qr.Q(qr(AOmega))
	for(i in seq_len(q)){
		tildeQ<-qr.Q(qr(tMmult(Q)))
		Q<-qr.Q(qr(Mmult(tildeQ)))
	}
	Q
	}

trace_mf<-function(Mmult,tMmult,k,cols,trace.full=NULL){
    Omega<-matrix(rnorm(k*cols), ncol=k)
	AOmega<-Mmult(Omega)
	AAOmega<-tMmult(AOmega)
	if (!is.null(trace.full))
		tr<-sum(rowSums(AOmega*AOmega))/k
	trsquared<-sum(rowSums(AAOmega*AAOmega))/k
	if (is.null(trace.full))
		trsquared
	else
		trsquared*(trace.full/tr)^2 
}


pchisqsum_smf<-function(x,Mmult,tMmult,ncolM,traceM, n=100,p=10,q=2, tr2.sample.size=300, method=c("saddlepoint","integration")){
	method<-match.arg(method)
	Q<-matfreeQ(Mmult,tMmult, n+p,ncolM,q=q)
	B<-tMmult(Q)
	ee<-svd(B,nu=0,nv=0)$d[1:n]^2
	tr <- traceM
	tr2<-trace_mf(Mmult,tMmult,k=tr2.sample.size,ncolM,trace.full=traceM)
	tr.small<-tr-sum(ee)
	tr2.small<-tr2-sum(ee^2)
	scale<-tr2.small/tr.small
	nu<-(tr.small^2)/tr2.small
    pchisqsum(x, c(rep(1,n), ceiling(nu)), c(ee, scale), method=method,lower.tail=FALSE)
}

