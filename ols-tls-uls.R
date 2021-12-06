rm(list=ls())
require('stable')
require('MASS')
alpha=1.5;beta=c(0.5,0.5)
n=50;p=0
BB=1
cmax=dmin=c()
mseteta11=mseteta12=mseteta21=mseteta22=mseteta31=mseteta32=teta11=teta12=teta21=teta22=teta31=teta32=c() 
mseteta11n=mseteta12n=mseteta21n=mseteta22n=mseteta31n=mseteta32n=teta11n=teta12n=teta21n=teta22n=teta31n=teta32n=c()
mseteta11a=mseteta12a=mseteta21a=mseteta22a=mseteta31a=mseteta32a=teta11a=teta12a=teta21a=teta22a=teta31a=teta32a=c() 

x1=runif(n,0,100)
x2=runif(n,100,300)
x=cbind(rep(1,n),x1,x2)
delta=c(0,0);
#dis=mvstable.isotropic(alpha,d=2,1,delta)
dis=mvstable.indep(alpha,beta,c(1,1),delta,param=1)
#dis=mvstable.elliptical(alpha, matrix(c(1,.5,.5,1),2,2), delta)
ee=rmvstable(dis,n)

e=matrix(ee,n,2,byrow=T)
teta=matrix(c(-5,2,10,7,.5,4),3,2,byrow=T)
y=x%*%teta+e
k11=lm(y~x-1)
k1=coefficients(k11)
rols=residuals(k11)
mseteta11=c(mseteta11,abs(k1[1,1]-teta[1,1]))
mseteta12=c(mseteta12,abs(k1[1,2]-teta[1,2]))
mseteta21=c(mseteta21,abs(k1[2,1]-teta[2,1]))
mseteta22=c(mseteta22,abs(k1[2,2]-teta[2,2]))
mseteta31=c(mseteta31,abs(k1[3,1]-teta[3,1]))
mseteta32=c(mseteta32,abs(k1[3,2]-teta[3,2]))

teta11=c(teta11,k1[1,1])
teta12=c(teta12,k1[1,2])
teta21=c(teta21,k1[2,1])
teta22=c(teta22,k1[2,2])
teta31=c(teta31,k1[3,1])
teta32=c(teta32,k1[3,2])
############################
############################
m=matrix(c(y,x,e),n)

mnew=matrix(c(y,x,e,(m[,6])+(m[,7])),n)
cutp=stable.fit.mle(mnew[,8], param=1)
if(round(cutp[2],4) < 1 & round(cutp[2],4)> -1){
cp1=ceiling(2/cutp[1]);dp1=floor(n+1-(2/cutp[1]))
}
if(round(cutp[1],4)>=1 & round(cutp[2],4)==1 |round(cutp[2],4)==-1){
cp1=ceiling(2/cutp[1]);dp1=floor(n+1-(2/cutp[1]))
}
if(round(cutp[1],4)<1 & round(cutp[2],4)==1 |round(cutp[2],4)==-1){
cp1=ceiling(2/cutp[1]);dp1=floor(n+1-(2/cutp[1]))
}
cp1;dp1
mm1=mnew[order(mnew[,8]),]
#mt1=data.frame(mm1[(max(c1,c2)+1):(min(d1,d2)-1),])
mt1=data.frame(mm1[(cp1):(dp1),])

yt=as.matrix(mt1[,1:2])
xt=as.matrix(mt1[,3:5])
kn=lm(yt~xt-1)$coef

mseteta11n=c(mseteta11n,abs(kn[1,1]-teta[1,1]))
mseteta12n=c(mseteta12n,abs(kn[1,2]-teta[1,2]))
mseteta21n=c(mseteta21n,abs(kn[2,1]-teta[2,1]))
mseteta22n=c(mseteta22n,abs(kn[2,2]-teta[2,2]))
mseteta31n=c(mseteta31n,abs(kn[3,1]-teta[3,1]))
mseteta32n=c(mseteta32n,abs(kn[3,2]-teta[3,2]))

teta11n=c(teta11n,kn[1,1])
teta12n=c(teta12n,kn[1,2])
teta21n=c(teta21n,kn[2,1])
teta22n=c(teta22n,kn[2,2])
teta31n=c(teta31n,kn[3,1])
teta32n=c(teta32n,kn[3,2])
#############################
################################

alpha;beta
mean(mseteta11);mean(mseteta12);mean(mseteta21);mean(mseteta22);mean(mseteta31);mean(mseteta32);mean(teta11);mean(teta12);mean(teta21);mean(teta22);mean(teta31);mean(teta32)
mean(mseteta11n);mean(mseteta12n);mean(mseteta21n);mean(mseteta22n);mean(mseteta31n);mean(mseteta32n);mean(teta11n);mean(teta12n);mean(teta21n);mean(teta22n);mean(teta31n);mean(teta32n)


#################################
c=cp1;d=dp1
for(j in c:d){
int3=function(ep1) {sapply(ep1, function(ep1) {
int2=function(ep2) {sapply(ep2, function(ep2) {
int1=function(kesi)

ep1*(dmvstable(dis,matrix(c(ep1,ep2),ncol=1))/dstable(kesi,cutp[1],cutp[2],cutp[3],cutp[4],1))*(j*choose(n,j)*(pstable(kesi,cutp[1],cutp[2],cutp[3],cutp[4],1)^(j-1))*((1-pstable(kesi,cutp[1],cutp[2],cutp[3],cutp[4],1))^(n-j))*dstable(kesi,cutp[1],cutp[2],cutp[3],cutp[4],1))

integrate(int1,0,4)$value
})
}
integrate(int2,-3,ep1)$value 
})
}
#integrate(int3,-2,2)#$value
k1=c(k1,integrate(int3,-3,3)$value)
}
#####################
####for e2
for(j in c:d){
int3=function(ep2) {sapply(ep2, function(ep2) {
int2=function(ep1) {sapply(ep1, function(ep1) {
int1=function(kesi)

ep2*(dmvstable(dis,matrix(c(ep1,ep2),ncol=1))/dstable(kesi,cutp[1],cutp[2],cutp[3],cutp[4],1))*(j*choose(n,j)*(pstable(kesi,cutp[1],cutp[2],cutp[3],cutp[4],1)^(j-1))*((1-pstable(kesi,cutp[1],cutp[2],cutp[3],cutp[4],1))^(n-j))*dstable(kesi,cutp[1],cutp[2],cutp[3],cutp[4],1))

integrate(int1,0,3)$value
})
}
integrate(int2,-3,3)$value 
})
}
#integrate(int3,-3,3)#$value
k11=c(k1,integrate(int3,-3,3)$value)
}

alpha;beta
################
mseteta11u=mseteta12u=mseteta21u=mseteta22u=mseteta31u=mseteta32u=teta11u=teta12u=teta21u=teta22u=teta31u=teta32u=c() 
B=1
n=50;
alpha=.7;beta=c(0,0)
for(i in 1:B){
x1=runif(n,0,100)
x2=runif(n,100,300)
x=cbind(rep(1,n),x1,x2)
delta=c(0,0);
#dis=mvstable.isotropic(alpha,d=2,1,delta)
dis=mvstable.indep(alpha,beta,c(1,1),delta,param=1)
#dis=mvstable.elliptical(alpha, matrix(c(1,.5,.5,1),2,2), delta)
#dis=mvstable.discrete.spec.meas2d(alpha,c(0,pi/2,pi,-pi/2),c(1/4,1/4,1/4,1/4),1,0,param=1,degrees=FALSE)
ee=rmvstable(dis,n)

e=matrix(ee,n,2,byrow=T)
teta=matrix(c(-5,2,10,7,.5,4),3,2,byrow=T)
y=x%*%teta+e
m=matrix(c(y,x,e),n)

mnew=matrix(c(y,x,e,(m[,6])+(m[,7])),n)

cp1=4;dp1=45
mm1=mnew[order(mnew[,8]),]
mt1=data.frame(mm1[(cp1+1):(dp1-1),])

yt=as.matrix(mt1[,1:2])
xt=as.matrix(mt1[,3:5])


vci=solve(vc)
k1u=solve(t(xt)%*%vci%*%xt)%*%t(xt)%*%vci%*%(yt-k1[5:44]-k11[5:44])
mseteta11u=c(mseteta11u,abs(k1u[1,1]-teta[1,1]))
mseteta12u=c(mseteta12u,abs(k1u[1,2]-teta[1,2]))
mseteta21u=c(mseteta21u,abs(k1u[2,1]-teta[2,1]))
mseteta22u=c(mseteta22u,abs(k1u[2,2]-teta[2,2]))
mseteta31u=c(mseteta31u,abs(k1u[3,1]-teta[3,1]))
mseteta32u=c(mseteta32u,abs(k1u[3,2]-teta[3,2]))

teta11u=c(teta11u,k1u[1,1])
teta12u=c(teta12u,k1u[1,2])
teta21u=c(teta21u,k1u[2,1])
teta22u=c(teta22u,k1u[2,2])
teta31u=c(teta31u,k1u[3,1])
teta32u=c(teta32u,k1u[3,2])
}

mean(mseteta11u);mean(mseteta12u);mean(mseteta21u);mean(mseteta22u);mean(mseteta31u);mean(mseteta32u);mean(teta11u);mean(teta12u);mean(teta21u);mean(teta22u);mean(teta31u);mean(teta32u)

