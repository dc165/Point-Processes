## First I took NSK.xls, saved it as text, 
## removed unnecessary columns, sorted by positive or negative, 
## and removed negatives. 

set.seed(8147) 
## load in data from mpox_cases_XY GPS_Sept13.xlsx 
setwd("C:/Users/14842/Desktop/Point Processes/mpox")
# x = read.table("NSK.txt",header=T) 
x = data_yahuma

## x is now 50 by 4. 
onset1 = as.Date(x[,1],"%m/%d/%y") 
sampdate1 = as.Date(x[,2],"%m/%d/%y") 
xcoord1 = x[,3] 
ycoord1 = x[,4] 

interv1 = as.numeric(sampdate1-onset1) 
mean(interv1,na.rm=T) ## 8.27 
hist(interv1,nclass=40) 

## If sampdate = NA, remove it. 
n = nrow(x)
ind = c(1:n)[!is.na(sampdate1)]  
onset2 = onset1[ind] 
sampdate2 = sampdate1[ind]
xcoord2 = xcoord1[ind] 
ycoord2 = ycoord1[ind] 
interv2 = interv1[ind] 
interv3 = interv2[!is.na(interv2)] 

n = length(ycoord2) 

## If onset = NA, sample from the distribution of sampdate - onset, to get 
## estimated onset times. 
for(i in 1:n){
if(is.na(onset2[i])){
  onset2[i] = sampdate2[i] - sample(interv3)[1] 
} 
} 

onset3 = julian(onset2) 
min(onset3) ## 19692. 

## Assume the first case is on day 1? It is unclear when surveillance started. 
## Assume the last day of surveillance is the last event + 1. 
## All the xs and ys are the same. 

onset4 = onset3 - (min(onset3)-1) 

z = list() 
z$t = sort(onset4) ## I sorted them here. 
z$lon = c(runif(z$n))
z$lat = c(runif(z$n))
z$n = length(z$t) 
T = max(z$t) + 1

system("R CMD SHLIB mysg.c")
dyn.load("mysg.dll")

X1 = 1
Y1 = 1
M0 = 3.5
m3 = function(x) signif(x,3) 

## Construct a list, wbin, where wbin[[17]] = c() if bin 17 is empty, and 
## if wbin[[17]] = c(1,2,10), then points 1,2,and 10 are in bin 17. 
## I will have 10 x 10 x 10 = 1000 bins. 
wbin = list()
for(i in 1:1000) wbin[[i]] = c(0)
for(m in 1:z$n) {
    gridindex = 10*10*floor(z$t[m]*10/T)+
    10*floor(z$lon[m]*10/X1)+ceiling(z$lat[m]*10/Y1)
    wbin[[gridindex]] = c(wbin[[gridindex]],m)
}
for(i in 1:1000) wbin[[i]] = wbin[[i]][-1] 

sumsqstoyan = function(theta,draw=0){
mu = theta[1]; K = theta[2]; tmean = theta[3]; tsd = theta[4]  
cat("\n mu = ",m3(mu),", K = ",m3(K),", tmean = ",m3(tmean),", tsd = ",m3(tsd),".\n") 
if(min(mu,K,tmean,tsd)<0.000000001) return(99999) 
if(K>.99999) return(99999)
if(draw){
r = seq(0,3,length=100)
t = alpha/pi * exp(-alpha * r^2)
lines(r,t,col="orange",lty=2) 
}
const = K
b = T*X1*Y1/10/10/10
mysum = rep(b,1000)
for(i in 1:1000){ ## i is the bin index. 
    if(length(wbin[[i]]) > .5){
       mysum[i] = 0
       for(j in wbin[[i]]){ ## j is the index of a point in bin i. 
           gkj = 0
           if(j>1) for(k in 1:(j-1)){ ## k are indices of previous points. 
              # r2 = (z$lon[j]-z$lon[k])^2+(z$lat[j]-z$lat[k])^2
              gkj = gkj + dnorm(z$t[j]-z$t[k],mean=tmean,sd=tsd)
           }
       lambdaj = mu/X1/Y1 + const*gkj
       if(lambdaj < 0){
	   cat("lambda ",j," is less than 0.")
	   return(99999)
       }
       mysum[i] = mysum[i] + 1/lambdaj
       }
    }
}
if(draw) lines(r,t,col="white",lty=2) 
sum((mysum-b)^2)
}

sumsqstoyan(c(1,.9,7,2))

theta1 = c(0.5,.9,7,2)
b1 = optim(theta1,sumsqstoyan)
b2 = optim(b1$par,sumsqstoyan,hessian=T)
theta2 = b2$par
sqrt(diag(solve(b2$hess))) ## for SEs 

## plot the fit. 
par(mfrow=c(1,1))
r = seq(0,20,length=100)
## s = theta0$beta * exp(-theta0$beta * r)
## plot(r,s,col="green",xlab="t",ylab="g(t)",type="l")
t = dnorm(r,mean=theta2[3],sd=theta2[4]) 
plot(r,t,col="blue",type="l",xlab="t (days)",ylab="estimated g(t)for NSK (cases/day)") 
## legend("topright",lty=c(1,1),c("real","estimated"),col=c("green","blue"))

expxy = function(n,m,theta){
    ## exponential triggering in space. f(r) = alpha/pi exp(-alpha r^2). 
    ## Here the density does not depend on magnitude of the mainshock. 
    ## To see that this is a density, 
    ## ∫f(x,y)dxdy = ∫f(r)rdrdø = 2π∫f(r)rdr 
    ## = 2alpha ∫ exp(-alpha r^2) r dr = -exp(-alpha r^2) ,r=0to∞, = 0+1, for alpha>0.  
    v = rexp(n,rate=theta$alpha)
    dist1 = sqrt(v)
    thet1 = runif(n)*2*pi
    x = cos(thet1)*dist1
    y = sin(thet1)*dist1
    cbind(x,y)
}

library(MASS) 
g4 = expxy(100000,1,theta=list(alpha=theta2[3]))
g5 = kde2d(g4[,1],g4[,2],lims=c(-.3,.3,-.3,.3))
image(g5,main="estimated planar triggering function g(x,y) for Equateur")
contour(g5,add=T) 
