    ##### This is for fitting 
    ##### lambda(t,x,y) = mu + K SUM gt(t-t_i)gxy(x-x_i,y-y_i) 
    ##### gt(t) = beta exp(-beta t). 
    ##### gxy(x,y) = f(r) = alpha/pi exp(-alpha r^2), with x^2+y^2=r^2,
    ##### The space S = [0,1] x [0,1] (km), in time [0,T]. 
    ##### In general the parameter vector theta = (mu, K, alpha, beta). 
    ##### However, to check it, we might just estimate one parameter 
    ##### and consider the others known. 

## Use simetas2022.r to simulate it. 
# source("simetasmay2017.r")
set.seed(43) 
T = 1000; mu = 8; theta_K = 0.3; theta_alpha = 2; theta_beta = 3
z = simhawk(T=T,gt = expgt, gxy = expxy, gmi = pointprod, mdensity = pointmag)

## Make sure the data are stored in z, and you define T and all the parameters externally. 

## Figure out which cell each point is in. Store a vector gridind and int m = number of cells. 
## For simplicity I will just use a 4x4 spatial grid of cells here. 
## So m = 16 and each cell has size T/m. 

m = 16
gridind = floor(z$lon*4)+floor(z$lat*4)*4 + 1

## Compile the C code to compute SG. 
system("R CMD SHLIB mysg.c")
dyn.load("mysg.so")

sgr1 = function(theta1,p){
    ## This tries one parameter at a time, keeping the others fixed at true values. 
    ## I've set it now on the pth parameter. 
    if(z$n > 100000) return(8e20) 
	## the memory allocation doesn't seem to work for really large n. 
    theta2 = c(8,.3,2,3) 
    theta2[p] = theta1 
    if(min(theta1) < 0.00000001) return(9e20)
    a3 = .C("sgc", as.double(z$lon), as.double(z$lat), as.double(z$t), as.integer(z$n), 
		as.double(T), as.double(theta2), as.integer(gridind), as.integer(m), double(1))
    points(theta1,a3[[9]],cex=.5,col=sample(20)[1]) 
    cat("\n ",theta1[1],", ",a3[[9]]) 
    a3[[9]]
}

par(mfrow=c(2,2)) 
### try it out
for(p in 1:4){
thetatrue = c(8,.3,2,3)
ymax = 100 
plot(c(0,2*thetatrue[p]),c(0,1000),type="n",xlab=paste("parameter ",as.character(p)),ylab="sg") 
abline(v = thetatrue[p]) 
for(i in 1:100) sgr1(i/100*2*thetatrue[p],1) 
} 

## theta0 = list(mu=8,K=.3,a=2,b=3)
## theta1 = as.vector(unlist(theta0))*1.5
## b1 = optim(theta0,sgr)
## theta1 = b1$par 
## sqrt(diag(solve(b2$hess))) ## for SEs 

########## For T = 1000, estimating one parameter at a time, when the truth was (8, .3, 2, 3), 
########## and starting values (12,.45,3,4.5), it gave me (7.71, 0.187, 0.790, 2.57). 
