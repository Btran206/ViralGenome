N = 229354
n = 296
gene.ind <- read.table("hcmv.txt", header=TRUE)

par(mfrow=c(3,2))
stripchart(gene.ind$location, pch=21, cex=0.25, main = "hcmv")

set.seed(100)
rand1 <- sample.int(N, size=n, replace=FALSE)
stripchart(rand1, pch=21, cex=0.25, main = "Random 1")
set.seed(200)
rand2 <- sample.int(N, size=n, replace=FALSE)
stripchart(rand2, pch=21, cex=0.25, main = "Random 2")
set.seed(300)
rand3 <- sample.int(N, size=n, replace=FALSE)
stripchart(rand3, pch=21, cex=0.25, main = "Random 3")
set.seed(400)
rand4 <- sample.int(N, size=n, replace=FALSE)
stripchart(rand4, pch=21, cex=0.25, main = "Random 4")
set.seed(500)
rand5 <- sample.int(N, size=n, replace=FALSE)
stripchart(rand5, pch=21, cex=0.25, main = "Random 5")
#Q2
get.pairs <- function(sample.m){
  pair <- c()
  for (i in seq(1,length(sample.m))){
    x = sample.m[i] + sample.m[i+1]
    pair <- append(pair, x)
  }
  return (pair)
}
get.tripl <- function(sample.m){
  trip <-c()
  for (i in seq(1,length(sample.m))){
    x = sample.m[i] + sample.m[i+1] + sample.m[i+1]
    trip <- append(trip, x)
  }
  return (trip)
}
get.diff <- function(sample.pair){
  diff <- c()
  for(i in seq(1,length(sample.pair))){
    x = abs(sample.pair[i] - sample.pair[i+1])
    diff <- append(diff, x)
  }
  return (diff)
}

par(mfrow=c(1,1))

x = get.pairs(rand1)
g = get.tripl(rand1)
y = get.diff(x)
y.trip = get.diff(g)
plot(ecdf(y), lwd = 2, col = rgb(1,0,0,0.5), main = "ecdf")
set.seed(100)
theor <- rgamma(296, shape = 1, rate = 1/mean(y, na.rm = TRUE))
lines(ecdf(theor), lwd = 2, col = rgb(0,0,1,0.5))
legend(x = 0, y = 0.9, legend = c("Rand1 pair", "Gamma"), lty = c(1,1), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))
plot(ecdf(y.trip), lwd = 2, col = rgb(1,0,0,0.5), main = "ecdf")
set.seed(100)
theor <- rgamma(296, shape = 1, rate = 1/mean(y.trip, na.rm = TRUE))
lines(ecdf(theor), lwd = 2,  col = rgb(0,0,1,0.5))
legend(x = 0, y = 0.9, legend = c("Rand1 triple", "Gamma"), lty = c(1,1), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))

x = get.pairs(gene.ind$location)
g = get.tripl(gene.ind$location)
y = get.diff(x)
y.trip = get.diff(g)
plot(ecdf(y), lwd = 2, col = rgb(1,0,0,0.5), main = "ecdf")
set.seed(100)
theor <- rgamma(296, shape = 1, rate = 1/mean(y, na.rm = TRUE))
lines(ecdf(theor), lwd = 2,  col = rgb(0,0,1,0.5))
legend(x = 0, y = 0.9, legend = c("hcmv pair", "Gamma"), lty = c(1,1), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))
plot(ecdf(y.trip), lwd = 2, col = rgb(1,0,0,0.5), main = "ecdf")
set.seed((100))
theor <- rgamma(296, shape = 1, rate = 1/mean(y.trip, na.rm = TRUE))
lines(ecdf(theor), lwd = 2,  col = rgb(0,0,1,0.5))
legend(x = 0, y = 0.9, legend = c("hcmv triple", "Gamma"), lty = c(1,1), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))
### end of q2
regionsplit <- function(n.region, gene, site){
  count.int <- table(cut(site, breaks = seq(1, length(gene), length.out=n.region+1), include.lowest=TRUE))
  count.vector <- as.vector(count.int)
  count.tab <- table(factor(count.vector,levels=0:max(count.vector)))
  return (count.tab)
}
n.region <- 57
gene <- seq(1,N)
rand1.tab=regionsplit(n.region, gene, rand1)
rand2.tab=regionsplit(n.region, gene, rand2)
rand3.tab=regionsplit(n.region, gene, rand3)
rand4.tab=regionsplit(n.region, gene, rand4)
rand5.tab=regionsplit(n.region, gene, rand5)
real.tab = regionsplit(n.region, gene, gene.ind$location)
trunc=7
lvls=factor(c(0:(trunc-1),paste(">=",trunc,sep="")),levels=c(0:(trunc-1),paste(">=",trunc,sep="")))
O1=as.vector(rand1.tab)
O2=as.vector(rand2.tab)
O3=as.vector(rand3.tab)
O4=as.vector(rand4.tab)
O5=as.vector(rand5.tab)
real=as.vector(real.tab)
O1.trunc=c(O1[1:trunc],sum(O1[-(1:trunc)]))
O2.trunc=c(O2[1:trunc],sum(O2[-(1:trunc)]))
O3.trunc=c(O3[1:trunc],sum(O3[-(1:trunc)]))
O4.trunc=c(O1[1:trunc],sum(O4[-(1:trunc)]))
O5.trunc=c(O1[1:trunc],sum(O5[-(1:trunc)]))
R4.trunc=c(O1[1:trunc],sum(real[-(1:trunc)]))
lambda=n/n.region
p=c(dpois(0:(trunc-1),lambda),1-sum(dpois(0:(trunc-1),lambda)))
E=p*n.region
chisq.test(O1.trunc,p=p,simulate.p.value=TRUE)
chisq.test(O2.trunc,p=p,simulate.p.value=TRUE)
chisq.test(O3.trunc,p=p,simulate.p.value=TRUE)
chisq.test(O4.trunc,p=p,simulate.p.value=TRUE)
chisq.test(O5.trunc,p=p,simulate.p.value=TRUE)
chisq.test(R4.trunc,p=p,simulate.p.value=TRUE)

par(mfrow=c(2,3))
plot((R4.trunc - E) / sqrt(E), type = 'h', ylab = "standardized residuals for hcmv", xlab = "interval index")

plot((O1.trunc - E) / sqrt(E), type = 'h', ylab = "standardized residuals for Rand1", xlab = "interval index")
plot((O2.trunc - E) / sqrt(E), type = 'h', ylab = "standardized residuals for Rand2", xlab = "interval index")
plot((O3.trunc - E) / sqrt(E), type = 'h', ylab = "standardized residuals for Rand3", xlab = "interval index")
plot((O4.trunc - E) / sqrt(E), type = 'h', ylab = "standardized residuals for Rand4", xlab = "interval index")
plot((O5.trunc - E) / sqrt(E), type = 'h', ylab = "standardized residuals for Rand5", xlab = "interval index")


nsim = 2000   # number of simulations

# Simulate data
max.count = rep(NULL, nsim)
for(i in 1:nsim) {
  x = sample(1:N, size = n, replace = F)
  y = regionsplit(n.region,gene, x)
  y = max(y)
  max.count <- append(max.count,y)
}
max.count
real.tab
hist(max.count, main = "Distribution of max counts for randomly generated samples")
lines(x = max(real) + .5,y = 0, type = "p", col = rgb(1,0,0,0.5), pch = 19, cex = 1.5)
legend(35, 400, legend = "Maximum Count in hcmv",
       pch = 19,
       col = rgb(1,0,0,0.5), text.col = 1)

