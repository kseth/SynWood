library(testthat)
source("precFunctions.R")
x <- runif(10000, -10, 10)
y <- dnorm(x)

test_intervals <- 1:99/100 
out <- density_crI2(x, y, spar = 0.3, prI = test_intervals)
out2 <- density_crI(x, y)

cat("old prediction: ", out2$ll, " new prediction: ", out$cutoff_ll["0.95"], " true: ", dnorm(qnorm(0.5+0.95/2)), "\n")

preds <- c()

for(t in test_intervals)
{
	true <- dnorm(qnorm(0.5+t/2))
	pred <- out$cutoff_ll[as.character(t)]
	err <- abs(pred-true)/mean(pred, true)
	preds <- c(preds, err)
	cat(t, " ", err, "\n")
}

dev.new()
plot(test_intervals, preds, pch = "*", col = "grey", ylab = "% error", xlab = "precision interval")
abline(v = 0.95, col = "red")
## conclusions:
## as size sample (length(x)) + steps increases, % error goes to 0
## better to take the min of the ypairs than the mean of the ypairs

library(mvtnorm)
x1 <- seq(-5, 5, length.out=100)
x2 <- seq(-5, 5, length.out=100)
z <- mat.or.vec(length(x1), length(x2))
I <- matrix(rep(0, 4), ncol = 2)
diag(I) <- c(1, 1) 
for(i1 in 1:length(x1))
	for(i2 in 1:length(x2))
		z[i1, i2] <- dmvnorm(c(x1[i1], x2[i2]), c(0, 0), I)

dev.new()
start <- Sys.time()
out <- twoDim_precI(x1, x2, z, prI=c(0.5, 0.75, 0.95, 0.99), col=colorRampPalette(c("red", "orange", "yellow", "green"))(length(x1)))
# the below polygons show the 1D positions of the 99, 95, and 50% conf. intervals
# polygon(2.58/2*c(-2, -2, 2, 2), 2.58/2*c(-2, 2, 2, -2))
# polygon(1.96/2*c(-2, -2, 2, 2), 1.96/2*c(-2, 2, 2, -2))
# polygon(0.67/2*c(-2, -2, 2, 2), 0.67/2*c(-2, 2, 2, -2))
end <- Sys.time()
print(end - start)
