# Swiss Roll Example

library(foreach)
library(doParallel)
library(plgp)
library(splines)

source("R/functions.R")
ROLL <- get_swiss_roll(1000, 20000, 0.01)
y <- ROLL$y[1:500]
X <- ROLL$X[1:500, ]

ystar <- ROLL$y[501:600]
Xstar <- ROLL$X[501:600, ]

X <- scale(X, center = TRUE, scale = FALSE)
y <- scale(y, center = TRUE, scale = FALSE)

Xstar <- scale(Xstar, center = TRUE, scale = FALSE)
ystar <- scale(ystar, center = TRUE, scale = FALSE)

# Get screening scores
screening.scores <- get_NIS_scores(X, y)
par(mfrow =c(1,1))
plot(sort(screening.scores))

# Save the first n.keep covariates with the highest scores
n.keep = 10000
covars.save <- order(screening.scores, decreasing = TRUE)[1:n.keep]

# Get the screened X's
X.s <- X[,covars.save]
Xstar.s <- Xstar[,covars.save]

source("R/functions.R")
# Run the SkGP procedure
# K = 20 Models
# Sampling snrs and thetas on a grid
#   you can also set specific snrs using snr.method = "set"
out <- sketched_GP(y, X.s, ystar, Xstar.s, m = 60, K = 20,
                   SNRs = c(0.1, 0.5, 1), n.theta = 100, n.snrs = 25,prediction = TRUE,
                   snr.method = "set", snr.max = 1 ,stacking.method = "kfold", n.folds = 20)


stack.weights <- out$stack.weights
mu.pred <- t(out$mu.pred)
sig.pred <- t(out$sig.pred)
mu.fit <- t(out$mu.fit)
sig.fit <- t(out$sig.fit)

stack.mod.pred <- rowSums(mu.pred*matrix(data = unname(stack.weights), nrow = nrow(mu.pred), ncol = ncol(mu.pred), byrow = TRUE))
stack.sd.pred <-  sqrt(rowSums(sig.pred*matrix(data = unname(stack.weights), nrow = nrow(sig.pred), ncol = ncol(sig.pred), byrow = TRUE)))

stack.mod.fit <- rowSums(mu.fit*matrix(data = unname(stack.weights), nrow = nrow(mu.fit), ncol = ncol(mu.fit), byrow = TRUE))
stack.sd.fit <-  sqrt(rowSums(sig.fit*matrix(data = unname(stack.weights), nrow = nrow(sig.fit), ncol = ncol(sig.fit), byrow = TRUE)))


pmse.pred <- sum((stack.mod.pred - ystar)^2)/length(ystar)
ord = order(Xstar[,2])
plot(Xstar[ord,2], ystar[ord], main = paste0("Fit to Test Data: PMSE =  ", round(pmse.pred, 3)))
lines(Xstar[ord,2], stack.mod.pred[ord], col = "Maroon", lwd = 2)
lines(Xstar[ord,2], stack.mod.pred[ord] + stack.sd.pred[ord], col = "lightblue", lwd = 1)
lines(Xstar[ord,2], stack.mod.pred[ord] - stack.sd.pred[ord], col = "lightblue", lwd = 1)
