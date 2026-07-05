
# simule a arvore
tree <- phytools::pbtree(n = 100)

# distorça os ramos para forçar alta covariancia filogenética (não independencia)
tree_distorted <- compute.brlen(tree, method = "Grafen", power = 3)
tree_distorted$edge.length <- tree_distorted$edge.length / max(nodeHeights(tree_distorted))


random_fit <- 
  function(){
  x <- rnorm(n=100, sqrt(0.02))
  y <- rnorm(n=100, sqrt(0.02))
  fit <- lm(y~x)
  setNames(c(fit$coefficients[1], 
             fit$coefficients[2],
             anova(fit)[["Pr(>F)"]][1]),
           c("alfa", "beta", "p"))
}
X <- t(replicate(1000, random_fit()))

hist(X[,2], xlab= "beta (inclinação)", ylab = "Frequência", main = "beta")
hist(X[,3],xlab= "valor de p", ylab = "Frequência", main = "valor de p")
text(x = 0.8, y = 110, labels = "0.04", font = 2, cex = 1.2)

length(which(X[, 3]<= 0.05))/1000 # valor de p

brownian_lm <- 
  function(){
  x <- fastBM(tree_distorted)
  y <- fastBM(tree_distorted)
  fit <- lm(y~x)
  setNames(c(fit$coefficients[2],anova(fit)[["Pr(>F)"]][1]),c("beta","p"))
}

X_brownian <- t(replicate(1000, brownian_lm()))

hist(X_brownian[,1], xlab = "beta (inclinação)", ylab = "Frequência", main = "beta")
hist(X_brownian[,2], xlab="valor de p", ylab = "Frequência",  main="valor de p")
length(which(X_brownian[, 2] <= 0.05))/1000 # valor de p
text(x = 0.8, y = 110, labels = "0.84", font = 2, cex = 1.2)
plot(1, type = "n", xlim = c(0, 6), ylim = c(-4, 6),
     xlab = "X", ylab = "Y", main = "")



# Loop through the list and add lines
for (i in seq_along(vec_lm_1)) {
  abline(vec_lm_1[[i]], lwd = 1)
}
