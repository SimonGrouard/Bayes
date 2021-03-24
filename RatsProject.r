N <- 30
T <- 5
y <-structure(c(151, 145, 147, 155, 135, 159, 141, 159, 177, 134, 
              160, 143, 154, 171, 163, 160, 142, 156, 157, 152, 154, 139, 146, 
              157, 132, 160, 169, 157, 137, 153, 199, 199, 214, 200, 188, 210, 
              189, 201, 236, 182, 208, 188, 200, 221, 216, 207, 187, 203, 212, 
              203, 205, 190, 191, 211, 185, 207, 216, 205, 180, 200, 246, 249, 
              263, 237, 230, 252, 231, 248, 285, 220, 261, 220, 244, 270, 242, 
              248, 234, 243, 259, 246, 253, 225, 229, 250, 237, 257, 261, 248, 
              219, 244, 283, 293, 312, 272, 280, 298, 275, 297, 350, 260, 313, 
              273, 289, 326, 281, 288, 280, 283, 307, 286, 298, 267, 272, 285, 
              286, 303, 295, 289, 258, 286, 320, 354, 328, 297, 323, 331, 305, 
              338, 376, 296, 352, 314, 325, 358, 312, 324, 316, 317, 336, 321, 
              334, 302, 302, 323, 331, 345, 333, 316, 291, 324), .Dim = c(30, 
                                                                          5))
x <- c(8.0, 15.0, 22.0, 29.0, 36.0)
xbar <- 22
R <- structure(c(0.005, 0, 0, 5), .Dim = c(2, 2)) 
s <- 10^6
A <- t(structure(c(c(1,1,1,1,1),x), .Dim = c(5,2)))
B <- solve(t(A)%*%A)%*%t(A)
a <- b <- 0.001

mu.beta <- c(0,0)
tau <- 1
beta <- structure(c(100,100,100,100,100,
                    100,100,100,100,100,
                    100,100,100,100,100,
                    100,100,100,100,100,
                    100,100,100,100,100,
                    100,100,100,100,100,
                    6,6,6,6,6,
                    6,6,6,6,6,
                    6,6,6,6,6,
                    6,6,6,6,6,
                    6,6,6,6,6,
                    6,6,6,6,6),
                  .Dim = c(30, 2))
Omega <- structure(c(1,0,0,1),
                        .Dim=c(2,2))


Nchain <- 10000

chain.mu.beta <- matrix(NA, Nchain+1, 2)
chain.mu.beta[1,] <- mu.beta

chain.tau <- matrix(NA, Nchain+1)
chain.tau[1] <- tau

chain.beta <- matrix(NA, Nchain+1, 30*2)
chain.beta[1,] <- beta



chain.Omega <- matrix(NA, Nchain+1, 4)
chain.Omega[1,] <- Omega


for (iter in 1:Nchain){
  
  new.mu.beta <- chain.mu.beta[iter,]
  new.tau <- chain.tau[iter]
  new.beta <- structure(chain.beta[iter,],.Dim=c(30,2))
  new.Omega <- structure(chain.Omega[iter,],.Dim=c(2,2))
  
  # nouveau mu.beta
  
  mean.mu.beta <- solve(1/s*diag(1,2) + N*solve(new.Omega))%*%solve(new.Omega)%*%colSums(new.beta)
  cov.mu.beta <- solve(1/s*diag(1,2)+N*solve(new.Omega))
  
  new.mu.beta <- mvrnorm(1, mean.mu.beta, cov.mu.beta)                     
  
  # nouveau tau 
  mu <- matrix(NA, N, T)
  for (i in 1:N){
    mu[i,] <- t(A)%*%new.beta[i,]
  }
  
  tau.alpha <- a + N*T/2
  tau.beta <- b + sum((y-mu)^2)/2
  
  new.tau <- 1/rgamma(1, tau.alpha, tau.beta)
  
  #nouveau beta
  
  for (i in 1:N){
    sigma.beta <- solve(solve(new.Omega) + 1/new.tau*A%*%t(A))
    m.beta <- sigma.beta%*%(solve(new.Omega)%*%new.mu.beta + 1/new.tau*A%*%y[i,])
    
    new.beta[i,] <- mvrnorm(1, m.beta, sigma.beta)
  }
  
  
  # nouveau Omega
  
  ro <- 2 + 30 
  psi <- R 
  for (i in 1:N){
    psi <- psi + (new.beta[i,]-new.mu.beta)%*%t(new.beta[i,]-new.mu.beta)
  }
  
  new.Omega <- solve(structure(rWishart(1, ro, solve(psi)),.Dim=c(2,2)))
  
  
  chain.mu.beta[iter+1,] <- new.mu.beta
  chain.tau[iter+1] <- new.tau
  chain.beta[iter+1,] <- new.beta
  chain.Omega[iter+1,] <- new.Omega
}

#color = c('red','blue','green','purple','pink')
#plot(y[27,],type = 'o', lwd=3, ylim=c(100,400))
#for (i in 1:5){
#  lines(y[10+i,], type='o', lwd=3, col=color[i])
#}

#plot(seq(0,1000,0.1), dgamma(1/seq(0,1000,0.1),shape=0.001,rate=0.001))
