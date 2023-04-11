rm(list= ls())

p <- 3
lambda <- matrix(
  1:9,
  nrow = p
)
lambda
k <- 3
helper <- matrix(rep(0,k**2),ncol = p)
big_matrix <- Matrix::bdiag(lapply(1:k, function(x) lambda))


### 2. Marts - øvelse
M <- t(matrix(
  c(-0.7,0.1,0.1,0.0,0.5,
    0,-0.5,0,0,0.5,
    0,0.1,-0.7,0.1,0.5,
    0,0.1,0,-0.6,0.5,
    0,0,0,0,0), ncol = 5,nrow = 5
))
M
rho <- -1
rate <- c(rho, 1, 1, rho,0)
lump_sum <- t(
  matrix(
    c(0,2,0,0,0,
      0,0,0,0,0,
      0,2,0,0,0,
      0,2,0,0,0,
      0,0,0,0,0), ncol = 5
  )
)
lump_sum
delta_rate <- diag(rate)
delta_rate
lambda_1 <- t(
  matrix(
    c(0,0.05,0,0,0,
      0,0,0,0,0,
      0,0.1,0,0,0,
      0,0.05,0,0,0,
      0,0,0,0,0), ncol = 5
  )
)
lambda_1
lamda_0 <- M - lambda_1
r <- 0.04
R <- lambda_1 * lump_sum + delta_rate

big_matrix <- as.matrix(Matrix::bdiag(M-r*diag(rep(1,5),5),M))
big_matrix[1:5,6:10] <- R
T <- 1
van_loan <- Matrix::expm(big_matrix*T)
V <- van_loan[1:5,6:10] %*% rep(1,5)
V
## Equivalence premium

R_derivative <- t(matrix(
  c(1,0,0,0,0,
    0,0,0,0,0,
    0,0,0,0,0,
    0,0,0,1,0,
    0,0,0,0,0)
))
big_matrix_top <- as.matrix(Matrix::bdiag(M-r*diag(rep(1,5),5),M))
big_matrix_top[1:5,6:10] <- R
top <- -c(1,rep(0,4)) %*% Matrix::expm(big_matrix_top*T)[1:5,6:10] %*% rep(1,5)
big_matrix_buttom <- as.matrix(Matrix::bdiag(M-r*diag(rep(1,5),5),M))
big_matrix_buttom[1:5,6:10] <- R_derivative
buttom <- c(1,rep(0,4)) %*% Matrix::expm(big_matrix_buttom*T)[1:5,6:10] %*% rep(1,5)

rho <- as.numeric(top/buttom)

m <- list(M,R)

bToeplitz <- function(m) {
  n <- length(m) #number of diagonals
  p <- unique(as.numeric(sapply(1:n, function(i) dim(m[[i]])))) #dimensions
  if (length(p) > 1) {
    print("Aborting: Dimensions is not equal!")
  } else {
    require(Matrix)
    result <- matrix(rep(0,p**2*n**2),ncol = n*p)
    for (i in 1:n) {
      m_diag <- m[[i]]
      increment <- as.matrix(Matrix::bdiag(
        lapply(1:(n-i+1), function(j){
          m_diag
        })
      ))
      result[1:((n-i+1)*p),(p*(i-1)+1):(n*p)] <- result[1:((n-i+1)*p),(p*(i-1)+1):(n*p)]+ increment
    }
    return(result)
  }
}
bToeplitz(list(M,R))
big_matrix <- bToeplitz(
  list(M,
       lambda_1*lump_sum + diag(rate),
       (1/2)*lambda_1*lump_sum*lump_sum)
) - as.matrix(Matrix::bdiag(2*r*diag(1,5),
                            r*diag(1,5),
                            0*diag(1,5)))

big_matrix <- big_matrix - as.matrix(Matrix::bdiag(diag(1,5)))
tmp <- expm::expm(big_matrix)
2*tmp[1:5,11:15] %*% rep(1,5)
tmp[6:10,11:15] %*% rep(1,5)
