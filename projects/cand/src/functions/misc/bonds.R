# bonds.R
bondPrice <- function(
  t,T,r,
  par, # The parameters for the rate model in alphabetic order
  # for instance c(a,b,sigma) for the Vasicek model
  model = "Vasicek"
  # Vasicek:  dr = b(a-r)dt + sigma dW
  # CIR:      dr = b(a-r)dt + sigma sqrt(r) dW
) {
  const <- switch (model,
    Vasicek = c(1),
    CIR = sqrt(par[2]**2+2*par[3]**2) # gamma = (b^2 + 2*sigma^2)^(1/2)
  )
  A0 <- switch (model,
    Vasicek = function(x) {1},
    CIR = function(x) {
      (
        (2*const*exp((par[2]+const)*x/2))
          / ((const + par[2])*(exp(const*x)-1)+2*const)
      )**(2*par[1]*par[2]**2/(par[3]**2))
    }
  )
  B <- switch (model,
    Vasicek = function(x) {
      (1/par[2])*(1-exp(-par[2]*x))
    },
    CIR = function(x) {
      (2*(exp(const*x)-1)) / ((const + par[2])*(exp(const*x)-1)+2*const)
    }
  )
  A <- switch (model,
               Vasicek = function(x) {
                 (B(x)-x)*(par[1]-0.5*par[3]**2*par[2]**(-2)) -
                   B(x)**2*par[3]**2/(4*par[2])
               },
               CIR = function(x) {0}
  )
  A0(T-t)*exp(A(T-t)-B(T-t)*r)
}