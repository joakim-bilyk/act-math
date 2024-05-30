# Appendix B
## Numerical Studies
source("src/functions/reserveMC.R")

# Appendix C
## Simulations
### SDE
#### Weiner sampler
source("src/functions/simIncr.R")
source("src/functions/simBridge.R")
#### One-dim
source("src/functions/simSDE.R")
source("src/functions/Euler.R")
source("src/functions/Milstein.R")
source("src/functions/RungeKuttaStage.R")
source("src/misc/simSDEMisc.R")
#### Multi-dim
source("src/functions/Euler_multi.R")
source("src/functions/RK_multi.R")
### Pure-Jump
source("src/functions/simPureJump.R")
source("src/functions/RK_solve.R")
source("src/functions/RK.R")
source("src/functions/RK_epoch.R")
source("src/misc/simPureJumpMisc.R")
### Jump-Diffusion
source("src/functions/simJumpDiff.R")
source("src/misc/jumpDiffMisc.R")

# Appendix D
## Partial Differential Equations
### Two-Spacial Equations
source("src/functions/eulerOneDim.R")
source("src/functions/grid.R")
source("src/functions/eulerMatrixOneDim.R")
source("src/functions/eulerEpochOneDim.R")
source("src/functions/approxFun.R")
source("src/misc/PDEMisc.R")

### Three-Spacial Equations
source("src/functions/eulerTwoDim.R")
source("src/functions/eulerEpochTwoDim.R")
source("src/functions/eulerMatrixTwoDim.R")
source("src/functions/eulerPMatrix.R")
source("src/functions/approxFun3.R")

# Misc functions
source("src/functions/misc/bonds.R")
source("src/functions/misc/rungekutta.R")