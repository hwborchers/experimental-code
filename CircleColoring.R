## Find largest subcircle and continue

## Libraries
# library(DEoptim)
  library(DEoptimR)
# library(pso)
# library(GenSA)    # too slows, stops with error

#-- Graphical Utilities
plotCircle <- function(x, y, r, n = 100, ...) {
    # C = (x, y, r)
    z  <- seq(0, 2*pi, length.out = n)
    xs <- x + r * cos(z)
    ys <- y + r * sin(z)
    # border color: border = NULL; fill color: col = NA
    polygon(xs, ys, ...)
    invisible(NULL)
}

#-- Generate random points
# n = 40
# x = 2*runif(n)-1.0
# y = 2*runif(n)-1.0
# r = numeric(n)
# o = which(x^2+y^2<=1)
# x = x[o]; y = y[o]

#-- Read random points and set radius to 0.05
Pts = read.table("RandomPoints.dat", header =TRUE)
n = nrow(Pts)
# Pts[, 3] = 0.05

#-- Set plot options and plot unit circle and points
par(mar = c(1,1,1,1))                           # small plot outlet
plot(c(-1,1), c(-1,1), type='n',asp = 1)        # plot geometry
plotCircle(0,0,1, col="gray95")                 # plot unit circle
for (i in 1:n) {                                # plot points in color
    plotCircle(Pts[i,1], Pts[i, 2], Pts[i, 3],
               border = FALSE, col = "red")
}
points(Pts$x, Pts$y, pch=20, col="red")  # , cex=2.0

pracma::tic()
K = 48
for (k in 1:K) {
#-- Objective function
objfun = function(p) {
    xp=p[1]; yp=p[2]
    if (xp^2 + yp^2 >= 1) return(0.0)
    d0 = 1.0 - sqrt(xp^2 + yp^2)
    for (i in 1:nrow(Pts)) {
        di = sqrt((Pts[i,1]-xp)^2 + (Pts[i,2]-yp)^2) - Pts[i,3]
        if (di <= 0.0) {
            return(0.0)
        } else {
            if (di < d0) {
                d0 = di
            }
        }
    }
    return(-d0)
}

lb = c(-1, -1); ub = c(1, 1)

#-- Call the global optimization solver
# library(DEoptim)
# sol = DEoptim(objfun, lb, ub,
#               DEoptim.control(trace=FALSE, itermax=500))
# s = c(unname(sol$optim$bestmem[1]), unname(sol$optim$bestmem[2]), 
#              -sol$optim$bestval)
# sol = optim(par=c(s[1], s[2]), fn=objfun, method="BFGS")
# s = c(sol$par[1], sol$par[2], -sol$value)

# library(pso)
# sol = psoptim(par=c(0,0), fn=objfun, lower=-1, upper=1)
# s = c(sol$par[1], sol$par[2], -sol$value)

# library(DEoptimR)
  sol = JDEoptim(lb, ub, objfun)
  s = c(sol$par[1], sol$par[2], -sol$value)

# library(GenSA)
# sol = GenSA(par=NULL, objfun, lower=lb, upper=ub)
# s = c(sol$par[1], sol$par[2], -sol$value)

# library(Rmalschains)
# sol = malschains(objfun, lower=lb, upper=ub,
#                  maxEvals=400, verbosity=0,
#         control=malschains.control(popsize=100, ls="sw", istep=100))
# s = c(sol$sol[1], sol$sol[2], -sol$fitness)

# library(ABCoptim)
# sol = abc_optim(runif(2), objfun, lb=lb, ub=ub)
# s = c(sol$par[1], sol$par[2], -sol$value)

#-- Plot the resulting circle 
plotCircle(s[1], s[2], s[3], 
           col=k, border=k)  # "white"
text(s[1], s[2], k, col = "white")

#-- Add the circle to the existing data
Pts = rbind(Pts, c(s[1], s[2], s[3]))
}

points(Pts$x[1:n], Pts$y[1:n], pch=20, col="red")
for (i in 1:K) {
    text(Pts[i+n,1], Pts[i+n,2], i, col="white", cex=0.75)
}
pracma::toc()

## Plot the series of radii
par(mar=c(2,2,2,1))
plot(Pts[31:(30+K), 3], type="h", main="Sequence of 'minimal' radii")
lines(Pts[31:(30+K), 3], col="blue", lwd=2)
grid()

r = Pts[31:(30+K), 3]
rs = sort(r,decreasing = TRUE)
lines(rs, col = 2, lty = 2, lwd=2)

cat("coverage =",sum(r^2), '\n')
dev = pracma::trapz(1:K, abs(r-rs))
cat("deviance =", dev, '\n')

## Clean up

par(opar)

