# Find largest subcircle and continue

## Libraries

library(DEoptim)

## Graphical utilities

plotCircle <- function(x, y, r, n = 100, ...) {
    # C = (x, y, r)
    z  <- seq(0, 2*pi, length.out = n)
    xs <- x + r * cos(z)
    ys <- y + r * sin(z)
    # border color: border = NULL; fill color: col = NA
    polygon(xs, ys, ...)
    invisible(NULL)
}

## Random points

Pts = read.table("RandomPoints.dat", header =TRUE)
n = nrow(Pts)

## Plot unit circle

opar = par(mar=c(1,1,1,1))                      # small plot outlet
plot(c(-1,1), c(-1,1), type='n',asp = 1)        # plot geometry
plotCircle(0,0,1, col="gray95")                 # plot unit circle
for (i in 1:n) {                                # plot points in color
    if (Pts[i,3] == 0){
        points(Pts[i,1], Pts[i, 2], col = "red", pch=20)
    } else {
        plotCircle(Pts[i,1], Pts[i,2], Pts[i,3], border = FALSE, col = "red")
    }
}

## Generate K circles /w max radius

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
    
    #-- Call the global optimization solver
    # library(DEoptim)
    lb = c(-1, -1); ub = c(1, 1)
    sol = DEoptim::DEoptim(objfun, lb, ub,
                   DEoptim::DEoptim.control(trace=FALSE, itermax=500))
    s = c(unname(sol$optim$bestmem[1]), unname(sol$optim$bestmem[2]), 
          -sol$optim$bestval)
    # Refine with local optimizer
    # sol = optim(par=c(s[1], s[2]), fn=objfun, method="BFGS")
    # s = c(sol$par[1], sol$par[2], -sol$value)
    q
    #-- Plot the resulting circle 
    plotCircle(s[1], s[2], s[3], 
               col=k, border=k)  # border = "white"
    text(s[1], s[2], k, col = "white", cex = 0.75)
    
    #-- Add the circle to the existing data
    Pts = rbind(Pts, c(s[1], s[2], s[3]))
}

## Plot the sequence of radii

plot(Pts[31:(30+24), 3], type="h", main="Sequence of 'minimal' radii")
lines(Pts[31:(30+24), 3], col="blue", lwd=2)
grid()

## Clean up

par(opar)

