include("brent1.jl")
include("derivs.jl")
using Formatting, Printf
# ax < bx < cx, f(ax) > f(bx), and f(cx) > f(bx).


function f(x)
	return -(log(x) - x)
	# return x^2
end

println("\n"*"-"^80*"\n    Problem 2 - minimization \n"*"-"^80*"\n")


fmin,xmin,niter = brent(0,2,5, f,1.0e-8,1.0e-8)
fmax = -1*fmin
println("\n"*"-"^40*"\n    using brent's method: \n"*"-"^40*"\n")
# println("using brent's method: ")
# println("Function max obtained at $xmin with function value $fmax. converged after $niter iterations")
@printf("Function max obtained at %.4f with function value %s. converged after %d iterations\n\n", xmin, fmax, niter)

fmin,xmin,niter = newton_raphson(f,2;toler=1.0e-4,N=100)
fmax = -1*fmin
# println("using Newton-Raphson method: ")
println("\n"*"-"^40*"\n    using Newton-Raphson method:  \n"*"-"^40*"\n")

# println("Function max obtained at $xmin with function value $fmax. converged after $niter iterations")
@printf("Function max obtained at %.4f with function value %s. converged after %d iterations\n\n", xmin, fmax, niter)

