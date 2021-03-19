include("brent1.jl")
include("derivs.jl")
using LinearAlgebra, Printf
# import LinearAlgebra: norm, inv
# ax < bx < cx, f(ax) > f(bx), and f(cx) > f(bx).

# a_A = a_B = 0.1

# b_A = b_B = 1
# e_A = e_B = 1

function supply(p, a)
    return exp.(a.*p) .- 1.0
end

# demand function
function demand(p, b, ϵ)
    pbar = sum(p) / length(p)
    return b.*(p./pbar).^(-ϵ)
end

z(p, a, b, ϵ) = demand(p, b, ϵ) .- supply(p, a) # excess demand function


objective(p, a, b, ϵ) =  sum(z(p, a, b, ϵ).^2)


# a = objective(p, 0.1, 1.0, 1.0)

# fmin,xmin,niter = brent(ax=0,bx=2,cx=5,f,1.0e-8,1.0e-8)

println("\n"*"-"^80*"\n    Problem 3 - minimization \n"*"-"^80*"\n")

fmin,xmin,niter = newton_raphson(p -> objective(p, 0.1, 1.0, 1.0), [3.0; 3.0])

# println("Function min obtained at $xmin with minimized value $fmin. converged after $niter iterations")
println("\n"*"-"^40*"\n    using Newton-Raphson method:  \n"*"-"^40*"\n")

@printf("Function min obtained at %s with function value %s. converged after %d iterations\n\n", xmin, fmin, niter)
