include("../../gauher1.jl")

μ = 1
σ = 1

E(μ, σ) = exp(μ + (σ^2)/2)

g(z) = exp(z)

function gaussianHermite(g::Function, n)
	w, x = gauher(n)
	return π^(-1/2)*sum(summand.(g, x,w))
end

function summand(g::Function, x, w)
	return w*g(σ * x *sqrt(2) + μ)
end


print("Analytic Answer is:")
println(E(μ, σ))

println("numerical solution using Gaussian-Hermitian method")
println("for n=2 $(gaussianHermite(g,2))")
println("for n=10 $(gaussianHermite(g,10))")
