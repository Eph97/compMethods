include("../../brent1.jl")
include("interpolate.jl")


β = 0.9; α=0.4; δ=0.1; p = 0.9
k_bar = ((1/β-(1-δ))/(α))^(1/(α-1)) # analytic steady state (solution to F'(k)=1/β)
# include("derivs.jl")
# using Plots
using PyPlot
using LaTeXStrings

# returns grid of N equally-spaced points on interval [xmin, xmax]
function createGrid(xmin, xmax, N)
	step = (xmax-xmin) / (N-1)
	return collect(xmin:step:xmax)
end

function bellman(U, k_i, k_j,spline)
	return U(F(k_i) - k_j) + β*interp(k_j,spline)[1]
end

# returns (discrete) value and policy functions on grid defined by the vector a
function valueFuncIter(U, k; β=0.9 , toler=1e-6, maxiter=1000, verbose=1)
	N = length(k)
	v0=zeros(length(k))
	g = fill(0.0, N) # g[i] will be j<=i that maximizes U(a_i-a_j) + βv_j (this is the optimal decision rule/policy function)
	# U_ar = [U(F(k[i])-k[j]) for i in 1:N, j in 1:N] # NxN array of all possible values of U, given grid a
	err = 1+toler; niter = 0;
	tol = zeps = 1.0e-8
	while err > toler && niter < maxiter
		niter += 1
		v1 = copy(v0)

		v0Spline = makespline(k,v0)

		# define new value function by iterating over all grid points
		for i in 1:N
			fmin,xmin,niter = brent(0, ((F(k[i]))/2), F(k[i]), k_j -> (-1*bellman(U, k[i], k_j,v0Spline)), tol, zeps)
			v1[i] = -1*fmin
			g[i] = xmin
		end

		# rhs = U_ar[i,:] .+ β*v0[:] # RHS of Bellman equation (for all a_j<=a_i)
		# v1[i], g[i] = findmax(rhs) # simultaneously update value function and get the optimizing j (for given i)
		# (-1*rhs), 

		err = maximum(abs.(v1 -v0)) # error defined using the sup norm
		v0 = v1
	end

	if err <= toler && verbose >= 1
		println("Value function iteration converged in $niter iterations (maxerr = $err)")
	elseif err > toler
		println("Value function iteration failed to converge in $niter iterations (maxerr = $err)")
	end

	println("number of iterations", niter)
	return v0, g
end

# returns the indices (in a) corresponding to the optimal path of savings, when a0 = a[idx0]
function getPath(g, idx0, T)
	idxs = fill(idx0, T) # optimal path of cake stock (expressed as an index of a)
	for t in 2:T
		idxs[t] = g[idxs[t-1]]
	end
	return idxs
end

function f(k; alpha=0.4)
	return k^alpha
end

function F(k; alpha=0.4, delta=0.1)
	return f(k) + (1 - delta)*k
end


# the consumer's felicity function
function U(c)
	# we don't allow negative consumption, even if c > -ϵ
	if c < 0
		return -Inf
	else
		return log(c)
	end
end

N = 15 # number of grid points
a = 0.5
b = 1.5
# k = createGrid(a*k_bar, b*k_bar, N)
#
k = createGrid(a*k_bar,b*k_bar, N)
vsol, g = valueFuncIter(U, k)


function getContinuousPath(spline, initial::Float64, T)
	# klo, khi = findind(initial, spline)
	path = fill(initial, T) # optimal path of cake stock (expressed as an index of k)
	for t in 2:T
		path[t] = interp(path[t-1],spline)[1]
	end
	return path
end

yspline = makespline(k,g)


fig1 = figure()
subplot(121)
plot(k, vsol) 
title("Value Function ( N= $N)") 
xlabel(L"k")
ylabel(L"v(k)")

subplot(122)
plot(k, g, label="decision rule") 
plot(k, k, label="45 deg line") 
title("Optimal Decision Rule") 
xlabel(L"k")
ylabel(L"k'")
legend()
fig1.savefig("fig1.png")


T = 41
initial = 0.25
k_t = getContinuousPath(yspline, initial, T) # path of cake stock; consumer starts off with a_0=1, which has index N

fig2 = figure()
suptitle("Dynamic Paths", fontsize=16)
subplot(121)
plot(0:(T-1), k_t) 
title("Dynamic Path (initial = $initial)") 
xlabel(L"t")
ylabel(L"k")

initial = 5.0
k_t = getContinuousPath(yspline, initial, T) # path of cake stock; consumer starts off with a_0=1, which has index N

subplot(122)
plot(0:(T-1), k_t) 
title("Dynamic Path (initial = $initial)") 
xlabel(L"t")
# ylabel(L"k'")
fig2.savefig("Q3DynamicPaths.png")

