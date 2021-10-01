# include("../../derivs.jl")
include("../../random1.jl")

using Plots
using LaTeXStrings
using Statistics

pyplot() 

β = 0.9; α=0.4; δ=0.1; p = 0.9
# k_bar = ((0.9*0.4)/(1 - 0.9+ 0.9* 0.1))^(1/(1-0.4))
states = [1.1 0.9]
Z = length(states)
k_ss(z) = ((1/β-(1-δ))/(z*α))^(1/(α-1)) # analytic steady state (solution to F'(k)=1/β)
k_bar = k_ss(1)

Ev = [(i==j ? p : 1 - p) for i in 1:2, j in 1:2] 
as_ints(a::AbstractArray{CartesianIndex{L}}) where L = reshape(reinterpret(Int, a), (L, size(a)...))


# returns grid of N equally-spaced points on interval [xmin, xmax]
function createGrid(xmin, xmax, N)
	step = (xmax-xmin) / (N-1)
	return collect(xmin:step:xmax)
end

# returns (discrete) value and policy functions on grid defined by the vector a
function valueFuncIter(U, k; β=0.9, N = length(k[:,1]),states = states, v0=zeros(N, length(states)),  toler=1e-6, maxiter=1000, verbose=1)
	# seed = 333; n=1
	Z = length(states)
	# N = length(k[:,1])
	g = fill(0, N, Z) # g[i] will be j<=i that maximizes U(a_i-a_j) + βv_j (this is the optimal decision rule/policy function)
	U_ar = [U(F(k[i], states[z])-k[j]) for i in 1:N, j in 1:N, z in 1:Z] # NxN array of all possible values of U, given grid a
	err = 1+toler; niter = 0;
	while err > toler && niter < maxiter
		niter += 1
		v1 = copy(v0)
		for i in 1:N
		    rhs = U_ar[i,:,:] .+ β*v0*Ev # RHS of Bellman equation (for all feasible j)
		    v1[i,:], ind = findmax(rhs, dims=1) # update value function and get the optimizing j
		    g[i,:] = as_ints(ind)[1,:,:]
		end

		err = maximum(abs.(v1-v0)) # error defined using the sup norm
		v0 = v1
	end

	if err <= toler && verbose >= 1
		println("Value function iteration converged in $niter iterations (maxerr = $err)")
	elseif err > toler
		println("Value function iteration failed to converge in $niter iterations (maxerr = $err)")
	end

	return v0, g
end


function f(k, z; alpha=0.4, testing = 0)
	if testing == 1
		println( z)
	end
	return z*k^alpha
end

function F(k, z; alpha=0.4, delta=0.1)
	return f(k, z) + (1 - delta)*k
end

# the consumer's felicity function
function U(c)
	if c < 0
		return -Inf
	else
		return log(c)
	end
end

function lifetimeU(c; beta=0.9)
	total = 0
	for t in 1:(length(c))
		total += beta^(t-1) * U(c[t])
	end
	return total
end

# keeps in same state with probability π = 0.9 and changes states with probability 
index(r, z) = (r > 0.9) ? ((z % Z) + 1) : z

# returns the indices (in a) corresponding to the optimal path of savings, when a0 = a[idx0]
function getPath(g, idx0, T; n=1, z = 1, seed = 333)
	# r,seed = randuniform(seed,n)
	# z = index(r, z)
	idxs = [fill(idx0, T) fill(z, T)] # optimal path of cake stock (expressed as an index of a)
	for t in 2:T
		idxs[t, 1] = g[idxs[t-1, 1], z]
		idxs[t, 2] = z
		r,seed = randuniform(seed,n)
		z = index(r, z)
	end
	return idxs
end
N = 101 # number of grid points
k = createGrid(0.25, 1.75, N) .*k_bar
vsol, g = valueFuncIter(U, k)


T = 1000

k0 = 4 # initial stock of capital
z = 1
i0 = findfirst(k .>= k0) # index of k (roughly) corresponding to k0
simulatedPath = getPath(g, i0, T+2, z=1)
k_t = k[simulatedPath[:,1]] # path of capital
f_t = f.(k_t[1:(T+1)],  states[simulatedPath[1:(T+1),2]])  # output
c_t = F.(k_t[1:(T+1)],  states[simulatedPath[1:(T+1),2]]) - k_t[2:(T+2)] # path of consumption
I_t = k_t[2:(T+2)] - (1 - δ) .* k_t[1:(T+1)] # path of investment ( recall k_{t+1} = (1 - δ)k_t + I_t where I_t can also be written s Y_t

t=100

# # (a)
plot_a1 = plot(k, vsol, label=["value function good state" "value function bad state"], title="Value function", xlabel=L"k", ylabel=L"V(k)", legend=:topleft)

plot_a2 = plot(k, [k[g] k], label=["Policy rule, good state" "Policy rule, bad state" "45-degree line"], title="Policy rule", xlabel=L"k", ylabel=L"k'", legend=:topleft)

# # (b)
# # capital stock and consumption
plot_b1 = plot(0:t, [k_t[1:t+1]],  title=L"Simulated path of $k$ ($k_0=%$k0$)", 
	       xlabel=L"t", label=L"k_t", marker=3)

# output
plot_b2 = plot(0:t, [f_t[1:t+1]],  title=L"Simulated path of output $f$ ($k_0=%$k0$)", 
	       xlabel=L"t", label=L"f_t ", marker=3, legend=:topright)

# # hline!([k_bar], label=L"k_{ss}^g", linestyle=:dash) # add steady-state k
# # (b2) output, saving
plot_b3 = plot(0:t, [c_t[1:t+1] ],  title=L"consumption ($k_0=%$k0$)", 
	       xlabel=L"t", label=L"c_t", marker=3, legend=:bottomright)

plot_b4 = plot(0:t, [I_t[1:t+1] ],  title=L"Investment($k_0=%$k0$)", 
	       xlabel=L"t", label=L"I_t", marker=3, legend=:topright) 
# legend!(bbox_to_anchor=[1.05,1],loc=2,borderaxespad=0)

# plot_a = plot(plot_a1, plot_a2, layout=2, size = (800,600))
# savefig(plot_a, "plot_a")

# plot_b = plot(plot_b1, plot_b2, plot_b3, plot_b4, layout=4, size = (800,600))
# savefig(plot_b, "plot_b")

plot(plot_a1, plot_a2, plot_b1, plot_b2, plot_b3, plot_b4, layout=(3,2), size = (800,600))
println("\n k_t: mean = $(mean(k_t)) \n std = $(stdm(k_t, mean(k_t))) \n std/mean $(stdm(k_t, mean(k_t))/ mean(k_t))")
println("\n f_t: mean = $(mean(f_t)) \n std = $(stdm(f_t, mean(f_t))) \n std/mean $(stdm(f_t, mean(f_t))/ mean(f_t))")
println("\n c_t: mean = $(mean(c_t)) \n std = $(stdm(c_t, mean(c_t))) \n std/mean $(stdm(c_t, mean(c_t))/ mean(c_t))")
println("\n I_t: mean = $(mean(I_t)) \n std = $(stdm(I_t, mean(I_t))) \n std/mean $(stdm(I_t, mean(I_t))/ mean(I_t))")
