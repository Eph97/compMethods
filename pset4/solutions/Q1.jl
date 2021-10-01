include("derivs.jl")
# using Plots
using LaTeXStrings

# beta = 0.9
# k_bar = twoSideDeriv(F,inv(beta))
# k_bar = ((0.9*0.4)/(1 - 0.9+ 0.9* 0.1))^(1/(1-0.4))
# returns grid of N equally-spaced points on interval [xmin, xmax]
function createGrid(xmin, xmax, N)
	# xmax *= k_bar
	# xmin *= k_bar
    step = (xmax-xmin) / (N-1)
    return collect(xmin:step:xmax)
end

# returns (discrete) value and policy functions on grid defined by the vector a
function valueFuncIter(U, a; β=0.9, v0=zeros(length(a)), toler=1e-6, maxiter=1000, verbose=1)
    N = length(a)
    g = fill(0, N) # g[i] will be j<=i that maximizes U(a_i-a_j) + βv_j (this is the optimal decision rule/policy function)
    U_ar = [U(F(a[i])-a[j]) for i in 1:N, j in 1:N] # NxN array of all possible values of U, given grid a
    err = 1+toler; niter = 0;
    while err > toler && niter < maxiter
        niter += 1
        v1 = copy(v0)

        # define new value function by iterating over all grid points
        for i in 1:N
            rhs = U_ar[i,:] .+ β*v0[:] # RHS of Bellman equation (for all a_j<=a_i)
            v1[i], g[i] = findmax(rhs) # simultaneously update value function and get the optimizing j (for given i)
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

# function objective(a_i, a_j, v_j)
# 	return utility(F(a_i) - a_j) + (beta* v_j)
# end

# the consumer's felicity function
function U(c)
    # we don't allow negative consumption, even if c > -ϵ
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



# N = 31 # number of grid points
# k = createGrid(0.5, 1.5, N)
# vsol, g = valueFuncIter(U, k)




# T = 41
# k_t = k[getPath(g, N, T)] # path of cake stock; consumer starts off with a_0=1, which has index N
# # c_t = f.(a_t[1:(T-1)]) .- (1 - delta) .* a_t[2:T]  # path of consumption
# c_t = F.(k_t[1:(T-1)]) .- k_t[2:T]

# lifetimeUtil_t = lifetimeU(c_t)

# difference = lifetimeUtil_t - vsol[N]

# dist = maximum(abs.(lifetimeUtil_t - vsol

# plot_k1 = plot(0:(T-2), [k_t[1:(T-1)] k_t[2:T]], title="N = $N",label=[L"k_t" L"k_{t+1}"], xlabel="a" ,ylabel="a'", legend=:topleft)
# savefig(plot_k1, "plot_k1")
# plot_a1 = plot(a, vsol, title="Value function (N = $N)", xlabel=L"k", ylabel=L"V(k)", legend=false)
# plot_a2 = plot(k, [k[g] k], label=["Optimal decision rule" "45-degree line"], title="Optimal decision rule (N = $N)", xlabel=L"k", ylabel=L"k'", legend=:topleft)


# plot_b = plot(0:(T-2), [k_t[1:(T-1)] c_t[1:(T-1)]], label=[L"k_t" L"c_t"], title="Simulation (N = $N)", xlabel=L"t", marker=3)
# plot_c = plot(0:(T-2), c_t[1:(T-1)], label=L"c_t", title="consumption (N = $N)", xlabel=L"t", marker=3)
# savefig(plot_c, "plot_c")

#
# plot(plot_k1, plot_c , layout=2, size = (800,600))
