import LinearAlgebra: norm

using Plots


eps = 0.01
beta = 0.9

function Q2(N; toler = 1.0e-3, Maxiter = 500)
	A = LinRange(0, 1, N)
	v = fill(5.0, 1,N)
	v_1 = fill(1.0, 1,N)
	A_prime = fill(0.0, N, 2)
	for M in 1:Maxiter
		for i in 1:N
			mxval, mxindx = g(A, v, i)
			v_1[i] = mxval
			A_prime[i,1] = A[mxindx]
			A_prime[i,2] =mxindx
			# println(mxval, " ", mxindx)
		end
		if norm(v_1 - v, Inf)<toler
			println("converged after $M iterations")
			# println(A_prime)
			println("a'   (index)")
			# println("v", v)
			return A, A_prime, v
		else
			v = copy(v_1)
		end
	end
	println("Did Not Converge")
end

function g(A, v, i)
	mxval, mxindx = findmax(utility.(A[i], A[1:i], v[1:i]))
	return mxval, mxindx
end

function utility(a_i, a_j, v_j)
	return log(a_i - a_j + eps) + (beta* v_j)
end

function consumption(N, a, a_p)
	c_t = fill(0.0, N, 2)
	n = copy(N)
	for i in 1:N
		# println("n is $n, a[n] is $(a[n]) a_p[n,1] is $(a_p[n, 1])")
		c_t[i, 1] = a[n]
		c_t[i, 2] = a[n] - a_p[n, 1]
		n = convert(Int64, a_p[n, 2])
	end
	return c_t
end


N=21
A, A_p, v = Q2(N)
c_t = consumption(N, A, A_p)

# plot_a1 = plot(A, [A_p[:,1] A], title="N = $N",label=["uniform distribution" "optimal path"], xlabel="a" ,ylabel="a'", legend=:topleft)
# savefig(plot_a1, "plot_a1")

# plot_v1 = plot(A, v[:],  title="value function N = $N", xlabel="a" ,ylabel="v(a)",label="v(a)", legend=:topleft)
# savefig(plot_v1, "plot_v1")

# plot_c1 = plot(1:N, [c_t[:,1] c_t[:, 2]],  title="consumption N = $N", legend=:topright, label=["a_t" "c_t"], marker=3)
# savefig(plot_c1, "plot_c1")


# N=51
# A, A_p, v = Q2(N)
# c_t = consumption(N, A, A_p)

# plot_a2 = plot(A, [A_p[:,1] A],  title="N = $N", xlabel="a" ,ylabel="a'", label=["uniform distribution" "optimal path"], legend=:topleft)
# savefig(plot_a2, "plot_a2")

# plot_v2 = plot(A, v[:],  title="value function N = $N", xlabel="a" ,ylabel="v(a)", label="v(a)", legend=:topleft)
# savefig(plot_v2, "plot_v2")

# plot_c2 = plot(1:convert(Int64,(N+1)/2), [c_t[1:convert(Int64,(N+1)/2),1] c_t[1:convert(Int64, (N+1)/2), 2]],  title="consumption N = $N", label=["a_t" "c_t"], xlabel="t", legend=:topright, marker=3)

# savefig(plot_c2, "plot_c2")


# fn = plot(plot_a1, plot_a2, plot_v1, plot_v2, plot_c1, plot_c2, layout=(3, 2))
# savefig(fn)
