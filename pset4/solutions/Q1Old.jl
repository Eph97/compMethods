include("derivs.jl")
# using Plots


# eps = 0.01
alpha = 0.4
beta = 0.9
delta = 0.1

function Q1(N; toler = 1.0e-3, Maxiter = 500)
	A = LinRange(0.5, 1.5, N)
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
			println("a'   (index)")
			println(A_prime)
			return A, A_prime, v
		else
			v = copy(v_1)
		end
	end
	println("Did Not Converge")
end

# optimal mapping
function g(A, v, i)
	mxval, mxindx = findmax(objective.(A[i], A[1:i], v[1:i]))
	return mxval, mxindx
end

function f(k)
	return k^alpha
end

function F(k)
	return f(k) + (1 - delta)*k
end

function objective(a_i, a_j, v_j)
	return utility(F(a_i) - a_j) + (beta* v_j)
end

function utility(c_t)
	if c_t < 0
		return -Inf
	else
		return log(c_t)
	end
end

function consumption(N, a, a_p)
	c_t = fill(0.0, N, 2)
	n = copy(N)
	for i in 1:N
		println("n is $n, a[$n] is $(a[N]) a_p[$N,1] is $(a_p[N, 1])")
		c_t[i, 1] = a[n]
		c_t[i, 2] = a[n] - a_p[n, 1]
		n = convert(Int64, a_p[n, 2])
	end
	return c_t
end


N=21
A, A_P, V = Q1(N)
C_T = consumption(N, A, A_P)

println("check steady state condition β ")
x = 2.0
abs(twoSideDeriv(F,x) - inv(beta)) < 1.0e-3


twoSideDeriv(F,inv(beta))
