using Formatting, LinearAlgebra

function second_deriv(f,x; delta=1.0e-6)
	secondDeriv = (twoSideDeriv(f, x+delta) - twoSideDeriv(f, x - delta)) / (2 * delta)
   return secondDeriv
end

function twoSideDeriv(f,x;delta=1.0e-6)

	twoSideDeriv = (f(x+delta) - f(x-delta))/(2*delta)

   return twoSideDeriv

end


function newton_raphson(f::Function,x_0;toler=1.0e-4,N=100)
	n=1
	x=0. # for scoping purposes
	while n<=N
		# x=x_0-(twoSideDeriv(f, x_0)/second_deriv(f, x_0))
		x=x_0 - inv(H(f,x_0))*grad(f,x_0)
		if f(x)==0 || maximum(abs.(x-x_0))<toler
			# println("x is $x and the iteration number is $n")
			return f(x), x, n
		end
		x_0=x
		n=n+1
	end
	y=f(x)
	println("Method did not converge. The last iteration gives $x with function value $y")
end


# computes the numerical Jacobian for a multidimensional function f: R^N -> R^M
function numjacobian(f::Function, x; h=1e-6)
    N = length(x); M = length(f(x))
    if (N == 1 && M == 1) # f: R -> R, just take two-sided numerical derivative
        return (f(x+h) - f(x-h)) / (2h)
    else
        vv = zeros(N,N) + Diagonal(fill(h,N)) # matrix of basis column vectors (identity matrix) scaled by h
        return [(f(x+vv[:,i])[j] - f(x-vv[:,i])[j]) / (2h) for i in 1:N, j in 1:N]
    end
end

function grad(f::Function, x; h=1e-6)
    N = length(x); M = length(f(x))
    if (N == 1 && M == 1) # f: R^n -> R, just take two-sided numerical derivative
        return (f(x+h) - f(x-h)) / (2h)
    else
        vv = zeros(N,N) + Diagonal(fill(h,N)) # matrix of basis column vectors (identity matrix) scaled by h
        return [(f(x+vv[:,i]) - f(x-vv[:,i])) / (2h) for i in 1:N]
    end
end

function H(f::Function, x; h=1e-6)
    N = length(x); M = length(f(x))
    if (N == 1 && M == 1) # f: R^n -> R, just take two-sided numerical derivative
        return (twoSideDeriv(f, x+h) - twoSideDeriv(f, x - h)) / (2 * h)
    else
        vv = zeros(N,N) + Diagonal(fill(h,N)) # matrix of basis column vectors (identity matrix) scaled by h
	return [(grad(f, (x+vv[:,i]))[j] - grad(f, (x-vv[:,i]))[j]) / (2h) for i in 1:N,  j in 1:N]
    end
end
