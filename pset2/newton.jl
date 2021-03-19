include("numderiv1.jl")
import LinearAlgebra: norm, pinv, Transpose, inv

# analytic Jacobian
function J(p)
	x,y = p
	J = [(y/(2*(x^2)) - 0.1*exp(0.1*x)) 1/(2y);
	     1/(2x) (x/(2*(y^2)) - 0.1*exp(0.1*y))]
end


function newton(f::Function,x_0;toler=1.0e-4,N=100)
	n=1
	x=0. # for scoping purposes
	while n<=N
		x=x_0-(f(x_0)/numderiv(f, x_0))
		if f(x)==0 || abs(x-x_0)<toler
			# println("x is $x and the iteration number is $n")
			return x
		end
		x_0=x
		n=n+1
	end
	y=f(x)
	println("Method did not converge. The last iteration gives $x with function value $y")
end


function multiDNewton(f::Function,x_0;toler=1.0e-4,N=100, lambda=0.5)
	n=1
	x=0. # for scoping purposes
	while n<=N
		x=x_0 - (1-lambda)*(inv(J(x_0))*f(x_0))
		if f(x)==0 || norm(x - x_0, Inf)<toler
			return println("x is $x and the iteration number is $n")
		end
		x_0=x
		n=n+1
	end
	y=f(x)
	println("Method did not converge. The last iteration gives $x with function value $y")
end


