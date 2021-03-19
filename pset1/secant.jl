function secant(f::Function,x_0,x_1;toler=1.0e-4,N=100)
	n=1
	x=0. # for scoping purposes
	while n<=N
		x=x_1-f(x_1)*(x_1-x_0)/(f(x_1)-f(x_0))
		if f(x)==0 || abs(x-x_1)<toler
			return println("x is $x and the iteration number is $n")
		end
		x_0=x_1
		x_1=x
		n=n+1
	end
	y=f(x)
	println("Method did not converge. The last iteration gives $x with function value $y")
end
