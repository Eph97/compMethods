import LinearAlgebra: norm, inv

function g(f, x)
	g = f(x) .+ x

end

function MulitIter(f::Function, x; toler=1.0e-6, N=100)
	n=1
	while n < N
		if norm((g(f,x) .- x), Inf) < toler
			return println("p is $x and the iteration number is $n")
		end
		x = g(f,x)
		n += 1
	end
	y = f(x)
	println("Method did not converge. The last iteration ($n) gives $x with function value $y")
end
