include("zbrent1.jl")
include("newton.jl")

beta=0.9
a_0=5
y_0 =20
y_1 =5
r=0.1

function felicity(c)
	return log(c)
end

function foc_felicity(c)
	return 1/(c)
end


function g(a_1)
	return -foc_felicity(a_0 + y_0 - a_1) + (1+r)*(beta)*foc_felicity((1+r)*a_1 +y_1)
end


function analytic_a_1() 
	analytic_a_1= (beta * (1+r)*a_0 + beta*(1+r)*y_0 - y_1) / ((1+r)*(1+beta))
end

ans =  analytic_a_1()

println("Analytic answer $ans")



println("Using Newton's method.")
println(newton(g,2))

println("Using Brent's method.")
zbrent(g,0,10,1.0e-4,1.0e-4)


