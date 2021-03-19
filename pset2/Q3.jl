include("zbrent1.jl")
include("newton.jl")

prob = .5
yH = 10
yL = 0

beta=0.9
y_0 =20
# y_1 =5 # is now yH and yL
r=0.1


function felicity(c)
	return log(c)
end

function foc_felicity(c)
	return 1/(c)
end

function z(a_1; a_0=5)
	return -foc_felicity(a_0 + y_0 - a_1) + (beta)*prob*(1+r)foc_felicity((1+r)*a_1 +yH) + (beta)*(1-prob)*(1+r)foc_felicity((1+r)*a_1 +yL)
end

function g(a_1; a_0=5)
	y_1 =5
	return -foc_felicity(a_0 + y_0 - a_1) + (1+r)*(beta)*foc_felicity((1+r)*a_1 +y_1)
end

function value(a_0)
	a_1_optimal = newton(a_1 ->g(a_1, a_0=a_0),2)
	y_1 = 5
	value = felicity(a_0 + y_0 - a_1_optimal) + beta*(1+r)*felicity((1+r)*a_1_optimal +y_1) + beta*(1+r)foc_felicity((1+r)*a_1_optimal + y_1)
	println("a_0: $a_0, optimal a_1 is $a_1_optimal, with value: $value")
	return value
end

function value_with_uncertainty(a_0)
	a_1_optimal = newton(a_1 -> z(a_1, a_0=a_0),2)
	value = felicity(a_0 + y_0 - a_1_optimal) + beta*prob*(1+r)*felicity((1+r)*a_1_optimal +yH) + beta*(1-prob)*(1+r)foc_felicity((1+r)*a_1_optimal +yL)
	println("a_0: $a_0, optimal a_1 is $a_1_optimal, with value: $value")
	return value

end

As = [0 1 2 3 4 5 6 7 8 9 10]

println("Using Newton's method.")

a_1_optimal = newton(z,2)

println("optimal a_1 is $a_1_optimal")
newton(z,2)

println("Using Brent's method.")
zbrent(z,0,13,1.0e-4,1.0e-4)


# for part c
println("values without uncertainty")
value.(As)

println("With uncertainty")
value_with_uncertainty.(As)



