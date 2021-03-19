include("FuncIter.jl")
include("newton.jl")

a_A = a_B = 0.1

b_A = b_B = 1
e_A = e_B = 1

# p_A = p_B = 1

# Supply functions
function y_A_s(p_A)
	y_A_s = exp(a_A*p_A) - 1
end
function y_B_s(p_B)
	y_B_s = exp(a_B * p_B) - 1
end

# Demand functions
function y_A_d(p_A, p_B)
	p_bar = (p_A + p_B)/2
	y_A_d = b_A * (p_A / p_bar)^(-e_A)
end

function y_B_d(p_B, p_A)
	p_bar = (p_A + p_B)/2
	y_B_d = b_B * (p_B / p_bar)^(-e_B)
end

function supply(p)
	p_A, p_B = p
	return [y_A_s(p_A); y_B_s(p_B)]
end

function demand(p)
	p_A, p_B = p
	return [y_A_d(p_A, p_B); y_B_d(p_B, p_A)]
end

function excess(p)
	excess = demand(p) .- supply(p)
end

# secant(excess,-1,3,N=20, toler=1.0e-4)
MulitIter(excess, [ 10;10 ], toler=1.0e-4)

println("Now increasing demand by 10%")

b_A = 1.1
MulitIter(excess, [10;10 ], toler=1.0e-4)

println("Using Newton's method (with analytic jacobian)")

b_A = 1
multiDNewton(excess, [ 10; 10 ], toler=1.0e-4)

println("Now increasing demand by 10%")

b_A = 1.1

multiDNewton(excess, [ 10; 10 ], toler=1.0e-4)
