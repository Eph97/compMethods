include("interpolate.jl")
include("brent1.jl")
include("gauher1.jl")
include("random2.jl")

# using Random
using PyPlot
using LaTeXStrings
# BC
# (1+r)*a_t + y_t = c_t+ a_t+1

# y_0 = b

β = 0.7; r = 0.2; b= 1.0; σ = 0.5; μ=0; tol = zeps = 1.0e-8

function createGrid(xmin, xmax, N)
    step = (xmax-xmin) / (N-1)
    return collect(xmin:step:xmax)
end


function bellman(U, s_t, a_t, v::Function, g::Function=identity; gauss=false)
	if (gauss ==true)
		return U(s_t - a_t) + β*v(g, a_t)
	else
		return U(s_t - a_t) + β*v(a_t)
	end

end


# the consumer's felicity function
function U(c)
    # we don't allow negative consumption
    if c < 0
        return -Inf
    else
        return log(c)
    end
end


function summand(g::Function, x ,w, a)
	# ((1+r)*a + b + exp(σ*x*sqrt(2)))
	return w*g(((1+r)*a + b + exp(σ*x*sqrt(2)))
)
end

function G(g::Function, a; n=3)
	w, x = gauher(n)
	return π^(-1/2)*sum(summand.(g, x,w, a))
end

N = 21 # number of grid points
s = createGrid(b,  b +  exp(3*σ*sqrt(2)), N)


v3(s) = U(s)
g3(s) = 0
# v2(s) = bellman(U, s, a_t, v3)
v2(s) = brent(0, (s/2), s , a -> (-1*bellman(U, s, a, v3)), tol,zeps)[1]
g2(s) = brent(0, (s/2), s , a -> (-1*bellman(U, s, a, v3)), tol,zeps)[2]

v1(s) =  brent(0, (s/2), s , a -> (-1*bellman(U, s, a, G, v2, gauss=true)), tol,zeps)[1]
g1(s)  =  brent(0, (s/2), s , a -> (-1*bellman(U, s, a, G, v2, gauss=true)), tol,zeps)[2]

# v0(s) = brent(0, (s/2), s , a -> (-1*bellman(U, s, a, G, v1, gauss=true)), tol,zeps)[1]
# g0(s) = brent(0, (s/2), s , a -> (-1*bellman(U, s, a, G, v1, gauss=true)), tol,zeps)[2]
# fmin,xmin,niter = v1(s)
# vsol, g = valueFuncIter(U, a)
#
#
v1_sol = -1*v1.(s)
v2_sol = -1*v2.(s)
v3_sol = v3.(s)

g1_sol = g1.(s)
g2_sol = g2.(s)


vspline = makespline(s,v1_sol)
v1_interp(s_t) = interp(s_t,vspline)[1]

v0(s) = (-1*brent(0, (s/2), s , a -> (-1*bellman(U, s, a, G, v1_interp, gauss=true)), tol,zeps)[1])
g0(s) = brent(0, (s/2), s , a -> (-1*bellman(U, s, a, G, v1_interp, gauss=true)), tol,zeps)[2]


#2.a
println("g_0 is ", g0(b))
println("v_0 is ", v0(b))

fig1 = figure(figsize=(8, 8))
subplot(321)
plot(s, v1_sol, label=L"v_1")  
plot(s, v2_sol, label=L"v_2")  
plot(s, v3_sol, label=L"v_3")  
title("Value Function ( N= $N)") 
xlabel(L"S")
ylabel(L"v(S)")
legend()

# subplot(222)
# title("Value Function ( N= $N)") 
# xlabel(L"S")
# ylabel(L"v_2(S)")


# subplot(223)
# title("Value Function ( N= $N)") 
# xlabel(L"S")
# ylabel(L"v_3(S)")

# fig2 = figure()
subplot(325)
plot(s, g1_sol, label="decision rule")  
plot(s,s, label="45 deg line") 
title("Decision Rule") 
xlabel(L"S")
ylabel(L"g_1(S)")
legend()

subplot(326)
plot(s, g2_sol, label="decision rule")  
plot(s,s, label="45 deg line") 
title("Decision Rule") 
xlabel(L"S")
ylabel(L"g_2(S)")
legend()
# fig1.savefig("2_a.png")


#2.b 
# simulate 

function s_t(a_t; seed=333, n=1)
	z, seed = randomnormal(seed,n)
	u = σ*z + μ
	return ((1 + r)*a_t + b + exp(σ*u*sqrt(2))), seed
end

function sim(seed)
	a_t = fill(0.0, 4)
	s_sim = fill(b, 4) 

	a_t[1] = g0(s_sim[1])
	s_sim[2], seed = s_t(a_t[1], seed=seed)

	a_t[2] = g1(s_sim[2])
	s_sim[3],seed = s_t(a_t[2], seed=seed)

	a_t[3] = g2(s_sim[3])
	s_sim[4] = a_t[3]*(1+r)

	a_t[4] = g3(s_sim[4])

	c_t = s_sim .- a_t
	return c_t, a_t, seed
end

fig3 = figure(figsize=(8,8))

ax1 = subplot(121)
title(L"c_t simulation") 
xlabel(L"t")
ylabel(L"c_t")

ax2 = subplot(122)
title(L"a_{t+1} simulation") 
xlabel(L"t")
ylabel(L"a_{t+1}")

function simulate(n; seed = 333)
	for i in 1:n
		c_t, a_t, seed = sim(seed)
		ax1.plot(0:3, c_t) 

		ax2.plot(0:3, a_t) 
	end
end

simulate(10)

# fig3.savefig("2_b.png")

