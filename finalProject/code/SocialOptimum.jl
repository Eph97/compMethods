include("interpolate.jl")
include("brent1.jl")

# using OhMyREPL #Just For Fun
using Plots, LaTeXStrings, Parameters
using Interpolations

# parameters {{{
@with_kw struct params
	ρ::Float64 #discount the future at rate ρ
	pi::Float64 # infection fatality rate
	γ::Float64 # Rate that infectiousness ends
	ν::Float64 # value of a statistical life
	δ::Float64 # rate at which a cure for the disease is found (assuming cure is 100% effective)
	β::Float64 # β = 0.3 + γ by equation (8) if Nₛ(T) = A(T) = 1
	k::Float64 # Expected cost of infection (ie π*ν)
	η::Float64 # Expected cost of infection (ie π*ν)
end

# setting parameter Values
δ=(0.67/365); ρ=(0.05/365); β=(0.3 + (1/7)); γ=(1/7.0); pi=0.0062; k=31755*0.0062
# η = 1/(ρ + δ)*(1-exp(-(ρ+δ)))
η = exp(-(ρ+δ))
p =params(ρ= ρ, pi=pi, γ=γ , ν=31755, δ=δ , β=β, k=k, η = η)
R₀ = β/γ # Assuming β > γ because otherwise the disease never breaks out.
# }}}
# R(t) = R₀*A(t)^2*(N_s(t)) # The effective reproduction number
#
# the consumer's felicity function. We allow for a+epsilon to permit activity levels of zero. Similar to our crumbs example from class.
function U(a; ϵ=0.00001)
	if a < 0 return -Inf # don't allow negative activity
	else return log(a+ϵ) - a + 1 # utility function adding epsilon to account for if a = 0
	end
end


# S_init = 0.9999223 , I_init=0.0000527
# Full SIRD inititialization and evolution {{{
function SIRD_Init(T; params = p, S_init = 0.9999223 , I_init=0.0000527, R_init = 0.47, today=false)
	S = zeros(T); I = zeros(T); R = zeros(T); D = zeros(T);
	S[1] = S_init #no social distancing before t = 0
	I[1] = I_init # Deaths before March 13, 2020
	if today == true
		R[1] = R_init + (1 - pi)*γ*I[1]
	else
		R[1] = (1 - pi)*γ*I[1]
	end
	D[1] =  pi*γ*I[1]
	return S, I, R, D
end

function SIRD_evolution(S, I,R, D, a; params=p )
	@unpack ρ, pi, γ, ν, δ, β, k, η = params
	S_new =S - β*(a^2)*S*I
	I_new =I + β*(a^2)*S*I - γ*I
	R_new = R + (1 - pi)*γ*I
	D_new = D + pi*γ*I
	return S_new, I_new, R_new, D_new
end
# }}}

# SIRD model init and evolution S, I only {{{
function SI_Init(T; params = p, S_init = 0.9999223 , I_init=0.0000527)
	S = zeros(T); I = zeros(T);
	S[1] = S_init #no social distancing before t = 0
	I[1] = I_init # Deaths before March 13, 2020
	return S, I
end

function evolution(S, I, a; params=p )
	@unpack ρ, pi, γ, ν, δ, β, k, η = params
	S_new =S - β*(a^2)*S*I
	I_new =I + β*(a^2)*S*I - γ*I
	return S_new, I_new
end
# }}}

function bellman(a, S, I, spline; params=p)
	@unpack ρ, pi, γ, ν, δ, β, k, η= params

	S_prime, I_prime = evolution(S, I, a)

	return ((S + I)*U(a) - γ*I*k) + η*spline(S_prime, I_prime)
end


# returns (discrete) value and policy functions on a grid {{{
toler=1e-6; 
tol = zeps = 1.0e-8
function valueFuncIter(A; N = length(A), toler=1e-6, maxiter=1000, verbose=1, params=p)
	v0=zeros(N, N)
	g = fill(0.0, N, N) 
	grid1 = interpGrid(0, 1, N)
	S = creategrid(0, 1, N); I = creategrid(0,1,N);
	err = 1+toler; niter = 0;
	tol = zeps = 1.0e-8
	while err > toler && niter < maxiter
		niter += 1
		v1 = copy(v0)
		interp_cubic = CubicSplineInterpolation((grid1, grid1), v0, extrapolation_bc=Line())
		# define new value function by iterating over all grid points
		for i in 1:N
			for j in 1:N
				fmin,xmin,niter = brent(0, 1/2, 1, a_j -> (-1*bellman(a_j, S[i], I[j], interp_cubic)), tol, zeps)

				v1[i, j] = -1*fmin
				g[i, j] = xmin
			end

		end
		err = maximum(abs.(v1 - v0)) # error defined using the sup norm
		v0 = v1
	end

	if err <= toler && verbose >= 1
		println("Value function iteration converged in $niter iterations (maxerr = $err)")
	elseif err > toler
		println("Value function iteration failed to converge in $niter iterations (maxerr = $err)")
	end
	return v0, g
end
# }}}


N=30
A_bar = creategrid(0, 1, N)

vsol, g = valueFuncIter(A_bar)

grid2 = interpGrid(0 , 1, N)

g_interp = CubicSplineInterpolation((grid2, grid2), g)

# functions for computing the paths {{{
function SIPath(T)
	S, I = SI_Init(T) # SI path given some Activity path.
	A = zeros(T)
	for i in 1:(T-1)
		A[i] = g_interp(S[i], I[i])
		S[i+1], I[i+1] = evolution(S[i], I[i], A[i])
	end
	A[T] = g_interp(S[T], I[T])
	return A, S, I
end

# R_init is for if we're looking at the path from today
function SIRDPath(T;S_init = 0.9999223 , I_init=0.0000527, R_init = 0.47, today=false)
	S, I, R, D = SIRD_Init(T, S_init=S_init , I_init=I_init, R_init=R_init, today=today) # SI path given some Activity path.
	A = zeros(T)
	for i in 1:(T-1)
		A[i] = g_interp(S[i], I[i])
		S[i+1], I[i+1], R[i+1], D[i+1] = SIRD_evolution(S[i], I[i],R[i], D[i], A[i])
	end
	A[T] = g_interp(S[T], I[T])
	return A, S, I, R, D
end


function PureSIRDPath(T;S_init = 0.9999223 , I_init=0.0000527, R_init = 0.47, today=false)
	S, I, R, D = SIRD_Init(T, S_init=S_init , I_init=I_init, R_init=R_init, today=today) # SI path given some Activity path.
	# S, I, R, D = SIRD_Init(T) # SI path given some Activity path.
	A = ones(T)
	for i in 1:(T-1)
		# A[i] = g_interp(S[i], I[i])
		S[i+1], I[i+1], R[i+1], D[i+1] = SIRD_evolution(S[i], I[i],R[i], D[i], A[i])
	end
	# A[T] = g_interp(S[T], I[T])
	return A, S, I, R, D
end
# }}}

# plots {{{
# plots dynamic paths of pandemic variables
function plots(T; today=false)
	if today == true
		A, S, I, R,D = SIRDPath(450, today=true, S_init=(1 - (31000/320000000 +0.47)), R_init=0.47, I_init= 31000/320000000)
	else
		A, S, I, R, D = SIRDPath(T)
	end
	RepoNumber = [R₀*(A[i]^2)*S[i] for i in 1:T]

	A_pure, S_pure, I_pure, R_pure, D_pure = PureSIRDPath(T)

	RepoNumber_pure = [R₀*(A_pure[i]^2)*S_pure[i] for i in 1:T]

	plt_A = plot(0:(T-1), A,  label=L"A_t", xlabel="days", ylabel = "A", xformatter=_->"")
	plot!(0:(T-1), A_pure,  label="SIRD", xlabel="days", xformatter=_->"")
	# plt_A = plot(0:(T-1), A,  label=L"A_t", title="Social Activity A", xlabel="days", ylabel = "A", xformatter=_->"")

	plt_S = plot(0:(T-1), S, label=L"N_s", xlabel="days", ylabel = L"N_s", xformatter=_->"")
	plot!(0:(T-1), S_pure,  label="SIRD", xlabel="days", xformatter=_->"")
	# plt_S = plot(0:(T-1), S, label=L"N_s", title="share susceptible N_s", xlabel="days", ylabel = "N_s(t)", xformatter=_->"")

	plt_I = plot(0:(T-1), I, label=L"N_i", xlabel="days", ylabel = L"N_i", xformatter=_->"")
	plot!(0:(T-1), I_pure,  label="SIRD", xlabel="days", xformatter=_->"")
	# plt_I = plot(0:(T-1), I, label=L"N_i", title="share infected N_i", xlabel="days", ylabel = "N_i(t)", xformatter=_->"")

	plt_D = plot(0:(T-1), D, label=L"N_d", xlabel="days", ylabel = L"N_d", xformatter=_->"")
	plot!(0:(T-1), D_pure,  label="SIRD", xlabel="days", xformatter=_->"")
	# plt_D = plot(0:(T-1), D, label=L"N_d", title="Deaths", xlabel="days", ylabel = "N_d(t)", xformatter=_->"")
	# R_t = R₀.*S # The effective reproduction number
	#
	plt_R = plot(0:(T-1), R, label=L"N_r", xlabel="days", ylabel = L"N_r", xformatter=_->"")
	plot!(0:(T-1), R_pure,  label="SIRD", xlabel="days", xformatter=_->"")

	plt_repo = plot(0:(T-1), RepoNumber, label=L"R_t", xlabel="days", ylabel = "Effective Reproduction Number", xformatter=_->"")
	plot!(0:(T-1), RepoNumber_pure,  label="SIRD", xlabel="days", xformatter=_->"")
	# plt_repo = plot(0:(T-1), RepoNumber, label="Rt", title="Effective Reproduction Number R", xlabel="days", ylabel = "R", xformatter=_->"")

	plot_grid = plot(plt_A, plt_S, plt_I, plt_D, plt_R, plt_repo, layout=(2,3), size=(600,600))

	return plot_grid
end


###Plots used are below and are commented out.
#
# plot1 = plots(450)
# savefig(plot1, "../plots/plots1")

# plot2 = plots(450, today=true)
# savefig(plot2, "../plots/plots2")


# V_sol_surface = surface(A_bar, A_bar, vsol, xlabel = L"N_s", ylabel = L"N_i", zlabel = "Value Function", cmap="viridis")
# savefig(V_sol_surface, "../plots/valueFunction_Surface")

# V_sol_contour = contour(A_bar, A_bar, vsol,  xlabel = L"N_s", ylabel = L"N_i", zlabel = "Value Function Contours", cmap="viridis", title="Contour Map")
# savefig(V_sol_contour, "../plots/valueFunction_Contour")

#
# g_surface = surface(A_bar, A_bar, g, xlabel = L"N_s", ylabel = L"N_i", zlabel = "decision rule", cmap="viridis")
# savefig(g_surface, "../plots/g_surface")

# g_contour = contour(A_bar, A_bar, g, xlabel = L"N_s", ylabel = L"N_i", zlabel = "decision rule", cmap="viridis", title="Contour map") 
# savefig(g_contour, "../plots/g_contour")

# }}}

