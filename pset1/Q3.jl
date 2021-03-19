include("bisect2.jl")
include("secant.jl")
include("FuncIter.jl")

a = 0.1
function supply(p)
	supply = exp(a*p) - 1
end

b = 1.0
e = 1.0
function demand(p)
	demand = (b)*(p^(-e))
end


function excess(p)
	excess = demand(p) - supply(p)
end

xroot,fxroot,xlow,xhigh, niter = bisect(excess,0,5;toler=1.0e-4,mxiter=30)
println("root is $xroot found after $niter iterations")

secant(excess,-1,3,N=20, toler=1.0e-4)

Iter(excess, 1, toler=1.0e-4)


