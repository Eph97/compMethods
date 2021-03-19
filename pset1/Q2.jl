include("numderiv1.jl")
# This is for problem 2a-c
function g(x)
   g = 0.5*(x^(-0.5))+0.5*(x^(-0.2))
end

function deriv(x)
   deriv = -0.25*(x^(-1.5))-0.1*(x^(-1.2))
end

println("problem 2.a")

B = Array{Float64}(undef, 10)

for i in 1:10
	B[i] = 10.0^(-i)
end
println("running over epsilon in array of values:")
println(B)

# print("Input an integer: ")
# x = parse(Float64, readline())
x = 1.5

oneSideCalcDeriv = Array{Float64}(undef, 10)
for e in zip(1:10, B)

	oneSideCalcDeriv[e[1]] = numderiv(g, x, delta=x*e[2])
	# println(oneSideCalcDeriv)
end

print("the computed derivatives are")
println(oneSideCalcDeriv)

actualDeriv = deriv(x)
print("The Actual Derivative is: ")
println(actualDeriv)
# test = Array{Float64}(actualDeriv, 10)
#
actualDerivArray = Array{Float64}(undef, 10)
for i in 1:10
	actualDerivArray[i] = actualDeriv
end
diff1 = Array{Float64}(undef, 10)
for i in 1:10
	diff1[i]= abs(oneSideCalcDeriv[i] - actualDerivArray[i])
end

minvalue = B[findmin(diff1)[2]]
println("The most accurate value is found at epsilon $minvalue")

println("2.b")
println("Runing Twosided Derivative Now")

twoSideCalcDeriv = Array{Float64}(undef, 10)
for e in zip(1:10, B)
	twoSideCalcDeriv[e[1]] =twoSideDeriv(g, x, delta=x*e[2])
	# println(CalcDeriv)
end

print("the computed derivatives are")
println(twoSideCalcDeriv)

print("The Actual Derivative is: ")
println(actualDeriv)
diff2 = Array{Float64}(undef, 10)
for i in 1:10
	diff2[i] = abs(twoSideCalcDeriv[i] - actualDerivArray[i])
end


minvalue = B[findmin(diff2)[2]]
println("The most accurate value is found at epsilon $minvalue")


println("2.c")
println(" most accurate of all 20")

combined = vcat(diff1, diff2)
combinedMin = findmin(combined)

combinedminValue = twoSideCalcDeriv[6]
minvalueEpislon =  B[6]


println("The most accurate approximate derivative of all of the ones computed in parts a and b is $combinedminValue and is found with epislon value $minvalueEpislon")
