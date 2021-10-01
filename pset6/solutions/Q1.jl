include("interpolate.jl")

using PyPlot
using Statistics

function f(x)
   return sin(x)
end   

  
a = 0.0
b = 2.0*π

N = 5
# Create the grid.
x = creategrid(a,b,N)
# Calculate the function values on the grid.
y = f.(x)
# Compute the cubic spline using x and y.
yspline = makespline(x,y)

#
nxfine = 200
xfine = creategrid(a,b,nxfine)
ysplinefine = broadcast(x->interp(x,yspline)[1],xfine)
ylinfine = broadcast(x->linearinterp(x,yspline)[1],xfine)
yfine = f.(xfine)


fig1 = figure()
# suptitle("Dynamic Paths", fontsize=16)
subplot(321)
plot(x,y, linestyle="", marker="o")
# actual values
plot(xfine, yfine, label="actual values")
# interpolated values
plot(xfine, ysplinefine, label="cubic spline interpolation")
plot(xfine, ylinfine, label="linear interpolation")

title("interpolated values (N = $N)") 
xlabel(L"x")
ylabel(L"y")

splinediff = ysplinefine - yfine
maxdiff = maximum(abs.(splinediff))
meandiff = mean(abs.(splinediff))

println("At N = $N, maxmimum difference for cubic spline is $maxdiff and mean difference is $meandiff")

lindiff = ylinfine - yfine
maxdiff = maximum(abs.(lindiff))
meandiff = mean(abs.(lindiff))

println("At N = $N, maxmimum difference for linear interpolation is $maxdiff and mean difference is $meandiff")

subplot(322)
plot(xfine, splinediff, label="cubic spline difference")
plot(xfine, lindiff, label="Linear Interp difference")
title("differences")

N = 11
x = creategrid(a,b,N)
y = f.(x)
yspline = makespline(x,y)


nxfine = 200
xfine = creategrid(a,b,nxfine)
ysplinefine = broadcast(x->interp(x,yspline)[1],xfine)
ylinfine = broadcast(x->linearinterp(x,yspline)[1],xfine)
yfine = f.(xfine)

subplot(323)
plot(x,y, linestyle="", marker="o")
# actual values
plot(xfine, yfine, label="actual values")
# interpolated values
plot(xfine, ysplinefine, label="cubic spline interpolation")
plot(xfine, ylinfine, label="linear interpolation")

title("interpolated values (N = $N)") 
xlabel(L"x")
ylabel(L"y")


splinediff = ysplinefine - yfine
maxdiff = maximum(abs.(splinediff))
meandiff = mean(abs.(splinediff))

println("At N = $N, maxmimum difference for cubic spline is $maxdiff and mean difference is $meandiff")

lindiff = ylinfine - yfine
maxdiff = maximum(abs.(lindiff))
meandiff = mean(abs.(lindiff))

println("At N = $N, maxmimum difference for linear interpolation is $maxdiff and mean difference is $meandiff")

subplot(324)
plot(xfine, splinediff, label="cubic spline difference")
plot(xfine, lindiff, label="Linear Interp difference")
title("differences")

N = 21
# Create the grid.
x = creategrid(a,b,N)
# Calculate the function values on the grid.
y = f.(x)
# Compute the cubic spline using x and y.
yspline = makespline(x,y)

nxfine = 200
xfine = creategrid(a,b,nxfine)
ysplinefine = broadcast(x->interp(x,yspline)[1],xfine)
ylinfine = broadcast(x->linearinterp(x,yspline)[1],xfine)
yfine = f.(xfine)


# fig3 = figure()
# suptitle("Dynamic Paths", fontsize=16)
subplot(325)
plot(x,y, linestyle="", marker="o")
# actual values
plot(xfine, yfine, label="actual values")
# interpolated values
plot(xfine, ysplinefine, label="cubic spline interpolation")
plot(xfine, ylinfine, label="linear interpolation")

title("interpolated values (N = $N)") 
xlabel(L"x")
ylabel(L"y")

legend()
# close()

splinediff = ysplinefine - yfine
maxdiff = maximum(abs.(splinediff))
meandiff = mean(abs.(splinediff))

println("At N = $N, maxmimum difference for cubic spline is $maxdiff and mean difference is $meandiff")

lindiff = ylinfine - yfine
maxdiff = maximum(abs.(lindiff))
meandiff = mean(abs.(lindiff))

println("At N = $N, maxmimum difference for linear interpolation is $maxdiff and mean difference is $meandiff")
subplot(326)
plot(xfine, splinediff, label="cubic spline difference")
plot(xfine, lindiff, label="Linear Interp difference")
title("differences")
legend()
# plot(xfine, diff)
