function numderiv(f,x;delta=1.0e-6)

   deriv = (f(x+delta) - f(x))/delta 

   return deriv

end

function twoSideDeriv(f,x;delta=1.0e-6)

	twoSideDeriv = (f(x+delta) - f(x-delta))/(2*delta)

   return twoSideDeriv

end

# function f(x)

#    f = exp(x)

# end

# deriv = numderiv(f,1)
# println(deriv)


