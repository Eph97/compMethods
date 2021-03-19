include("io1.jl")

# Implement Brent's method for finding the root of a one-dimensional function.  See 
# Chapter 9.3 in Numerical Recipes in Fortran 77 for a detailed description of the function 
# zbrent.

function zbrent(fct,x1,x2,rtol,ftol)
# x1 and x2 must bracket the root.  rtol and ftol are convergence tolerances, the 
# first for the root itself and the second for the function value.  zbrent returns
# the value of the root.

   global d,e

   itmax = 100

   toler = 1.0e-12

   tol1 = rtol

   a = x1
   b = x2
   fa = fct(a)
   fb = fct(b)
   if (fb*fa > 0.0)
      writeio(stdout,("Root must be bracketed in zbrent: ",(4,15.8)),a,b,fa,fb,callwait=true)
   end
   
   c = b
   fc = fb

   for iter in 1:itmax
   
      if (fb*fc > 0.0) 
         c = a
         fc = fa
         d = b - a
         e = d
      end
      
      if (abs(fc) < abs(fb))
         a = b
         b = c
         c = a
         fa = fb
         fb = fc
         fc = fa
      end
      
      xm = 0.5*(c-b)
      
      if ((abs(xm) <= tol1) | (abs(fb) <= ftol))
         return b
      end

      if ( (abs(e) >= tol1) & (abs(fa) > abs(fb)) )
         s = fb/fa
         if (abs(c-a) <= toler) 
            p = 2.0*xm*s
            q = 1.0 - s
         else
            q = fa/fc
            r = fb/fc
            p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
            q = (q-1.0)*(r-1.0)*(s-1.0)
         end
         if (p > 0) 
            q = -q
         end   
         p = abs(p)
         term1 = 2.0*p
         term2 = min(3.0*xm*q-abs(tol1*q),abs(e*q))
         if (term1 < term2) 
            e = d
            d = p/q
         else
            d = xm
            e = d
         end
      else
         d = xm
         e = d
      end
      a = b
      fa = fb
      if (abs(d) > tol1) 
         b = b + d
      else
         b = b + abs(tol1)*sign(xm)
      end
      fb = fct(b)
         
   end

   wait("zbrent exceeding maximum number of iterations")

end    

