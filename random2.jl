using Random

# This function generates random variables uniformly distributed on the [0,1] interval.
# The function returns (r,seed), where seed is the updated seed for the random number
# generator and r is either a real number (if n = 1) or an array (if n > 1).
function randuniform(seed,n)
 
   d2p31m = 2147483647.0
   d2p31 = 2147483711.0
   dmultx = 16807.0
   
   r = zeros(n)
   
   for i in 1:n
      seed = dmultx*seed % d2p31m
      r[i] = seed/d2p31
   end
   
   if (n == 1)
      return r[1],seed
   else
      return r,seed   
   end   
   
end   

# This function generates random variables with a N(0,1) distribution.  The function returns 
# (r,seed), where seed is the updated seed for the random number generator and r is 
# either a real number (if n = 1) or an array (if n > 1). 
function randomnormal(seed,n)

   d2p31m = 2147483647.0
   d2p31 = 2147483711.0
   dmultx = 16807.0
   
   r = randn(MersenneTwister(Int(seed)),n)
   
   seed = dmultx*seed % d2p31m
    
   if (n == 1)
      return r[1],seed
   else
      return r,seed   
   end   
   
end   
   
   



      
