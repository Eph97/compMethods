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
      
