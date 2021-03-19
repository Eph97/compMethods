include("io1.jl")

d = 1.0

function f(x)
   f = x^3 - d
end  

function bisect(f,a,b;toler=1.0e-6,mxiter=30)

   fa = f(a)
   fb = f(b)
   
   if (fa*fb > 0) 
      wait("Root not bracketed.")
   elseif (a > b) 
      wait("Brackets out of order.")   
   else
   
      diff = b - a
   
      downwardsloping = (fa > 0)
      
      xlow = a
      xhigh = b
      niter = 0
      
      while (niter < mxiter)
      
         niter += 1

         xcur = (xlow + xhigh)/2

         fxcur = f(xcur)

         if downwardsloping
            if (fxcur > 0) 
               xlow = xcur
            else
               xhigh = xcur
            end
         else
           if (fxcur > 0) 
              xhigh = xcur
           else 
              xlow = xcur
           end
        end               
          
        diff = abs(xhigh-xlow)  
        
        writeio(stdout,(4,(5,15.8)),niter,xlow,xhigh,xcur,fxcur,diff,callwait=false)        
        
        if (diff < toler) 
           break
        end   
        
     end

     if (diff > toler)
        writeio(stdout,("Did not converge: ",(3,15.8)),xlow,xhigh,diff,callwait=false)             
     end
     
  end
  
  bisect = (xlow+xhigh)/2
   
end  
     
                 
      
      
         
   

   
