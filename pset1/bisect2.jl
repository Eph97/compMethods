include("io1.jl")

function f(x,d)
   f = x^3 - d
end  

function bisect(f,xlow,xhigh;toler=1.0e-6,mxiter=30)

   fxlow = f(xlow)
   fxhigh = f(xhigh)
   
   if (fxlow*fxhigh > 0) 
      wait("Root not bracketed.")
   elseif (xlow > xhigh) 
      wait("Brackets out of order.")   
   else
   
      diff = xhigh - xlow
   
      downwardsloping = (fxlow > 0)
      
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
        writeio(stdout,("Did not converge: ",(3,15.8)),xlow,xhigh,diff,callwait=true)             
     end
     
  end

  xroot = (xlow + xhigh)/2
  fxroot = f(xroot)
  
  bisect = xroot,fxroot,xlow,xhigh, niter
   
end  

# xlow = -2.0
# xhigh = 3.0
# d = 2.0 
# xroot,fxroot,xlow2,xhigh2 = bisect(x -> f(x,d),xlow,xhigh,toler=1.0e-6,mxiter=30);
# writeio(stdout,("xroot,f(xroot),xlow2,xhigh2: ",(4,15.8)),xroot,f(xroot,d),xlow2,xhigh2)

     
                 
      
      
         
   

   
