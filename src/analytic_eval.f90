double precision function analytic_eval(mld, z, nml, n1, ndep, w) 
  !      -----------------------------------------------------            
  !      DESCRIPTION                                                      
                                                                        
  !      Evaluates background density profile according                   
  !      to analytic buoyancy frequency function.                         
                                                                        
  !      __________________                                               
  !      ------------------                                               
  !      SET VARIABLES                                                    
                                                                        
  real*8, parameter ::PI=3.14159265358979323846
  real*8 :: mld, z, nml, n1, ndep, w, a, b, c, e, d,f,g,h 
  !      __________________                                               
  !      ------------------                                               
  !      EVALUATE FUNCTION                                                
                                                                        
       h = (3.14159265358979323846**(0.5)) 
       g = (mld + 3*w - z)/w 
       a = z*ndep 
       b = w*(log(tanh(g) + 1) + (2*z)/w) 
       c = ndep/2 - nml/2 
       d = h*n1*erf(g) 
       e = (2/w) 
                                                                        
       f = a - (b * c) - (d / e) 
                                                                        
                                                                        
       analytic_eval = (1026.0d0/9.8d0) * f + 1026.0 
!      ___________________                                              
                                                                        
       return 
      END                                           
