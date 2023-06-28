subroutine numflux_x(ConL, ConR, PrimL, PrimR, lam, flux)
   use comvar
   implicit none

   real :: ConL(nvar), ConR(nvar)
   real :: PrimL(nvar), PrimR(nvar)
   real :: lam, flux(nvar)
   real :: fluxL(nvar), fluxR(nvar)
   real :: con_diff(nvar)
   
   con_diff = ConR - ConL
   
   call flux_x(ConL, PrimL, fluxL)
   call flux_x(ConR, PrimR, fluxR)
  
   !write(*,*) con_diff
   flux = 0.5*(fluxL + fluxR) - 0.5*lam*con_diff
   
 end subroutine numflux_x
 
 
 
 subroutine flux_x(co, primb,flux)
    use comvar
    implicit none
    
    real :: co(nvar), prim(nvar), flux(nvar), primb(nvar)
    real :: p_dash, c_bar
     
    c_bar  = sqrt(gamma*primb(2)/primb(1))  
    p_dash = c_bar**2/primb(5)*co(5)
    
    call cons2prims(co, prim, primb)
    write(*,*) prim(1), primb(1)
    flux(1) = co(2)
    flux(2) = co(2)*prim(2) + p_dash
    flux(3) = co(2)*prim(3)
    flux(4) = co(2)*prim(4)
    flux(5) = (prim(5) + primb(5))*co(2)
  
  end subroutine flux_x
