subroutine numflux_y(ConL, ConR, PrimL, PrimR, lam, flux)
   use comvar
   implicit none

   real :: ConL(nvar), ConR(nvar)
   real :: PrimL(nvar), PrimR(nvar)

   real :: lam, flux(nvar)
   real :: fluxL(nvar), fluxR(nvar)
   real :: con_diff(nvar)
   
   con_diff = ConR - ConL
   
   call flux_y(ConL, PrimL, fluxL)
   call flux_y(ConR, PrimR, fluxR)
  
   
   
   flux = 0.5*(fluxL + fluxR) - 0.5*lam*con_diff
   
 end subroutine numflux_y
 
 
  subroutine flux_y(co, primb, flux)
    use comvar
    implicit none
    
    real :: co(nvar), prim(nvar), primb(nvar),flux(nvar)
    real :: p_dash, c_bar
    
    c_bar  = sqrt(gamma*Pre_bar/primb(1))  ! need to take car of p_bar
    p_dash = c_bar**2/primb(5)*co(5)
    
    call cons2prims(co, prim, primb)
    flux(1) = co(3)
    flux(2) = co(3)*prim(2) 
    flux(3) = co(3)*prim(3) + p_dash
    flux(4) = co(3)*prim(4) 
    flux(5) = (primb(5) + prim(5))*co(3)
  
  end subroutine flux_y
