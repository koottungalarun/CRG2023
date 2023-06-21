subroutine numflux_z(ConL, ConR, PrimL, PrimR, lam, flux)
   use comvar
   implicit none

   real :: ConL(nvar), ConR(nvar)
   real :: PrimL(nvar), PrimR(nvar)
   real :: lam, flux(nvar)
   real :: fluxL(nvar), fluxR(nvar)
   real :: con_diff(nvar)
   
   con_diff = ConR - ConL
   
   call flux_z(ConL, PrimL, fluxL)
   call flux_z(ConR, PrimR, fluxR)
   
   flux = 0.5*(fluxL + fluxR) - 0.5*lam*con_diff
   
 end subroutine numflux_z
 
 
 subroutine flux_z(co, primb, flux)
    use comvar
    implicit none
    
    real :: co(nvar), primb(nvar), flux(nvar)
    real :: prim(nvar)
  
    
    call cons2prims(co, prim, primb)
    flux(1) = 0.0
    flux(2) = co(4)*prim(2) 
    flux(3) = co(4)*prim(3) 
    flux(4) = co(4)*prim(4) 
    flux(5) = prim(5)*co(4)
  
 end subroutine flux_z
