subroutine setprob()

    use qinit_module, only: qinit_style
    use vp2_module

    implicit none

    integer :: iunit
    character(len=25) fname

    iunit = 7
    fname = 'setprob.data'
!   # open the unit with new routine from Clawpack 4.4 to skip over
!   # comment lines starting with #:
    call opendatafile(iunit, fname)

    read(7,*) rho_a
    read(7,*) rho_s
    read(7,*) n_param
    read(7,*) gamma_r
    read(7,*) c_mass
    read(7,*) hydrodrag
    read(7,*) cF_hyd 
    read(7,*) cP_hyd    
    read(7,*) remolding 
    read(7,*) tauy_i 
    read(7,*) tauy_r 
    read(7,*) remold_coeff 
    read(7,*) qinit_style

    visc_kin = visc_dyn/rho_s
    alpha_1 = (1.d0/n_param+1.d0)/(1.d0/n_param+2.d0) 
    alpha_2 = 1.d0 - 2.d0/(1.d0/n_param + 2.d0) + 1.d0/(2.d0/n_param+3.d0) 
    beta = (1.d0+1.d0/n_param)**n_param 
    cm_coeff = 1.d0+c_mass*rho_a/rho_s 

end subroutine setprob
