! ============================================
subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux,actualstep)
! ============================================
!
! # called before each call to step
! # use to set time-dependent aux arrays or perform other tasks.
!
! This particular routine sets negative values of q(1,i,j) to zero,
! as well as the corresponding q(m,i,j) for m=1,meqn.
! This is for problems where q(1,i,j) is a depth.
! This should occur only because of rounding error.
!
! Also calls movetopo if topography might be moving.

    use geoclaw_module, only: dry_tolerance
    use geoclaw_module, only: grav
    use geoclaw_module, only: spherical_distance
    use geoclaw_module, only: coordinate_system, earth_radius, deg2rad
    use topo_module, only: num_dtopo,topotime
    use topo_module, only: aux_finalized
    use topo_module, only: xlowdtopo,xhidtopo,ylowdtopo,yhidtopo

    use amr_module, only: xlowdomain => xlower
    use amr_module, only: ylowdomain => ylower
    use amr_module, only: xhidomain => xupper
    use amr_module, only: yhidomain => yupper
    use amr_module, only: xperdom,yperdom,spheredom,NEEDS_TO_BE_SET
    use vp2_module

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: meqn
    integer, intent(inout) :: mbc,mx,my,maux
    real(kind=8), intent(inout) :: xlower, ylower, dx, dy, t, dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    logical, intent (in) :: actualstep

    ! Local storage
    integer :: index,i,j,k
    real(kind=8) :: h,u,v
    real(kind=8) :: xm,xc,xp,ym,yc,yp,dx_meters,dy_meters

    ! Bing model
    real(kind=8) :: gxmod,gymod,upmax
    real(kind=8) :: up,vp,tau,g,tau_y,tau1,tau2

    ! Retrogression
    real(kind=8) :: xr,ret_vel

    ! check for h < 0 and reset to zero
    ! check for h < drytolerance
    ! set hu = hv = 0 in all these cells
    forall(i=1-mbc:mx+mbc, j=1-mbc:my+mbc,q(1,i,j) < dry_tolerance)
        q(1,i,j) = max(q(1,i,j),0.d0)
        q(2:meqn,i,j) = 0.d0
    end forall

    ! Check for NaNs in the solution
    call check4nans(meqn,mbc,mx,my,q,t,1)

    if (aux_finalized < 2) then
        ! topo arrays might have been updated by dtopo more recently than
        ! aux arrays were set unless at least 1 step taken on all levels
        aux(1,:,:) = NEEDS_TO_BE_SET ! new system checks this val before setting
        call setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
        endif

    ! Determine the stopping case
    ! based on the earth pressure and the yield stress

    g=grav*(1.d0-rho_a/rho_s)
    aux(6:7,:,:) = 0.d0

    if (remolding) then
       do i=1,mx
          do j=1,my
             if (q(1,i,j)>0.d0) then
                 q(6,i,j) = max(0.d0,q(6,i,j))
             else 
                 q(6,i,j) = 0.d0
             endif
          enddo
       enddo
    endif

    do i=1,mx

       xm = xlower + (i - 1.d0) * dx
       xc = xlower + (i - 0.5d0) * dx
       xp = xlower + i * dx

       do j=1,my

          ym = ylower + (j - 1.d0) * dy
          yc = ylower + (j - 0.5d0) * dy
          yp = ylower + j * dy

          h  = q(1,i,j)

          if (h > dry_tolerance) then

             up = q(2,i,j)
             vp = q(3,i,j)

             tau = 0.d0

             tau_y=tauy_r+(tauy_i-tauy_r)/exp(q(6,i,j)*remold_coeff)
               
             if (coordinate_system == 2) then
                 ! Convert distance in lat-long to meters
                 dx_meters = spherical_distance(xp,yc,xm,yc)
                 dy_meters = spherical_distance(xc,yp,xc,ym)
             else
                 dx_meters = dx
                 dy_meters = dy
             endif

             if (aux(1,i,j)+q(1,i,j)<0.d0) then
                aux(4,i,j) = grav*(1.d0-rho_a/rho_s) &
                   *dcos(atan((aux(1,i+1,j)-aux(1,i-1,j))/2.d0/dx_meters))
                aux(5,i,j) = grav*(1.d0-rho_a/rho_s) &
                   *dcos(atan((aux(1,i,j+1)-aux(1,i,j-1))/2.d0/dy_meters))
             else
                aux(4,i,j) = grav* &
                   dcos(atan((aux(1,i+1,j)-aux(1,i-1,j))/2.d0/dx_meters))
                aux(5,i,j) = grav* &
                   dcos(atan((aux(1,i,j+1)-aux(1,i,j-1))/2.d0/dy_meters))
             endif

             gxmod = aux(4,i,j)
             gymod = aux(5,i,j)

             if (dsqrt(up**2+vp**2)<1.d-4) then

                tau1 = (q(1,i+1,j)-q(1,i-1,j)+aux(1,i+1,j)-aux(1,i-1,j)) &
                      /dx_meters/2.d0*gxmod

                if (abs(tau1)<tau_y/rho_s/h) then
                   aux(6,i,j) = 1.d0
                endif

                tau2 = (q(1,i,j+1)-q(1,i,j-1)+aux(1,i,j+1)-aux(1,i,j-1)) &
                      /dy_meters/2.d0*gymod

                if (abs(tau2)<tau_y/rho_s/h) then
                   aux(7,i,j) = 1.d0
                endif

                tau = dsqrt(tau1**2 +tau2**2)

                if (tau < tau_y/rho_s/h) then
                   aux(6,i,j) = 1.d0
                   aux(7,i,j) = 1.d0
                endif

             endif
          else
             aux(6,i,j) = 1.d0
             aux(7,i,j) = 1.d0
          endif

       end do

    end do

    call checkval_bing(mx,my,meqn,mbc,q,alpha_1)
	
end subroutine b4step2
