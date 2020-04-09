#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine IST_WaveToy_testID( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS


  CCTK_INT  :: i,j,k
  CCTK_INT  :: ierr

  ! Second derivatives
  CCTK_REAL :: d2_lphi(3,3)

  CCTK_REAL :: RR2, RR, x1, y1, z1
  CCTK_REAL :: Rpt, Rmt
  CCTK_REAL :: IST_WaveToy_gaussian

  CCTK_REAL :: odx, ody, odz

  CCTK_INT, dimension(6) :: bndsize, is_ghostbnd, is_symbnd, is_physbnd
  CCTK_INT, dimension(3) :: istart(3), iend(3)

  CCTK_REAL :: tt

  ! SBP
  CCTK_REAL, dimension(:,:), allocatable :: qx2, qy2, qz2
  CCTK_INT,  dimension(:), allocatable   :: imin2, imax2, jmin2, &
                                            jmax2, kmin2, kmax2

  allocate(qx2(cctk_lsh(1),cctk_lsh(1)), &
           qy2(cctk_lsh(2),cctk_lsh(2)), &
           qz2(cctk_lsh(3),cctk_lsh(3)), &
           imin2(cctk_lsh(1)), imax2(cctk_lsh(1)), &
           jmin2(cctk_lsh(2)), jmax2(cctk_lsh(2)), &
           kmin2(cctk_lsh(3)), kmax2(cctk_lsh(3)))

  call Diff2_coeff( cctkGH, 0, cctk_lsh(1), imin2, imax2, qx2, -1 )
  call Diff2_coeff( cctkGH, 1, cctk_lsh(2), jmin2, jmax2, qy2, -1 )
  call Diff2_coeff( cctkGH, 2, cctk_lsh(3), kmin2, kmax2, qz2, -1 )


  ierr = GetBoundarySizesAndTypes( cctkGH, 6, bndsize, is_ghostbnd, is_symbnd, is_physbnd )

  ! set the limits for looping only in the interior points
  istart = 1 + bndsize(1:5:2)
  iend   = cctk_lsh(:) - bndsize(2:6:2)

  write(*,*) "bndsize     = ", bndsize
  write(*,*) "is_ghostbnd = ", is_ghostbnd
  write(*,*) "is_symbnd   = ", is_symbnd
  write(*,*) "is_physbnd  = ", is_physbnd
  write(*,*) "istart      = ", istart
  write(*,*) "iend        = ", iend


  odx = 1.0d0 / CCTK_DELTA_SPACE(1)
  ody = 1.0d0 / CCTK_DELTA_SPACE(2)
  odz = 1.0d0 / CCTK_DELTA_SPACE(3)


  tt = cctk_time

  do k = 1, cctk_lsh(3)
  do j = 1, cctk_lsh(2)
  do i = 1, cctk_lsh(1)

     x1  = x(i,j,k)
     y1  = y(i,j,k)
     z1  = z(i,j,k)

     RR2 = x1*x1 + y1*y1 + z1*z1
     RR  = sqrt(RR2)

     Rpt = RR + tt
     Rmt = RR - tt

     if (RR > 1.0d-8) then
        phi(i,j,k) =  0.5 * (   Rpt * IST_WaveToy_gaussian(Rpt,  R0pert_g, sigma_g)    &
                              + Rmt * IST_WaveToy_gaussian(Rmt, -R0pert_g, sigma_g) ) / RR
     else
        phi(i,j,k) = IST_WaveToy_gaussian(tt, R0pert_g, sigma_g)               &
             * (1.0d0 + tt * (R0pert_g - tt) / (sigma_g*sigma_g))
     end if

  end do
  end do
  end do


  ! loop only interior points

  do k = istart(3), iend(3)
  do j = istart(2), iend(2)
  do i = istart(1), iend(1)

     !------------- 2nd derivatives -----------
     d2_lphi(1,1) = odx * odx * sum( qx2(imin2(i):imax2(i),i) * phi(imin2(i):imax2(i),j,k) )
     d2_lphi(2,2) = ody * ody * sum( qy2(jmin2(j):jmax2(j),j) * phi(i,jmin2(j):jmax2(j),k) )
     d2_lphi(3,3) = odz * odz * sum( qz2(kmin2(k):kmax2(k),k) * phi(i,j,kmin2(k):kmax2(k)) )


     if(i == 5 .and. j == 5) then
        write(*,*) "x,y,z = ", x(i,j,k), y(i,j,k), z(i,j,k)
        write(*,*) "kmin2, k, kmax2 = ", kmin2(k), k, kmax2(k)
        write(*,*) "qz2 ", qz2(kmin2(k):kmax2(k),k)
        write(*,*) ""
     end if


     Kphi(i,j,k) = 666
  end do
  end do
  end do


  deallocate( qx2, qy2, qz2, &
              imin2, imax2, jmin2, jmax2, kmin2, kmax2)

end subroutine IST_WaveToy_testID



CCTK_REAL function IST_WaveToy_gaussian(RR, R0, sigma)
    implicit none
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL :: RR
    CCTK_REAL :: R0
    CCTK_REAL :: sigma

    IST_WaveToy_gaussian = exp(-0.5d0 * (RR - R0)**2 / (sigma * sigma))
    return
end function IST_WaveToy_gaussian
