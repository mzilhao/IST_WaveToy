
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine IST_WaveToy_calc_rhs( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL                 :: lphi, lKphi
  CCTK_REAL, dimension(3,3) :: d2_lphi
  CCTK_REAL                 :: rhs_lphi, rhs_lKphi
  CCTK_REAL                 :: odx, ody, odz
  CCTK_INT,  dimension(6)   :: bndsize, is_ghostbnd, is_symbnd, is_physbnd
  CCTK_INT,  dimension(3)   :: istart(3), iend(3)
  CCTK_INT                  :: i, j, k
  CCTK_INT                  :: ierr

  ! SBP
  CCTK_REAL, dimension(:,:), allocatable :: qx2, qy2, qz2
  CCTK_INT,  dimension(:),   allocatable :: imin2, imax2, jmin2, &
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

  ! write(*,*) "bndsize     = ", bndsize
  ! write(*,*) "is_ghostbnd = ", is_ghostbnd
  ! write(*,*) "is_symbnd   = ", is_symbnd
  ! write(*,*) "is_physbnd  = ", is_physbnd
  ! write(*,*) "istart      = ", istart
  ! write(*,*) "iend        = ", iend

  odx = 1.0d0 / CCTK_DELTA_SPACE(1)
  ody = 1.0d0 / CCTK_DELTA_SPACE(2)
  odz = 1.0d0 / CCTK_DELTA_SPACE(3)

  ! loop interior points

  do k = istart(3), iend(3)
  do j = istart(2), iend(2)
  do i = istart(1), iend(1)

     !------------ get local variables ---------
     lphi    = phi(i,j,k)
     lKphi   = Kphi(i,j,k)

     !------------ 2nd derivatives ------------
     d2_lphi(1,1) = odx * odx * sum( qx2(imin2(i):imax2(i),i) * phi(imin2(i):imax2(i),j,k) )
     d2_lphi(2,2) = ody * ody * sum( qy2(jmin2(j):jmax2(j),j) * phi(i,jmin2(j):jmax2(j),k) )
     d2_lphi(3,3) = odz * odz * sum( qz2(kmin2(k):kmax2(k),k) * phi(i,j,kmin2(k):kmax2(k)) )

     !------------ source terms -----------------
     rhs_lphi   = lKphi
     rhs_lKphi  =  d2_lphi(1,1) + d2_lphi(2,2) + d2_lphi(3,3)

     !------------ write to grid functions ------
     rhs_phi(i,j,k)  = rhs_lphi
     rhs_Kphi(i,j,k) = rhs_lKphi

     ! if (i == 4 .and. j == 4) then
     !    write(*,*) "bulk"
     !    write(*,*) "i,j,k = ", i,j,k
     !    write (*,*) ite(*,*) "x,y,z = ", x(i,j,k), y(i,j,k), z(i,j,k)
     !    write(*,*) ""
     ! end if

  end do
  end do
  end do

  deallocate( qx2, qy2, qz2, imin2, imax2, jmin2, jmax2, kmin2, kmax2)

end subroutine IST_WaveToy_calc_rhs
!
!=============================================================================
!
subroutine IST_WaveToy_calc_rhs_bdry( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_LOOP3_INTBND_DECLARE(bdry)

  CCTK_INT                :: i, j, k, ni, nj, nk
  CCTK_INT                :: reflevel
  CCTK_REAL               :: odx, ody, odz
  CCTK_REAL               :: lphi, lKphi
  CCTK_REAL, dimension(3) :: d1_lKphi(3)
  CCTK_REAL               :: rr, vec(3)

  ! SBP
  CCTK_REAL, dimension(:,:), allocatable :: qx, qy, qz
  CCTK_INT,  dimension(:),   allocatable :: imin, imax, jmin, &
                                            jmax, kmin, kmax

  allocate(qx(cctk_lsh(1),cctk_lsh(1)), &
           qy(cctk_lsh(2),cctk_lsh(2)), &
           qz(cctk_lsh(3),cctk_lsh(3)), &
           imin(cctk_lsh(1)), imax(cctk_lsh(1)), &
           jmin(cctk_lsh(2)), jmax(cctk_lsh(2)), &
           kmin(cctk_lsh(3)), kmax(cctk_lsh(3)))

  call Diff_coeff( cctkGH, 0, cctk_lsh(1), imin, imax, qx, -1 )
  call Diff_coeff( cctkGH, 1, cctk_lsh(2), jmin, jmax, qy, -1 )
  call Diff_coeff( cctkGH, 2, cctk_lsh(3), kmin, kmax, qz, -1 )

  reflevel = GetRefinementLevel(cctkGH)
  ! write(*,*) "reflevel = ", reflevel
  ! apply only on the coarsest level
  if (reflevel /= 0) return

  odx = 1.0d0 / (CCTK_DELTA_SPACE(1))
  ody = 1.0d0 / (CCTK_DELTA_SPACE(2))
  odz = 1.0d0 / (CCTK_DELTA_SPACE(3))

  ! loop boundary points

  CCTK_LOOP3_INTBND(bdry, i, j, k, ni, nj, nk)

     lphi   = phi(i,j,k)
     lKphi  = Kphi(i,j,k)

     rr     = sqrt( x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2 )
     vec(:) = (/ x(i,j,k) / rr, y(i,j,k) / rr, z(i,j,k) / rr /)

     d1_lKphi(1) = odx * sum( qx(imin(i):imax(i),i) * Kphi(imin(i):imax(i),j,k) )
     d1_lKphi(2) = ody * sum( qy(jmin(j):jmax(j),j) * Kphi(i,jmin(j):jmax(j),k) )
     d1_lKphi(3) = odz * sum( qz(kmin(k):kmax(k),k) * Kphi(i,j,kmin(k):kmax(k)) )

     rhs_Kphi(i,j,k) = -vec(1)*d1_lKphi(1) - vec(2)*d1_lKphi(2) - vec(3)*d1_lKphi(3)  &
                     - lKphi / rr

     rhs_phi(i,j,k)  = lKphi

     ! if (j == 4) then
     !    write(*,*) "bdry"
     !    write(*,*) "i,j,k    = ", i,j,k
     !    write(*,*) "x,y,z    = ", x(i,j,k), y(i,j,k), z(i,j,k)
     !    write(*,*) ""
     ! end if

  CCTK_ENDLOOP3_INTBND(bdry)

  deallocate( qx, qy, qz, imin, imax, jmin, jmax, kmin, kmax)

end subroutine IST_WaveToy_calc_rhs_bdry
!
!=============================================================================
!
