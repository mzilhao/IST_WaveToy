#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine IST_WaveToy_Boundaries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT ierr

  CCTK_INT, parameter :: one = 1

  ! The actual boundary conditions are being handled in the evolution file. Here
  ! we register all BCs as 'none', which enforces all the symmetry BCs.

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "IST_WaveToy::phi", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for IST_WaveToy::phi!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "IST_WaveToy::Kphi", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for IST_WaveToy::Kphi!")

end subroutine IST_WaveToy_Boundaries
