#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine IST_WaveToy_InitSymBound( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "IST_WaveToy::phi" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "IST_WaveToy::Kphi" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "IST_WaveToy::rhs_phi" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "IST_WaveToy::rhs_Kphi" )

end subroutine IST_WaveToy_InitSymBound
