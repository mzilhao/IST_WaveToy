
#include "cctk.h"
#include "cctk_Arguments.h"

subroutine IST_WaveToy_zero_rhs( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS

  rhs_phi    = 0
  rhs_Kphi   = 0

end subroutine IST_WaveToy_zero_rhs
