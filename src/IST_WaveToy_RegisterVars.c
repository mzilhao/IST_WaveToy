
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void IST_WaveToy_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, rhs;

  // register evolution and rhs gridfunction groups with MoL

  /* phi and rhs_phi */
  group = CCTK_GroupIndex("IST_WaveToy::phi");
  rhs   = CCTK_GroupIndex("IST_WaveToy::rhs_phi");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Kphi and rhs_Kphi */
  group = CCTK_GroupIndex("IST_WaveToy::Kphi");
  rhs   = CCTK_GroupIndex("IST_WaveToy::rhs_Kphi");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  if (ierr) CCTK_ERROR("Problems registering with MoL");

}
