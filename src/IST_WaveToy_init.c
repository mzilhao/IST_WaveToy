
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define SMALL (1.e-8)

static CCTK_REAL gaussian(CCTK_REAL RR, CCTK_REAL R0, CCTK_REAL sigma)
{
  const CCTK_REAL R0pert2 = (RR - R0)*(RR - R0);

  return exp(-0.5 * R0pert2 / (sigma * sigma));
}

static CCTK_REAL gaussian_r(CCTK_REAL RR, CCTK_REAL R0, CCTK_REAL sigma)
{
  return -gaussian(RR, R0, sigma) * (RR-R0) / (sigma * sigma);
}


void IST_WaveToy_init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL tt = cctk_time;

#pragma omp parallel for collapse(3)
  for (CCTK_INT k = 0; k < cctk_lsh[2]; ++k) {
    for (CCTK_INT j = 0; j < cctk_lsh[1]; ++j) {
      for (CCTK_INT i = 0; i < cctk_lsh[0]; ++i) {

        const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

        const CCTK_REAL x1 = x[ijk];
        const CCTK_REAL y1 = y[ijk];
        const CCTK_REAL z1 = z[ijk];

        const CCTK_REAL RR2 = x1 * x1 + y1 * y1 + z1 * z1;
        const CCTK_REAL RR  = sqrt(RR2);

        const CCTK_REAL Rpt = RR + tt;
        const CCTK_REAL Rmt = RR - tt;

        if (RR > SMALL) {
          phi[ijk]  = 0.5 * (   Rpt * gaussian(Rpt,  R0pert_g, sigma_g)
                              + Rmt * gaussian(Rmt, -R0pert_g, sigma_g) ) / RR;

          // Kphi = \partial_t phi
          Kphi[ijk] = 0.5 * (   Rpt * gaussian_r(Rpt,  R0pert_g, sigma_g)
                              - Rmt * gaussian_r(Rmt, -R0pert_g, sigma_g)
                              + gaussian(Rpt,  R0pert_g, sigma_g)
                              - gaussian(Rmt, -R0pert_g, sigma_g) ) / RR;
        } else {
          phi[ijk]  = gaussian(tt, R0pert_g, sigma_g)
            * (1. + tt * (R0pert_g - tt) / (sigma_g*sigma_g));

          Kphi[ijk] = gaussian_r(tt, R0pert_g, sigma_g)
            * (1. + tt * (R0pert_g - tt) / (sigma_g*sigma_g))
            + gaussian(tt, R0pert_g, sigma_g)
            * (R0pert_g - 2. * tt) / (sigma_g*sigma_g);
        }

      } /* for i */
    }   /* for j */
  }     /* for k */

    return;
}
