#include "header.h"
#include "rates.h"

void eval_spec_rates (const double * __restrict__ fwd_rates, const double * __restrict__ rev_rates, const double * __restrict__ pres_mod, double * __restrict__ sp_rates, double * __restrict__ dy_N) {
  //rxn 0
  //sp 1
  sp_rates[1] = -(fwd_rates[0] - rev_rates[0]);
  //sp 3
  sp_rates[3] = -(fwd_rates[0] - rev_rates[0]);
  //sp 4
  sp_rates[4] = (fwd_rates[0] - rev_rates[0]);
  //sp 5
  sp_rates[5] = (fwd_rates[0] - rev_rates[0]);

  //rxn 1
  //sp 0
  sp_rates[0] = -(fwd_rates[1] - rev_rates[1]);
  //sp 3
  sp_rates[3] += (fwd_rates[1] - rev_rates[1]);
  //sp 4
  sp_rates[4] -= (fwd_rates[1] - rev_rates[1]);
  //sp 5
  sp_rates[5] += (fwd_rates[1] - rev_rates[1]);

  //rxn 2
  //sp 0
  sp_rates[0] -= (fwd_rates[2] - rev_rates[2]);
  //sp 2
  sp_rates[2] = (fwd_rates[2] - rev_rates[2]);
  //sp 3
  sp_rates[3] += (fwd_rates[2] - rev_rates[2]);
  //sp 5
  sp_rates[5] -= (fwd_rates[2] - rev_rates[2]);

  //rxn 3
  //sp 2
  sp_rates[2] -= (fwd_rates[3] - rev_rates[3]);
  //sp 4
  sp_rates[4] -= (fwd_rates[3] - rev_rates[3]);
  //sp 5
  sp_rates[5] += 2.0 * (fwd_rates[3] - rev_rates[3]);

  //rxn 4
  //sp 0
  sp_rates[0] -= (fwd_rates[4] - rev_rates[4]) * pres_mod[0];
  //sp 3
  sp_rates[3] += 2.0 * (fwd_rates[4] - rev_rates[4]) * pres_mod[0];

  //rxn 5
  //sp 0
  sp_rates[0] -= (fwd_rates[5] - rev_rates[5]);
  //sp 3
  sp_rates[3] += 2.0 * (fwd_rates[5] - rev_rates[5]);

  //rxn 6
  //sp 0
  sp_rates[0] -= (fwd_rates[6] - rev_rates[6]);
  //sp 3
  sp_rates[3] += 2.0 * (fwd_rates[6] - rev_rates[6]);

  //rxn 7
  //sp 1
  sp_rates[1] += (fwd_rates[7] - rev_rates[7]) * pres_mod[1];
  //sp 4
  sp_rates[4] -= 2.0 * (fwd_rates[7] - rev_rates[7]) * pres_mod[1];

  //rxn 8
  //sp 1
  sp_rates[1] += (fwd_rates[8] - rev_rates[8]);
  //sp 4
  sp_rates[4] -= 2.0 * (fwd_rates[8] - rev_rates[8]);

  //rxn 9
  //sp 4
  sp_rates[4] -= 2.0 * (fwd_rates[9] - rev_rates[9]);
  //sp 1
  sp_rates[1] += (fwd_rates[9] - rev_rates[9]);

  //rxn 10
  //sp 3
  sp_rates[3] -= (fwd_rates[10] - rev_rates[10]) * pres_mod[2];
  //sp 4
  sp_rates[4] -= (fwd_rates[10] - rev_rates[10]) * pres_mod[2];
  //sp 5
  sp_rates[5] += (fwd_rates[10] - rev_rates[10]) * pres_mod[2];

  //rxn 11
  //sp 2
  sp_rates[2] += (fwd_rates[11] - rev_rates[11]) * pres_mod[3];
  //sp 3
  sp_rates[3] -= (fwd_rates[11] - rev_rates[11]) * pres_mod[3];
  //sp 5
  sp_rates[5] -= (fwd_rates[11] - rev_rates[11]) * pres_mod[3];

  //rxn 12
  //sp 1
  sp_rates[1] -= (fwd_rates[12] - rev_rates[12]) * pres_mod[4];
  //sp 3
  sp_rates[3] -= (fwd_rates[12] - rev_rates[12]) * pres_mod[4];
  //sp 6
  sp_rates[6] = (fwd_rates[12] - rev_rates[12]) * pres_mod[4];

  //rxn 13
  //sp 0
  sp_rates[0] += (fwd_rates[13] - rev_rates[13]);
  //sp 1
  sp_rates[1] += (fwd_rates[13] - rev_rates[13]);
  //sp 3
  sp_rates[3] -= (fwd_rates[13] - rev_rates[13]);
  //sp 6
  sp_rates[6] -= (fwd_rates[13] - rev_rates[13]);

  //rxn 14
  //sp 3
  sp_rates[3] -= (fwd_rates[14] - rev_rates[14]);
  //sp 5
  sp_rates[5] += 2.0 * (fwd_rates[14] - rev_rates[14]);
  //sp 6
  sp_rates[6] -= (fwd_rates[14] - rev_rates[14]);

  //rxn 15
  //sp 1
  sp_rates[1] += (fwd_rates[15] - rev_rates[15]);
  //sp 4
  sp_rates[4] -= (fwd_rates[15] - rev_rates[15]);
  //sp 5
  sp_rates[5] += (fwd_rates[15] - rev_rates[15]);
  //sp 6
  sp_rates[6] -= (fwd_rates[15] - rev_rates[15]);

  //rxn 16
  //sp 1
  sp_rates[1] += (fwd_rates[16] - rev_rates[16]);
  //sp 2
  sp_rates[2] += (fwd_rates[16] - rev_rates[16]);
  //sp 5
  sp_rates[5] -= (fwd_rates[16] - rev_rates[16]);
  //sp 6
  sp_rates[6] -= (fwd_rates[16] - rev_rates[16]);

  //rxn 17
  //sp 1
  sp_rates[1] += (fwd_rates[17] - rev_rates[17]);
  //sp 6
  sp_rates[6] -= 2.0 * (fwd_rates[17] - rev_rates[17]);
  //sp 7
  sp_rates[7] = (fwd_rates[17] - rev_rates[17]);

  //rxn 18
  //sp 1
  sp_rates[1] += (fwd_rates[18] - rev_rates[18]);
  //sp 6
  sp_rates[6] -= 2.0 * (fwd_rates[18] - rev_rates[18]);
  //sp 7
  sp_rates[7] += (fwd_rates[18] - rev_rates[18]);

  //rxn 19
  //sp 5
  sp_rates[5] += 2.0 * (fwd_rates[19] - rev_rates[19]) * pres_mod[5];
  //sp 7
  sp_rates[7] -= (fwd_rates[19] - rev_rates[19]) * pres_mod[5];

  //rxn 20
  //sp 2
  sp_rates[2] += (fwd_rates[20] - rev_rates[20]);
  //sp 3
  sp_rates[3] -= (fwd_rates[20] - rev_rates[20]);
  //sp 5
  sp_rates[5] += (fwd_rates[20] - rev_rates[20]);
  //sp 7
  sp_rates[7] -= (fwd_rates[20] - rev_rates[20]);

  //rxn 21
  //sp 0
  sp_rates[0] += (fwd_rates[21] - rev_rates[21]);
  //sp 3
  sp_rates[3] -= (fwd_rates[21] - rev_rates[21]);
  //sp 6
  sp_rates[6] += (fwd_rates[21] - rev_rates[21]);
  //sp 7
  sp_rates[7] -= (fwd_rates[21] - rev_rates[21]);

  //rxn 22
  //sp 4
  sp_rates[4] -= (fwd_rates[22] - rev_rates[22]);
  //sp 5
  sp_rates[5] += (fwd_rates[22] - rev_rates[22]);
  //sp 6
  sp_rates[6] += (fwd_rates[22] - rev_rates[22]);
  //sp 7
  sp_rates[7] -= (fwd_rates[22] - rev_rates[22]);

  //rxn 23
  //sp 2
  sp_rates[2] += (fwd_rates[23] - rev_rates[23]);
  //sp 5
  sp_rates[5] -= (fwd_rates[23] - rev_rates[23]);
  //sp 6
  sp_rates[6] += (fwd_rates[23] - rev_rates[23]);
  //sp 7
  sp_rates[7] -= (fwd_rates[23] - rev_rates[23]);

  //rxn 24
  //sp 2
  sp_rates[2] += (fwd_rates[24] - rev_rates[24]);
  //sp 5
  sp_rates[5] -= (fwd_rates[24] - rev_rates[24]);
  //sp 6
  sp_rates[6] += (fwd_rates[24] - rev_rates[24]);
  //sp 7
  sp_rates[7] -= (fwd_rates[24] - rev_rates[24]);

  //sp 8
  sp_rates[8] = 0.0;
  //sp 10
  sp_rates[9] = 0.0;
  //sp 11
  sp_rates[10] = 0.0;
  //sp 12
  sp_rates[11] = 0.0;
  //sp 9
  (*dy_N) = 0.0;
} // end eval_spec_rates

