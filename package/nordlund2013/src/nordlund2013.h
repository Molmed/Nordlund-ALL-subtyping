#ifndef _nordlund2013_NORDLUND2013_H
#define _nordlund2013_NORDLUND2013_H

#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>

extern "C" SEXP column_max_delta_beta(SEXP X, SEXP groups, SEXP feats);
extern "C" SEXP column_na_frac(SEXP X, SEXP thres, SEXP subset, SEXP feats);

#endif
