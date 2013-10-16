#ifndef _nordlund2014_NORDLUND2014_H
#define _nordlund2014_NORDLUND2014_H

#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>

extern "C" SEXP column_max_delta_beta(SEXP X, SEXP groups, SEXP feats);
extern "C" SEXP column_na_frac(SEXP X, SEXP thres, SEXP subset, SEXP feats);

#endif
