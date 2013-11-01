#include "nordlund2013.h"


SEXP column_max_delta_beta(SEXP X, SEXP groups, SEXP feats){
    int i, j, g, g_nonzero, *n, *X_dim;
    double *means, min, max;
    SEXP dbeta;

    X = AS_NUMERIC(X);
    groups = AS_INTEGER(groups);
    feats = AS_LOGICAL(feats);
    X_dim = INTEGER(GET_DIM(X));

    // Find number of groups
    g = 0;
    for(i = 0; i < X_dim[0]; i++)
        if(g < INTEGER(groups)[i])
            g = INTEGER(groups)[i];
    n = (int*)malloc(g*sizeof(int));
    means = (double*)malloc(g*sizeof(double));

    // Calculate delta beta
    PROTECT(dbeta = NEW_NUMERIC(X_dim[1]));
    for(j = 0; j < X_dim[1]; j++){
        if(LOGICAL(feats)[j]){
            n[0] = 0;
            n[1] = 0;
            means[0] = 0;
            means[1] = 0;
            // Calculate group means
            for(i = 0; i < X_dim[0]; i++){
                if(INTEGER(groups)[i] > 0 && !R_IsNA(REAL(X)[i+j*X_dim[0]])){
                    n[INTEGER(groups)[i]-1]++;
                    means[INTEGER(groups)[i]-1] += REAL(X)[i+j*X_dim[0]];
                }
            }
            g_nonzero = 0;
            for(i = 0; i < g; i++){
                if(n[i] > 0){
                    means[i] /= n[i];
                    g_nonzero++;
                }
            }
            if(g_nonzero < 2){
                REAL(dbeta)[j] = NA_REAL;
            } else {
                for(i = 0; n[i] == 0; i++);
                min = means[i];
                max = means[i];
                for(i++; i < g; i++){
                    if(n[i] > 0){
                        if(means[i] < min) min = means[i];
                        else if(means[i] > max) max = means[i];
                    }
                }
                REAL(dbeta)[j] = max - min;
            }
        } else {
            REAL(dbeta)[j] = NA_REAL;
        }
    }

    free(n);
    free(means);
    UNPROTECT(1);
    return dbeta;
}


SEXP column_na_frac(SEXP X, SEXP thres, SEXP subset, SEXP feats){
	int i, j, k, count;
    int *X_dim;
	SEXP ret;

    X_dim = INTEGER(GET_DIM(X));
    subset = AS_LOGICAL(subset);
    feats = AS_LOGICAL(feats);

    k = 0;
    for(i = 0; i < X_dim[0]; i++)
        if(LOGICAL(subset)[i]) k++;
    k = (int)ceil(REAL(thres)[0]*k);

	PROTECT(ret = NEW_LOGICAL(X_dim[1]));
	for(j = 0; j < X_dim[1]; j++){
        if(LOGICAL(feats)[j]){
            LOGICAL(ret)[j] = true;
            count = 0;
            for(i = 0; i < X_dim[0]; i++){
                if(LOGICAL(subset)[i] && R_IsNA(REAL(X)[i+X_dim[0]*j])){
                    count++;
                    if(count >= k){
                        LOGICAL(ret)[j] = false;
                        break;
                    }
                }
            }
        } else {
            LOGICAL(ret)[j] = NA_LOGICAL;
        }
	}
	UNPROTECT(1);
	return ret;
}
