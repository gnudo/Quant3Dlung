#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
                 int nrhs, const mxArray *prhs[]) /* Input variables */
{

	/*
     *   Define the macros for the output and input arguments
     */

  #define img_out   plhs[0]   /*  Macro for the output sinogram array  */
  #define img_in    prhs[0]   /*  Macro for the input image array  */
  #define iteration prhs[1]   /*  Macro for the number of iterations  */
  
  double *I_in;          /*  Operative double pointer to the input image */
  double *I_out;         /*  Operative double pointer to the output image */
  double *n_iter;        /*  Operative double the number of iterations */
  int M, N, m, n, k, l, ll, counter;
  
  /* Establish pointer connections */
  I_in = mxGetPr(img_in);
  n_iter = mxGetPr(iteration);
  l = (int) mxGetScalar(iteration);
  
  /* Get Matrix dimensions */
  M = mxGetM(img_in);
  N = mxGetN(img_in);
  
  /* Create output matrix */
  /* img_out = mxCreateNumericMatrix( M , N , mxUINT8_CLASS , mxREAL ); */
  img_out = mxCreateDoubleMatrix(M, N, mxREAL); /* Create the output matrix */
  I_out = mxGetPr(img_out); /* Get the pointer to the data of B */
  
  /* Assign the Input to Output */
    for(n = 0; n < N; n++) {
        for(m = 0; m < M; m++) {
            I_out[m + M*n] = I_in[m + M*n];
        }
    }
   
  /* Run first filter w/o Loop */
  for (ll = 0; ll < l; ll++) {
	for (n = 1; n < N-1; n++) {
		for(m = 1; m < M-1; m++) {
			if ( (int) I_in[m + m*N] == 1) {
				counter = 0;
				for (k = 0; k < 3; k++) {
					counter = counter + (int) I_in[m-1+k + M*(n-1)];
					counter = counter + (int) I_in[m-1+k + M*n];
					counter = counter + (int) I_in[m-1+k + M*(n+1)];
				}
				if (counter <3) {
					I_out[m + M*n] = 0.0;
				}
			}
		}
    }
   
  /* Assign the Input to Output */
    for(n = 0; n < N; n++) {
        for(m = 0; m < M; m++) {
            I_in[m + M*n] = I_out[m + M*n];
        }
    }
  }
   
  /* Run second part */
  for (ll = 0; ll < l; ll++) {
	for (n = 1; n < N-1; n++) {
		for(m = 1; m < M-1; m++) {
			if ( (int) I_in[m + M*n] == 1) {
				counter = 0;
                /* Struct El. 1 */
				for (k = 0; k < 2; k++) {
					counter = counter + (int) I_in[m-1+k + M*(n-1)];
					counter = counter + (int) I_in[m-1+k + M*n];
					counter = counter + (int) I_in[m-1+k + M*(n+1)];
				}
				if (counter <2) {
					I_out[m + M*n] = 0.0;
                    break;
				}
                counter = 0;
                /* Struct El. 2 */
                for (k = 0; k < 2; k++) {
					counter = counter + (int) I_in[m+k + M*(n-1)];
					counter = counter + (int) I_in[m+k + M*n];
					counter = counter + (int) I_in[m+k + M*(n+1)];
				}
				if (counter <2) {
					I_out[m + M*n] = 0.0;
                    break;
				}
                counter = 0;
                /* Struct El. 3 */
                for (k = 0; k < 3; k++) {
					counter = counter + (int) I_in[m-1+k + M*(n-1)];
					counter = counter + (int) I_in[m-1+k + M*n];
				}
				if (counter <2) {
					I_out[m + M*n] = 0.0;
                    break;
				}
                counter = 0;
                /* Struct El. 4 */
                for (k = 0; k < 3; k++) {
					counter = counter + (int) I_in[m-1+k + M*n];
					counter = counter + (int) I_in[m-1+k + M*(n+1)];
				}
				if (counter <2) {
					I_out[m + M*n] = 0.0;
                    break;
				}
			}
		}
	}
   
  /* Assign the Input to Output */
    for(n = 0; n < N; n++) {
        for(m = 0; m < M; m++) {
            I_in[m + M*n] = I_out[m + M*n];
        }
    }
  }
 
  /* mexPrintf("Hello, world!\n"); Do something interesting */
  return;
}
