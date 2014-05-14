//#!/bin/bash
/*/../bin/ls > /dev/null
# BEGIN BASH SCRIPT
ME=`basename $0`
MEBASE=`echo $ME | cut -d. -f 1`
TMP=/tmp/$MEBASE
printf "//" | cat - $0 > $TMP.cc

CFLAGS='-I/usr/include/eigen3 -std=c++11 -fPIC -O3 -march=native -flto'
LDFLAGS='-O3 -march=native -flto'

g++ -c $TMP.cc -o $TMP.o -I/opt/matlab/extern/include $CFLAGS
g++ -shared -Wl,-soname,$MEBASE.mexa64 -Wl,--no-undefined -Wl,--as-needed\
  -o $MEBASE.mexa64 $TMP.o $LDFLAGS -L/opt/matlab/bin/glnxa64 \
  -lmx -lmex -lmat

rm -f $TMP.cc $TMP.o
# END BASH SCRIPT
exit
*/

/*
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include <mex.h>

#include <Eigen/Dense>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  using namespace Eigen;

  if(nrhs < 1 || nlhs != 1) {
    mexWarnMsgTxt("Minimum one input and one output parameter required");
    return;
  }

  //Parse opts
  bool fast2D = false;

  if(nrhs > 1) {
    const mxArray *opt = prhs[1];

    const mxArray *opt_fast2D = mxGetField(opt, 0, "fast2D");
    if(opt_fast2D)
      fast2D = mxGetScalar(opt_fast2D);
  }

  const mwSize* dims = mxGetDimensions(prhs[0]);

  const mwSize m = dims[0], n = dims[1];

  if(n%m != 0)
    mexErrMsgTxt("Invalid input (inconsistent dimensions).");

  mxArray *R        = mxCreateNumericMatrix(m, n,
                                         mxDOUBLE_CLASS, mxREAL);

  const double *data_in  = mxGetPr(prhs[0]);
        double *data_out = mxGetPr(R);

  const mwSize nt = n / m;

  if(m == 2) {

    if(fast2D) {
    //special case closed form solution for 2x2

    for(mwSize t = 0; t < nt; ++t) {
      //get transformation
      Map<const Matrix2d> T(data_in + 4*t);

      //closed form solution
      double trace        = T.trace();
      double offdiff      = T(0,1) - T(1,0);
      double off_term     = offdiff / trace;
      double denom_global = 1 / sqrt(off_term*off_term + 1);

      Matrix2d R;
      R(0,0) = denom_global; R(0,1) = off_term * denom_global;
      R(1,0) = -R(0,1);      R(1,1) = R(0,0);

      if(R.determinant() < 0) //safety checks
        R *= 1;

      Map<Matrix2d>(data_out + 4*t) = R;
    }

    } else {
      //hope for better compiler optimization
      for(mwSize t = 0; t < nt; ++t) {
        //get transformation
        Map<const Matrix2d> T(data_in + 4*t);

        //get svd
        const JacobiSVD<Matrix2d> svd = T.jacobiSvd(ComputeFullU | ComputeFullV);

        Matrix2d R(svd.matrixU() * svd.matrixV().transpose());

        if(svd.singularValues().prod() < 0) //safety checks
          R *= -1;

        Map<Matrix2d>(data_out + 4*t) = R;
      }
    }

  } else if(m == 3) {
    //hope for better compiler optimization
    for(mwSize t = 0; t < nt; ++t) {
      //get transformation
      Map<const Matrix3d> T(data_in + 9*t);

      //get svd
      const JacobiSVD<Matrix3d> svd = T.jacobiSvd(ComputeFullU | ComputeFullV);

      Matrix3d R(svd.matrixU() * svd.matrixV().transpose());

      if(svd.singularValues().prod() < 0) //safety checks
        R *= -1;

      Map<Matrix3d>(data_out + 9*t) = R;
    }

  } else {
    //generic nd implementation
    for(mwSize t = 0; t < nt; ++t) {
      //get transformation
      Map<const MatrixXd> T(data_in + m*m*t, m, m);

      //get svd
      const JacobiSVD<MatrixXd> svd = T.jacobiSvd(ComputeFullU | ComputeFullV);

      MatrixXd R(svd.matrixU() * svd.matrixV().transpose());

      if(svd.singularValues().prod() < 0) //safety checks
        R *= -1;

      Map<MatrixXd>(data_out + m*m*t, m, m) = R;
    }
  }

  plhs[0] = R;
}
