#include "mex.h"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // validate arguments
    if (nrhs!=2)
        mexErrMsgIdAndTxt("mex:error", "Wrong number of arguments");
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]))
        mexErrMsgIdAndTxt("mex:error", "Input isnt real dense double array");
    if (mxGetNumberOfElements(prhs[0]) != mxGetNumberOfElements(prhs[1]))// || mxGetNumberOfElements(prhs[0]) != mxGetNumberOfElements(prhs[2]))
        mexErrMsgIdAndTxt("mex:error", "input must have same number of elements");
    
    
    const mwSize *ntpts=mxGetDimensions(prhs[0]);
    // allocate ATAput
    plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL); // [Fp Vp PS]
    plhs[1] = mxCreateDoubleMatrix(3, 1, mxREAL); // vectorB
    plhs[2] = mxCreateDoubleMatrix(3, 3, mxREAL); // ATA
    plhs[3] = mxCreateDoubleMatrix(ntpts[0], 3, mxREAL); // vectorA
    
    
    double *Param = mxGetPr(plhs[0]);
    double *Bvec = mxGetPr(plhs[1]);
    double *ATAinv = mxGetPr(plhs[2]);
    double *Amat = mxGetPr(plhs[3]);
    
    const double *C_vec = mxGetPr(prhs[0]);
    const double *Cp_vec = mxGetPr(prhs[1]);
    
    //populate Amat
    Amat[0]=0;
    Amat[ntpts[0]]=0;
    for (int i=1; i<ntpts[0]; i++){
        Amat[i+0]=Amat[(i-1)+0]-((C_vec[i]+C_vec[i-1])/2);
        Amat[i+ntpts[0]]=Amat[(i-1)+ntpts[0]]+(Cp_vec[i]+Cp_vec[i-1])/2;
    }
    
    Amat[2*ntpts[0]]=Amat[ntpts[0]];
    for (int i=1; i<ntpts[0]; i++){
        Amat[i+2*ntpts[0]]=Amat[(i-1)+2*ntpts[0]]+(Amat[i+1*ntpts[0]]+Amat[i-1+1*ntpts[0]])/2;
    }
    
    for (int row=0; row<3; row++){
        for (int col=0; col<3; col++){
            for (int i=0; i<ntpts[0]; i++){
                ATAinv[3*row+col]=ATAinv[3*row+col]+Amat[i+row*ntpts[0]]*Amat[i+col*ntpts[0]];
            }
        }
    }
    // 3x3 input matrix [a b c; d e f; g h i], and its determinant
    const double A = (ATAinv[4]*ATAinv[8]-ATAinv[5]*ATAinv[7]);
    const double B = -(ATAinv[3]*ATAinv[8]-ATAinv[5]*ATAinv[6]);
    const double C = (ATAinv[3]*ATAinv[7]-ATAinv[4]*ATAinv[6]);
    const double D = -(ATAinv[1]*ATAinv[8]-ATAinv[2]*ATAinv[7]);
    const double E = (ATAinv[0]*ATAinv[8]-ATAinv[2]*ATAinv[6]);
    const double F = -(ATAinv[0]*ATAinv[7]-ATAinv[1]*ATAinv[6]);
    const double G = (ATAinv[1]*ATAinv[5]-ATAinv[2]*ATAinv[4]);
    const double H = -(ATAinv[0]*ATAinv[5]-ATAinv[2]*ATAinv[3]);
    const double I = (ATAinv[0]*ATAinv[4]-ATAinv[1]*ATAinv[3]);    
    
    const double det = (ATAinv[0]*A + ATAinv[1]*B + ATAinv[2]*C);
    
    if (det != 0) {       
        ATAinv[0] =  A/det;
        ATAinv[1] =  D/det;
        ATAinv[2] =  G/det;
        ATAinv[3] =  B/det;
        ATAinv[4] =  E/det;
        ATAinv[5] =  H/det;
        ATAinv[6] =  C/det;
        ATAinv[7] =  F/det;
        ATAinv[8] =  I/det;
    }    
      
    for (int col=0; col<3; col++){
            for (int i=0; i<ntpts[0]; i++){
                Param[col]=Param[col]+Amat[i+ntpts[0]*col]*C_vec[i];
            }
    }

    for (int col=0; col<3; col++){            
        for (int row=0; row<3; row++){            
                Bvec[col]=Bvec[col]+ATAinv[row+col*3]*Param[row];
                
        }
    }
    Param[0]=Bvec[1];
    Param[1]=(Bvec[1]*Bvec[1])/(Bvec[0]*Bvec[1]-Bvec[2]);
    Param[2]=Param[1]*Bvec[2]/Bvec[1];
}