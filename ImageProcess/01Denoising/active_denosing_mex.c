#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "mex.h"
#include "matrix.h"


#define X(ix, iy) (ix)* iNy + (iy)
#define SQR(x) (x) * (x)
#define ABS(x) ( (x) > 0.0 ? x : (-x) )

float MatSum(float* A, float* B, int startB, int endB){
    float sum = 0.0;
    for (int j = startB; j <= endB-1; j++) {
        sum +=  A[j] * B[j];   
    }
    return sum;
}

float POW(float x, int y){
    float sum = 1;
    for(int j = 0; j < y; j++){
        sum *= x;
    }
    return sum;
}


float shrinkage(float x, float iMu) {
    // shrinkage operator
    if (x > iMu){
        return x - iMu;
    }
    else if (x < -iMu) {
        return x + iMu;
    }
    else {
        return (float)0;
    }   
}

void ddz(float* ddxyz, float* dz, int divx, int divy, int iNx, int iNy){
    float fDxv = 0.0;
    float fDyv = 0.0;
    int iX = 0;
    int diX = 0;
    for (int ix = 0; ix < iNx - 1; ix++){
        for (int iy = 0; iy < iNy - 1; iy++){
            iX = X(ix, iy);
            diX = X(ix+divx, iy+divy);
            ddxyz[iX] = dz[diX] - dz[iX];       
        }
    }
}


void ComputeR(float* r, float* pfIm0, float* u, float iMu, int iNx, int iNy){
    // Compute r^k
    for (int ix = 0; ix < iNx;ix++){
        for (int iy = 0; iy < iNy;iy++){
            r[X(ix,iy)] = shrinkage(pfIm0[X(ix,iy)] - u[X(ix,iy)], iMu);
        }
    }
}

void ComputeZ(float* dxz, float* dyz, float* v, float iEta, int iNx, int iNy){
    // Compute z^k
    float fDxv = 0.0;
    float fDyv = 0.0;
    int iX = 0;
    for (int ix = 0; ix < iNx - 1; ix++){
        for (int iy = 0; iy < iNy - 1; iy++){
            iX = X(ix, iy);
            fDxv = v[X(ix+1, iy)] - v[iX];
            fDyv = v[X(ix, iy+1)] - v[iX];  
            dxz[iX] = shrinkage(fDxv, iEta);
            dyz[iX] = shrinkage(fDyv, iEta);        
        }
    }
}

void ComputeU(float* u, float* r, float* pfIm0, float* v, float* w, float Lambda, float iMu, int iTheta, int iNx, int iNy){
    // Compute u 
    int iX;
    float temp = 0.0; 
    temp = Lambda + iMu * iTheta;
    for (int ix = 0; ix < iNx; ix++){
        for (int iy = 0; iy < iNy; iy++){
            iX = X(ix, iy);
            u[iX] = iMu * iTheta * (v[iX] - w[iX]) + Lambda * (pfIm0[iX] - r[iX]);
        }
    }
}

void Laplacian(float* Intermediate, float* DiffLapInter, int iNx, int iNy){
    int iX = 0;
    for (int ix = 1; ix < iNx - 1; ix++){
        for (int iy = 1; iy < iNy - 1; iy++){
            // Lap = f(i+1,j) + f(i-1,j) + f(i,j+1) + f(i,j-1) - 4f(i,j)
            iX = X(ix, iy);
            DiffLapInter[iX] = (float) (Intermediate[X(ix+1, iy)] + Intermediate[X(ix-1, iy)] + Intermediate[X(ix, iy+1)] + Intermediate[X(ix, iy-1)] - 4.0 * Intermediate[iX]);
        }
    }
}

void CreateCoef(float* A, float* DiffLapInter, int iNx, int iNy){
    for (int i = 1; i < iNx - 1; i++){
        A[i * iNy + i] = 5;  // diagonal element       
        A[i*iNy+i-1] = DiffLapInter[i*iNy+i-1];
        A[i*iNy+i+1] = DiffLapInter[i*iNy+i+1];
    }
}


void ComputeInterV(float* A, float* tempA, float Xi, float* vtemp, float* DiffLapInterV, float* IntermediateV, float* pfIm0, float* w, float* dxz, float* dyz, float* ddxz, float* ddyz,float Lambda, float iEta, int iTheta, int iNx, int iNy){
    float fDxz = 0.0;
    float fDyz = 0.0;
    int iX = 0;
    float sum_old = 0.0;
    float sum_update = 0.0;
    float error;
    float x_i_old;
    Xi = (float) ((float) (1.0 - Lambda) / (float) (POW(iEta, iTheta)));
    ddz(ddxz, dxz, 1, 0, iNx, iNy);
    ddz(ddyz, dyz, 0, 1, iNx, iNy);
    for (int ix = 0; ix < iNx; ix++){
        for (int iy = 0; iy < iNy; iy++){
            iX = X(ix, iy);            
            vtemp[iX] = pfIm0[iX] + w[iX] - Xi * (ddxz[iX] + ddyz[iX]);
        }
    }
    Laplacian(IntermediateV, DiffLapInterV, iNx, iNy);

    // // Gauss Seidel iterations
    // // Create coefficient matrix
    CreateCoef(A, DiffLapInterV, iNx, iNy);

    for (int k = 0; k < 1000; k++){      //  number of iteration 1000
        error = 0.0;
        x_i_old = 0.0;
        for (int i = 0; i < iNx ; i++){
            x_i_old = IntermediateV[i];
            for (int j = 0; j < iNy; j++){
                tempA[j] = A[i * iNx + j];
            }
            sum_update = MatSum(tempA, IntermediateV, 0, i);
            sum_old = MatSum(tempA, IntermediateV, i+1, iNy);
            IntermediateV[i] = 1 / A[i * iNy+i] * (vtemp[i] - sum_old - sum_update);   //?

            error += ABS(IntermediateV[i] - x_i_old);  // ?
        }
        if (error / (iNy) < 1e-8){
            break;
        }
        Laplacian(IntermediateV, DiffLapInterV, iNx, iNy);
        CreateCoef(A, DiffLapInterV, iNx, iNy);
    }    
}

void ComputeW(float* w, float* pfIm0, float* pfv, int iNx, int iNy){
    for (int ix = 0; ix < iNx; ix++){
        for (int iy = 0; iy < iNy; iy++){
            w[X(ix, iy)] += pfIm0[X(ix, iy)] - pfv[X(ix,iy)];
        }
    }
}


float ComputeUpsilon(float* pfIm0, float* pfu, float* r, float iMu, int iBeta,int iNx, int iNy){
    // compute Upsilon
    float temp = 0.0;
    float upsilon = 0.0;
    for (int ix = 0; ix < iNx; ix++){
        for (int iy = 0; iy < iNy; iy++){
            temp = (float) (ABS(r[X(ix, iy)]) + 1.0 / (2 * iMu) * SQR((pfIm0[X(ix, iy)] - pfu[X(ix, iy)] - r[X(ix, iy)])));  // ?
            temp = (float) (exp(- (temp * pfu[X(ix,iy)]) / iBeta));
            if (ix == 0 && iy == 0){
                upsilon = temp;
            }
            if (temp < upsilon){
                upsilon = temp;
            }
            // if (temp > upsilon){
            //     upsilon = temp;
            // }
        }
    } 
    return upsilon;
}

/***********************************/
/*********MAIN FUNCTION*************/
/***********************************/
void mexFunction(int iNbOut, mxArray *pmxOut[], int iNbIn, const mxArray *pmxIn[]){
    // iNbOut: number of  outputs
    // pmxOut: array of pointness to output arguments
    // iNbIn: number of  inputs
    // pmxIn: array of pointness to input arguments

    float *pfIm0, *pfVecParameters, *pfuj, Lambda, *r, *pfdxz, *pfdyz, *w;
    float *pfu, fSumImRef, *pfv, upsilon;
    float pfXi, *vtemp, *DiffLapInterV, *d, *z, *pfddxz, *pfddyz;
    float iMu, iEta, iAlpha;
    float *A, *tempA;   // In order to arrange neatly, ** is not used
    int iNy, iNx, iIter, iBeta, iTheta;

    /* Inputs */
    pfIm0 = mxGetData(pmxIn[0]);  // Given image
    iNx = (int) mxGetM(pmxIn[0]);
    iNy =(int) mxGetN(pmxIn[0]);
    pfv = mxGetData(pmxIn[1]);   // u
    pfVecParameters = mxGetData(pmxIn[2]);  // Vector of parameters


    iIter = (int) pfVecParameters[0];
    iMu = (float) pfVecParameters[1]; 
    iEta = (float) pfVecParameters[2]; 
    iAlpha = (float) pfVecParameters[3];
    iBeta = (int) pfVecParameters[4];
    iTheta = (int) pfVecParameters[5];


    /* Outputs */
    pmxOut[0] = mxCreateNumericArray(mxGetNumberOfDimensions(pmxIn[0]), mxGetDimensions(pmxIn[0]), mxSINGLE_CLASS, mxREAL);  // ? w 4133
    pfu =(float*) mxGetPr(pmxOut[0]);    

    /* Memory allocation */
    pfddxz = (float*) calloc(iNy * iNx, sizeof(float));
    if(!pfddxz){
        mexErrMsgTxt("Memory allocation failure\n");
    }
    pfddyz = (float*) calloc(iNy * iNx, sizeof(float));
    if(!pfddyz){
        mexErrMsgTxt("Memory allocation failure\n");
    }
    A = (float*) calloc(iNy * iNx, sizeof(float));
    if(!A){
        mexErrMsgTxt("Memory allocation failure\n");
    }
    pfv = (float*) calloc(iNy * iNx, sizeof(float));
    if(!pfv){
        mexErrMsgTxt("Memory allocation failure\n");
    }
    pfuj = (float*) calloc(iNy * iNx, sizeof(float));
    if(!pfuj){
        mexErrMsgTxt("Memory allocation failure\n");
    }
    r = (float*) calloc(iNy * iNx, sizeof(float));
    if(!r){
        mexErrMsgTxt("Memory allocation failure\n");
    }
    pfdxz = (float*) calloc(iNy * iNx, sizeof(float));
    if(!pfdxz){
        mexErrMsgTxt("Memory allocation failure\n");
    }
    pfdyz = (float*) calloc(iNy * iNx, sizeof(float));
    if(!pfdyz){
        mexErrMsgTxt("Memory allocation failure\n");
    }
    w = (float*) calloc(iNy * iNx, sizeof(float));
    if(!w){
        mexErrMsgTxt("Memory allocation failure\n");
    }
    vtemp = (float*) calloc(iNy * iNx, sizeof(float));
    if(!vtemp){
        mexErrMsgTxt("Memory allocation failure\n");
    }
    DiffLapInterV = (float*) calloc(iNy * iNx, sizeof(float));
    if(!DiffLapInterV){
        mexErrMsgTxt("Memory allocation failure\n");
    }
    d = (float*) calloc(iNy * iNx, sizeof(float));
    if(!d){
        mexErrMsgTxt("Memory allocation failure\n");
    }
    z = (float*) calloc(iNy * iNx, sizeof(float));
    if(!z){
        mexErrMsgTxt("Memory allocation failure\n");
    }
    tempA = (float*) calloc(iNy * iNx, sizeof(float));
    if(!tempA){
        mexErrMsgTxt("Memory allocation failure\n");
    }

    for (int ix = 0; ix < iNx; ix++){
        for (int iy = 0; iy < iNy; iy++){
            pfu[X(ix,iy)] = pfv[X(ix,iy)];
        }
    } 
    
    /* Parameters for the segmentation code */  
    for (int iter = 0;iter < iIter; iter++){
        upsilon = ComputeUpsilon(pfIm0, pfu, r, iMu, iBeta, iNx, iNy);
        // mexPrintf("upsilon = %f\n", upsilon);
        Lambda = shrinkage(upsilon, iAlpha);  // ? 2400
        // mexPrintf("Lambda = %f\n", Lambda);
        ComputeR(r, pfIm0, pfu, iMu, iNx, iNy);
        ComputeZ(pfdxz, pfdyz, pfv, iEta, iNx, iNy);
        ComputeU(pfu, r, pfIm0, pfv, w, Lambda, iMu, iTheta, iNx, iNy);
        ComputeInterV(A, tempA, pfXi, vtemp, DiffLapInterV, pfv, pfIm0, w, pfdxz, pfdyz, pfddxz, pfddyz, Lambda, iEta, iTheta, iNx, iNy); 
        ComputeW(w, pfIm0, pfv, iNx, iNy);                      
    }
    mexPrintf("Done!\n");
    for (int ix = 0; ix < iNx; ix++){
        for (int iy = 0; iy < iNy; iy++){
            pfu[X(ix,iy)] = ABS(pfu[X(ix,iy)]);
        }
    }  
    free(pfddxz);
    free(pfddyz);
    free(A);
    free(pfv);
    free(r);
    free(pfdxz);
    free(pfdyz);
    free(w);
    free(vtemp);
    free(DiffLapInterV);
    free(d);
    free(z);
    free(tempA);
}