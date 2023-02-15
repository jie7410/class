#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "mex.h"
#include "matrix.h"


#define X(ix, iy) (ix) * iNy + (iy)
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
        return 0;
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

float ComputeC(float Lambda, float* pfu, float* pfIm0, float* r, int iNx, int iNy){
    // Compute c^k
    float temp1 = 0.0;
    float temp2 = 0.0;
    int iX = 0;
    for (int ix = 0; ix < iNx; ix++){
        for (int iy = 0; iy < iNy; iy++){
            iX = X(ix, iy);
            temp1 += Lambda * (pfIm0[iX] - r[iX]) * pfu[iX];  // ? 2296
            temp2 += Lambda * pfu[iX];   // ?
        }
    }
    return temp1 / temp2;
}

float LimitU(float* pfu, int iNx, int iNy){
    // Limit u between 0 and 1
    float sumU = 0.0;
    for (int ix = 0; ix < iNx; ix++){
        for (int iy = 0; iy < iNy; iy++){
            sumU += pfu[X(ix,iy)];
        }
    }
    for (int ix = 0; ix < iNx; ix++){
        for (int iy = 0; iy < iNy; iy++){
            pfu[X(ix,iy)] /= sumU;
        }
    }
    return sumU;

}

void ComputeR(float* r, float* pfIm0, float c, float iMu, int iNx, int iNy){
    // Compute r^k
    for (int ix = 0; ix < iNx;ix++){
        for (int iy = 0; iy < iNy;iy++){
            r[X(ix,iy)] = shrinkage(pfIm0[X(ix,iy)] - c, iMu);
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

void ComputeU(float* d, float* IntermediateU, float* u, float* uj, float* r, float* pfIm0, float* v, float* w, float Lambda, float c, float iMu, int iTheta, float iTau, int iNx, int iNy){
    // Compute u
    /*  d^k+1 */
    int iX = 0;
    float SumUj = 0.0;
    for (int ix = 0; ix < iNx; ix++){
        for (int iy = 0; iy < iNy; iy++){
            SumUj += uj[X(ix, iy)];
        }
    }
    mexPrintf("SumUj = %f\n", SumUj);
    for (int ix = 0; ix < iNx; ix++){
        for (int iy = 0; iy < iNy; iy++){
            iX = X(ix, iy);
            d[iX] = (float) (ABS(r[iX]) + 1.0 / (2 * iMu) * SQR(pfIm0[iX] - c - r[iX]));
            IntermediateU[iX] = v[iX] - w[iX] - Lambda / iTheta * d[iX] - iTau / iTheta * SumUj;
            u[iX] = max(0, IntermediateU[iX]);
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


void ComputeInterV(float* A, float* tempA, float Xi, float* vtemp, float* DiffLapInterV, float* IntermediateV, float* pfIm0, float* w, float* dxz, float* dyz, float* ddxz, float* ddyz, float Lambda, float iEta, int iTheta, int iNx, int iNy){
    float fDxz = 0.0;
    float fDyz = 0.0;
    int iX = 0;
    float sum_old = 0.0;
    float sum_update = 0.0;
    float error;
    float x_i_old;
    ddz(ddxz, dxz, 1, 0, iNx, iNy);
    ddz(ddyz, dyz, 0, 1, iNx, iNy);
    Xi = (float) ((float) (1.0 - Lambda) / (float) (POW(iEta, iTheta)));
    for (int ix = 0; ix < iNx-1; ix++){
        for (int iy = 0; iy < iNy-1; iy++){
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
void ComputeV(float* pfv,float* IntermediateV, int iNx, int iNy){
    float V = 0.0;
    for (int i = 0; i < iNx * iNy; i++){
        V += IntermediateV[i];
    }
    for (int ix = 0; ix < iNx; ix++){
        for (int iy = 0; iy < iNy; iy++){
            pfv[X(ix, iy)] = (float) (IntermediateV[X(ix, iy)] - 1.0 / (1) * V);  // ?
        }
    }
}

void ComputeW(float* w, float* pfIm0, float* pfv, int iNx, int iNy){
    for (int ix = 0; ix < iNx; ix++){
        for (int iy = 0; iy < iNy; iy++){
            w[X(ix, iy)] += pfIm0[X(ix, iy)] - pfv[X(ix,iy)];
        }
    }
}

float ComputeUpsilon(float* pfIm0, float* pfu, float* r, float c, float iMu, int iBeta, int iNx, int iNy){
    // compute Upsilon
    float temp = 0.0;
    float upsilon = 0.0;
    for (int ix = 0; ix < iNx; ix++){
        for (int iy = 0; iy < iNy; iy++){
            temp = (float) (ABS(r[X(ix, iy)]) + 1.0 / (2 * iMu) * SQR((pfIm0[X(ix, iy)] - c - r[X(ix, iy)])));  // ?
            temp = (float) (exp(- (temp * pfu[X(ix,iy)]) / iBeta));
            if (upsilon < temp){
                upsilon = temp;
            }
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

    float *pfIm0, *pfIm1, *pfVecParameters, *pfuj, Lambda, *r, *pfdxz, *pfdyz, *w;
    float *BW, *pfu, fSumImRef, *pfv, upsilon, *IntermediateV, *pfddxz, *pfddyz;
    float pfXi, *vtemp, *DiffLapInterV, *d, *IntermediateU, *z;
    float iMu, iEta, iAlpha, iTau;
    float *A, *tempA, c, SumV, SumU, SumUj;   // In order to arrange neatly, ** is not used
    int iNy, iNx, iIter, iBeta, iTheta;
    int iNdim, iDim[2];

    /* Inputs */
    pfIm0 = mxGetData(pmxIn[0]);  // Given image
    iNx = (int) mxGetM(pmxIn[0]);
    iNy =(int) mxGetN(pmxIn[0]);
    BW = mxGetData(pmxIn[1]);   // u
    pfVecParameters = mxGetData(pmxIn[2]);  // Vector of parameters


    iIter = (int) pfVecParameters[0];
    iMu = (float) pfVecParameters[1]; 
    iEta = (float) pfVecParameters[2]; 
    iAlpha = (float) pfVecParameters[3];
    iBeta = (int) pfVecParameters[4];
    iTheta = (int) pfVecParameters[5];
    iTau = (float) pfVecParameters[6];

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
    pfIm1 = (float*) calloc(iNy * iNx, sizeof(float));
    if(!pfIm1){
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
    IntermediateV = (float*) calloc(iNy * iNx, sizeof(float));
    if(!IntermediateV){
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
    // pfXi = (float*) calloc(iNy * iNx, sizeof(float));
    // if(!pfXi){
    //     mexErrMsgTxt("Memory allocation failure\n");
    // }
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
    IntermediateU = (float*) calloc(iNy * iNx, sizeof(float));
    if(!IntermediateU){
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


    /* Compute init function u */
    fSumImRef = LimitU(pfIm0,iNx,iNy);
    fSumImRef = LimitU(pfIm0,iNx,iNy);
    mexPrintf("{%f}\n", fSumImRef);       
    for (int ix = 0; ix < iNy; ix++){
        for (int iy = 0; iy < iNy; iy++){
            if(BW[X(ix,iy)] == 1){
                pfu[X(ix, iy)] = pfIm0[X(ix,iy)];
                pfuj[X(ix, iy)] = 0;
            }
            else{
                pfu[X(ix, iy)] = 0;
                pfuj[X(ix, iy)] = pfIm0[X(ix,iy)];
            }
            pfv[X(ix, iy)] = pfu[X(ix, iy)];  //  u = v
            // pfIm0[X(ix,iy)] /= fSumImRef;
        }          
    }

    /* Parameters for the segmentation code */  
    // SumU = LimitU(pfu, iNx,iNy);
    // SumUj = LimitU(pfuj, iNx,iNy);
    // SumV = LimitU(pfv, iNx, iNy);
    c = 0;
    for (int iter = 0;iter < iIter; iter++){
        mexPrintf("{%d}th iterations,sum:{%d}th\n", iter+1, iIter);
        upsilon = ComputeUpsilon(pfIm0, pfu, r, c, iMu, iBeta, iNx, iNy);
        mexPrintf("upsilon={%f}\n", upsilon);
        Lambda = shrinkage(upsilon, iAlpha);  // ? 2400
        c = ComputeC(Lambda, pfu, pfIm0, r, iNx, iNy);
        mexPrintf("c={%f}\n", c);
        ComputeR(r, pfIm0, c, iMu, iNx, iNy);
        ComputeZ(pfdxz, pfdyz, pfv, iEta, iNx, iNy);
        ComputeU(d, IntermediateU, pfu, pfuj, r, pfIm0, pfv, w, Lambda, c, iMu, iTheta, iTau, iNx, iNy);
        ComputeInterV(A, tempA, pfXi, vtemp, DiffLapInterV, IntermediateV, pfIm0, w, pfdxz, pfdyz, pfddxz, pfddyz, Lambda, iEta, iTheta, iNx, iNy); 
        ComputeW(w, pfIm0, pfv, iNx, iNy);  
        ComputeV(pfv, IntermediateV, iNx, iNy);
    }
    

    for (int ix = 0; ix < iNx; ix++){
        for (int iy = 0; iy < iNy; iy++){
            // pfu[X(ix,iy)] = IntermediateU[X(ix,iy)];
        }
    }
        
    free(pfddxz);
    free(pfddyz);
    free(pfIm1);
    free(A);
    free(pfv);
    free(pfuj);
    free(IntermediateV);
    free(r);
    free(pfdxz);
    free(pfdyz);
    free(w);
    // free(pfXi);
    free(vtemp);
    free(DiffLapInterV);
    free(d);
    free(IntermediateU);
    free(z);
    free(tempA);
}