/**ASSIGN3DLBC This file contains C language functions implementing lower
 *          bounds for the  axial operations research 3D assignment
 *          problem. See the comments to the Matlab implementation of
 *          assign3DLB for more details on the algorithms.
 *
 *March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

//Prototypes for the function implemented here.
#include "assignAlgs3D.h"

//The header for the assignment algorithms in C.
#include "assignAlgs2D.h"

//For basic matrix operations that are simple in matlab, but more
//complicated in C.
#include "basicMatOps.h"

/*We need this for INFINITY to be defined, but if a compiler does not
 *support C99, then it must be explicitly defined.*/
#include <math.h>

//For memcpy and memset
#include <string.h>

//For uint64_t and other types
#include <stdint.h>

/*If a compiler does not support INFINITY in C99, then it must be
 * explicitly defined.*/
#ifndef INFINITY
static const uint64_t infVal=0x7ff0000000000000;
#define INFINITY (*(double*)(&infVal))
#endif

size_t assign3DLBDual0BufferSize(const size_t *nVals) {
/**ASSIGN3DLBDUAL0BUFFERSIZE Determine the minimum amount of memory (in
 *          bytes) needed for the tempSpace input of the
 *          assign3DLBDual0 function. Note that
 *          nVals[0]<=nVals[1]<=nVals[2].
 *
 *March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/

    const size_t n1=nVals[0];
    const size_t n2=nVals[1];
    
    return sizeof(double)*(n1*n2+n2+n1)+sizeof(ptrdiff_t)*(n1+n2)+assign2DCBufferSize(n2,n1);
}

double assign3DLBDual0(const size_t *nVals,const double *C, void *tempSpace) {
/**ASSIGN3DLBDUAL0 A C implementation of the lower bound of the axial
 *         operations research 3D assignment problem obtained by using the
 *         dual cost of the Lagrangian relaxation down to 2D assignment,
 *         computed setting the dual variables all to 0.
 *
 *INPUTS: nVals A length 3 array giving the side of each dimension of the
 *              3D array C. nVals[0]*nVals[1]*nVals[2]>=1. It is assumed
 *              that nVals[0]<=nVals[1]<=nVals[2];
 *            C A pointer to an nVals[0]XnVals[1]XnVals[2] 3D array stored
 *              by column as is the ordering in Matlab and Fortran.
 *    tempSpace A buffer of at least assign3DLBDual0BufferSize(nVals)
 *              bytes to be used for intermediate computations.
 *
 *OUTPUTS: The return value is a lower bound on the value of the axial
 *         operations research 3D assignment problem.
 *
 *March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    
    const size_t n1=nVals[0];
    const size_t n2=nVals[1];
    const size_t n3=nVals[2];
    const size_t n1n2=n1*n2;
    const size_t numElsInC=n1*n2*n3;
    double *d2, *u2D, *v2D;
    ptrdiff_t *row4col, *col4row;
    void *tempBuffer2DAssign;
    bool isInfeasible;
    size_t i1, i2, i3;
    double q;
    
    //If the problem is scalar.
    if(numElsInC==1) {
        return C[0];
    }

    d2=(double*)tempSpace;//Size n1*n2;
    u2D=d2+n1*n2;//Length n2
    v2D=u2D+n2;//Length n1
    row4col=(ptrdiff_t*)(v2D+n1);//Length n2
    col4row=row4col+n2;//Length n1.
    //assign2DC requires this buffer to be at least
    //assign2DCBufferSize(n2,n1) bytes in size.
    tempBuffer2DAssign=(void*)(col4row+n1);
    
//////
//DUAL COST WITH 0 DUAL VARIABLES.
//////
    //These loops essentially do:
    //d2=min(C,[],3);
    //d2=d2';
    //The transpose on d2 is necessary, because assign2DC requires that
    //the number of rows of the input be >= the number of columns and
    //we know that n2>=n1.
    for(i2=0;i2<n2;i2++) {
        for(i1=0;i1<n1;i1++) {
            double minVal=C[i1+n1*i2];

            for(i3=1;i3<n3;i3++) {
                double curVal=C[i1+n1*i2+n1n2*i3];

                if(curVal<minVal) {
                    minVal=curVal;
                }
            }

            //Store in a transposed order.
            d2[i2+n2*i1]=minVal;
        }
    }
    
    //This function modifies (the transposed) d2.
    isInfeasible=assign2DC(false, d2, &q, row4col,col4row, tempBuffer2DAssign, v2D, u2D, n2, n1);
    
    if(isInfeasible) {
        //If the 2D assignment problem is infeasible.
        return (double)INFINITY;
    }

    return q;
} 
    
size_t assign3DLBHungarianBufferSize(const size_t *nVals) {
/**ASSIGN3DLBHUNGARIANBUFFERSIZE Determine the minimum amount of memory (in
 *          bytes) needed for the tempSpace input of the
 *          assign3DLBHungarian function. Note that
 *          nVals[0]<=nVals[1]<=nVals[2].
 *
 *March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
    
    const size_t n1=nVals[0];
    const size_t n2=nVals[1];
    const size_t n3=nVals[2];
    const size_t numElsInC=n1*n2*n3;
    size_t bufferSize;
    
    bufferSize=sizeof(double)*(numElsInC+3*n3*n3+2*n3)+sizeof(ptrdiff_t)*(2*n3);
    return bufferSize+assign2DCBufferSize(n3,n3);
}

double assign3DLBHungarian(const size_t *nVals,const double *COrig, void *tempSpace) {
/**ASSIGN3DLBHUNGARIAN A C implementation of the Algorithm of [1] to obtain
 *      a lower bound on the axial operations research 3D assignment
 *      problem.
 *
 *INPUTS: nVals A length 3 array giving the side of each dimension of the
 *              3D array C. nVals[0]*nVals[1]*nVals[2]>=1. It is assumed
 *              that nVals[0]<=nVals[1]<=nVals[2];
 *            C A pointer to an nVals[0]XnVals[1]XnVals[2] 3D array stored
 *              by column as is the ordering in Matlab and Fortran.
 *    tempSpace A buffer of at least assign3DLBHungarianBufferSize(nVals)
 *              bytes to be used for intermediate computations.
 *
 *OUTPUTS: The return value is a lower bound on the value of the axial
 *         operations research 3D assignment problem.
 *   
 *REFERENCES:
 *[1] B.-J. Kim, W. L. Hightower, P. M. Hahn, Y.-R. Zhu, and L. Sun, "Lower
 *    bounds for the axial three-index assignment problem," European
 *    Journal of Operational Research, vol. 202, no. 3, pp. 654-668, May
 *    2010.
 *
 *March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */

    const size_t n1=nVals[0];
    const size_t n2=nVals[1];
    const size_t n3=nVals[2];
    const size_t numElsInC=n1*n2*n3;
    
    double lowerBound;
    double *C, *M, *MTemp, *MTrans;
    double *u,*v;//To hold dual variables.
    ptrdiff_t *col4row, *row4col;
    void *tempBuffer2D;
    size_t curIdx, i;

    //Allocate the buffer.
    C=(double*)tempSpace;//Size numElsInC.
    M=C+numElsInC;//Maximum size n3*n3.
    MTemp=M+n3*n3;//Maximum size n3*n3.
    MTrans=MTemp+n3*n3;//Maximum size n3*n3.
    u=MTrans+n3*n3;//Maximum length n3.
    v=u+n3;//Maximum length n3.
    col4row=(ptrdiff_t*)(v+n3);//Maximum length n3.
    row4col=col4row+n3;//Maximum length n3.
    //Maximum buffer size is assign2DCBufferSize(n3,n3).
    tempBuffer2D=(void*)(row4col+n3);

    //Copy COrig into C.
    memcpy(C,COrig,numElsInC*sizeof(double));

    lowerBound=0;
    for(curIdx=0;curIdx<3;curIdx++) {
        size_t numRow,numCol;//The dimensions of the M matrix.
        
        //curIdx is the index in C over which M holds the minimum values
        //over the other two indices.
        
        //The algorithm [1] is written for matrices with the same number
        //of elements in each dimensions, so we want to add dummy 
        //entries to the matrix if the matrix is rectangular to make it
        //n3Xn3Xn3. Thus, if n1<=n2<=n3, this means that that all of the
        //extra i values assigned have zero cost for all entries of the
        //matrix, and all of the extra j values added have infinite
        //cost for all of the existing i1 values so that they are never
        //assigned. However, it is not necessary to add all of those
        //values to the full matrix. By adding extra i values to
        //C(i,j,k), when curIdx=1, the M matrix is zero, C does not
        //change, and lowerBound=0, so that iteration is not necessary.
        //When curIdx=2, then the M matrix is used, but the infinite
        //values added to C for the extra j entries never show up in M,
        //(as it takes a minimum) and the extra i entries are all
        //zero. The extra i entries matter, because they lengthen the
        //dual vector, which changes the dual costs. With curIndex=3,
        //the dual costs do not matter and the extra entries in M could
        //have been omitted.
            
        //Thus, there is no need to actually add the extra elements to
        //C when n1~=n2 or n2~=n3 with n1<=n2<=n3. Rather, one can just
        //skip the first iteration, augment M with n3-n1 extra rows of
        //zeros for the second iteration, and then use a rectangular M
        //for the final iteration, which does not change any dual values.
        
        if((n1!=n2||n2!=n3)&&(curIdx==0)) {
            continue;
        }
        
       //Step 1: Create the M matrix and modify the C matrix.
        if(curIdx==0) {
            //For M:
            numRow=n2;
            numCol=n3;
            
            // M=reshape(min(C,[],curIdx),[numRow,numCol]);
            minMatOverDimCDouble(3,nVals,M,C,curIdx);
            
            for(i=0;i<n1;i++) {
                //C(i,:,:)=C(i,:,:)-reshape(M,[1,numRow,numCol]);
                fixedIdxMatSub(3, nVals,M, C,0, i);
            }
        } else if(curIdx==1) {
            numRow=n1;
            numCol=n3;
            
            //M is numColXnumCol. We zero-pad all rows after numRow.
            //The augmented M is used in the 2D assignment step.
            if(numCol!=numRow) {
                for(i=0;i<numCol;i++) {
                    memset(M+numRow+i*numCol,0,(numCol-numRow)*sizeof(double));
                }
            }

            //MTemp=reshape(min(C,[],curIdx),[numRow,numCol]);
            minMatOverDimCDouble(3,nVals,MTemp,C,curIdx);

            //M(1:numRow,:)=MTemp
            for(i=0;i<numCol;i++) {
                memcpy(M+i*numCol,MTemp+i*numRow,numRow*sizeof(double));
            }

            for(i=0;i<n2;i++) {
                //C(:,i,:)=C(:,i,:)-reshape(MTemp,[numRow,1,numCol]);
                //This is the same as
                //C(:,i,:)=C(:,i,:)-reshape(M(1:numRow,:),[numRow,1,numCol]);
                fixedIdxMatSub(3,nVals,MTemp, C,1, i);
            }

            //The non-augmented M, MTemp, is used again after the 2D
            //assignment step for modifying C.
        } else {//curIdx==2
            numRow=n1;
            numCol=n2;
            
            minMatOverDimCDouble(3,nVals,M,C,curIdx);
        }
        
        //Step 2: Find the minimum assignment and compute the reduced M
        //matrix.
        //The reduced cost matrix M will not be identical to that obtained
        //using the traditional Hungarian method ([ope1] suggested using the
        //Hungarian method) but it should be suitable to serve the same
        //purpose for reducing the C matrix.
        
        //Copy M into MTemp, because assign2DC will modify it. Also, the
        //copy is does transposed, because numCol>=numRow and assign2DC
        //requires the opposite relation.
        
        {
            bool isInfeasible;
            double gain;
            
            
            if(curIdx==1||numRow==numCol) {
                //M is a numColXnumCol augmented matrix if curIdx==1.
                isInfeasible=assign2DC(false, M, &gain, col4row, row4col, tempBuffer2D, u, v, numCol, numCol);
            } else { 
                size_t j;
                for(j=0;j<numCol;j++) {
                    for(i=0;i<numRow;i++) {
                        MTrans[j+numCol*i]=M[i+numRow*j];
                    }
                }
                //Since the matrix is transposed, we swap col4row, row4col
                //inputs as well as u and v inputs.
                isInfeasible=assign2DC(false, MTrans,&gain, row4col, col4row, tempBuffer2D, v, u, numCol, numRow);
            }

            //If the assignment problem is infeasible, return Inf cost.
            if(isInfeasible) {
                return (double)INFINITY;
            }
            lowerBound+=gain;
        }
        
        if(curIdx!=2) {
            double minM;
            size_t curRow, curCol;
            double *MCur;
            //The dual variables returned by assign2D are not the same as
            //as those returned by a true shortest augmenting path
            //algorithm (the Hungarian algorithm), because a preprocessing
            //step offsets them to make the algorithm work with C matrices
            //with negative entries. Adding the minimum element of M to the
            //u dual variables gets rid of this offset.
            
            //minM=min(0,min(M(:)));
            if(curIdx==1) {
                MCur=MTemp;//Don't use the augmented matrix.
            } else {
                MCur=M;
            }
            
            minM=vecMinD(MCur,numRow*numCol);
            if(minM>0) {
                minM=0;
            }
            
            //u=u+minM;
            for(curCol=0;curCol<numCol;curCol++) {
                u[curCol]+=minM;
            }

            //Create the reduced cost matrix of M using u and v.
            i=0;
            for(curCol=0;curCol<numCol;curCol++){
                for(curRow=0;curRow<numRow;curRow++) {
                    MCur[i]=-MCur[i]+u[curCol]+v[curRow];
                    i++;
                }
            }

            //Step 3: Modify the C matrix
            if(curIdx==0) {
                for(i=0;i<n1;i++) {
                    //C(i,:,:)=C(i,:,:)-reshape(MCur,[1,numRow,numCol]);
                    fixedIdxMatSub(3,nVals,MCur,C,0,i);
                }
            } else {//curIdx==1; it does not get here if curIdx==2.
                for(i=0;i<n2;i++) {
                    //C(:,i,:)=C(:,i,:)+reshape(MCur,[numRow,1,numCol]);
                    fixedIdxMatSub(3,nVals,MCur,C,1,i);
                }
            }
        }
    }
    
    return lowerBound;
}

double assign3DLBCPierskalla(const size_t *nVals,const double *C) {
/**ASSIGN3DLBPIERSKALLA This is a C implementation of the Pierskalla lower
 *          bound for the axial operations research 3D assignment given in
 *          [1].
 *
 *INPUTS: nVals A length 3 array giving the side of each dimension of the
 *              3D array C. nVals[0]*nVals[1]*nVals[2]>=1. It is assumed
 *              that nVals[0]<=nVals[1]<=nVals[2];
 *            C A pointer to an nVals[0]XnVals[1]XnVals[2] 3D array stored
 *              by column as is the ordering in Matlab and Fortran.
 *
 *OUTPUTS: The return value is a lower bound on the value of the axial
 *         operations research 3D assignment problem.
 *
 *REFERENCES:  
 *[1] W. P. Pierskalla, "The tri-substitution method for the three-
 *   dimensional assignment problem," Canadian Operational Research Society
 *   Journal, vol. 5, no. 2, pp. 71-81, Jul. 1967.
 *
 *March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    
    const size_t n1=nVals[0];
    const size_t n2=nVals[1];
    const size_t n3=nVals[2];
    const size_t n2n3=n2*n3;
    size_t i;
    double lowerBound;
    
    lowerBound=0;
    for(i=0;i<n1;i++) {
        double minVal;
        size_t j;
        
        minVal=C[i];
        for(j=1;j<n2n3;j++) {
            const double curVal=C[i+j*n1];
            if(curVal<minVal) {
                minVal=curVal;
            }
        }

        lowerBound+=minVal;
    }
    
    return lowerBound;
}

/*LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.*/
