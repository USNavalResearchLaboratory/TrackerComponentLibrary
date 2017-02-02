/*Code implementing the MMOSPA optimization algorithm in C++.
 *
*November 2013 David F. Crouse, Naval Research Laboratory, Washington D.C*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include <limits>
#include <algorithm>
#include "MMOSPAApproxCPP.hpp"

using namespace std;

//Prototypes for functions not prototyped in the headers.
void MMOSPAApproxForward(double *MMOSPAEst,
                         size_t *orderList,
                         MurtyHyp *problemSol,
                         ScratchSpace &workMem,
                         const double *x,
                         const double *w,
                         const size_t xDim,
                         const size_t numTar,
                         const size_t numHyp);

void doUpdate4Col(double *MMOSPAEst,
                  size_t *orderList,
                  MurtyHyp *problemSol,
                  ScratchSpace &workMem,
                  double *xOptCur,
                  const double *x,
                  const double *w,
                  const size_t xDim,
                  const size_t numTar,
                  const size_t varHyp);
                         
inline void multScalVec(double *result, const double scalarVal, const double* vectorVal, const size_t length){
    /* Multiply all of the elements in vectorVal by ScalarVal and save the
     * results in result*/
    size_t i;
    
    for(i=0;i<length;i++){
        result[i]=scalarVal*vectorVal[i];
    }
}

inline double dotProduct(const double *vec1, const double *vec2, const size_t length) {
    size_t i;
    double prodSum=0;
    
    for(i=0;i<length;i++){
        prodSum=prodSum+vec1[i]*vec2[i];   
    }
    
    return prodSum;
}

inline void diffScalVec(const double scalarVal, double *vec, const size_t length) {
    /*This subtracts the scalar value from all elements in the vector and saves the result back to the vector*/
    size_t i;
    
    for(i=0;i<length;i++) {
        vec[i]=scalarVal-vec[i];
    }
}

inline void vecScalMadd(double *vec1, const double scalarVal, const double *vec2, const size_t length) {
    /*vec1+=scalarVal*vec2*/
    size_t i;
    
    for(i=0;i<length;i++){
        vec1[i]+=scalarVal*vec2[i];
    }
}

void MMOSPAApproxForward(double *MMOSPAEst,
                         size_t *orderList,
                         MurtyHyp *problemSol,
                         ScratchSpace &workMem,
                         const double *x,
                         const double *w,
                         const size_t xDim,
                         const size_t numTar,
                         const size_t numHyp) {
/**MMOSPAAPPROXFORWARD  Perform a single forward step of the approximate
*                      MMOSPA algorithm without having any previous
*                      estimate.
*
*INPUTS:    x        An xDim X numTar X numHyp hypermatrix that holds
*                    numHyp hypotheses each consisting or numTar targets
*                    (or generic vectors) with xDim dimensions per target
*                    (per generic vector).
*           w        A numHyp X 1 vector of the probabilities of each of
*                    the numHyp hypotheses in x. The elements must all be
*                    positive and sum to one.
*
*OUTPUTS: MMOSPAEst The approximate xDim X numTar MMOSPA estimate.
*         orderList A numTarXnumHyp matrix specifying the ordering of the
*                   targets in each hypothesis that went into the
*                   approximate MMOSPA estimate. 
*
*October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.*/
    size_t curTar1,curTar2, curHyp;
    const size_t stackedTarDim=xDim*numTar;
    const size_t CSize=numTar*numTar;
    
   /* The initial hypothesis has a fixed ordering that determines the 
    * ultimate target ordering.*/
    multScalVec(MMOSPAEst,w[0],x,stackedTarDim);
    for(curTar1=0;curTar1<numTar;curTar1++){
        orderList[curTar1]=curTar1;
    }
    
    /*Enter a loop for all of the other hypotheses*/
    for(curHyp=1;curHyp<numHyp;curHyp++){
        size_t curHypOffset=curHyp*stackedTarDim;
        double cMax=-numeric_limits<double>::infinity();

        /* Fill the cost matrix using the partial MMOSPA estimate up to
         * this point. Rows are from xM (curTar1).
         * Columns are from x.*/
        for(curTar2=0;curTar2<numTar;curTar2++){
            size_t curOffset=curTar2*numTar;
            for(curTar1=0;curTar1<numTar;curTar1++){
                workMem.C[curTar1+curOffset]=dotProduct(MMOSPAEst+xDim*curTar1,x+xDim*curTar2+curHypOffset,xDim);//c(curRow,curCol)=sum(xM(:,curRow).*x(:,curCol,curHyp));

                if(workMem.C[curTar1+curOffset]>cMax){
                    cMax=workMem.C[curTar1+curOffset];
                }
            }
        }
        /* Now, evaluate c=cMax-c; so that we can use a 2D assignment
         * algorithm that performs minimization rather than one that
         * performs maximization.*/
        diffScalVec(cMax,workMem.C, CSize);
        
        /* Now, use the shortest path algorithm to find the optimal 2D
         * assignment. The result is placed in col4row in problemSol.*/
        shortestPathCPP(problemSol,workMem,numTar, numTar, 0);

        /*Record this ordering in orderList.*/
        copy(problemSol->col4row,problemSol->col4row+numTar,orderList+numTar*curHyp);

        /*Add the term with the now known ordering to the MMOSPA estimate.*/
        for(curTar1=0;curTar1<numTar;curTar1++){
            double *vec1Ptr=MMOSPAEst+curTar1*xDim;
            size_t idx=*(orderList+numTar*curHyp+curTar1);
            const double *vec2Ptr=x+curHypOffset+idx*xDim;

            //MMOSPAEst=MMOSPAEst+w(curHyp)*x(:,orderList(:,curHyp),curHyp);
            vecScalMadd(vec1Ptr, w[curHyp], vec2Ptr, xDim);
        }
    }
}

void doUpdate4Col(double *MMOSPAEst,
                  size_t *orderList,
                  MurtyHyp *problemSol,
                  ScratchSpace &workMem,
                  double *xOptCur,
                  const double *x,
                  const double *w,
                  const size_t xDim,
                  const size_t numTar,
                  const size_t varHyp) {
//xOptCur is xDim*numTar X 1 
    const size_t stackedTarDim = xDim*numTar;
    const size_t CSize=numTar*numTar;
    size_t cur1,cur2;
    size_t curHypOffset=varHyp*stackedTarDim;
    size_t curOrderOffset=numTar*varHyp;
    double cMax=-numeric_limits<double>::infinity();
    
    //xOptCur=(xMMOSPAEst-x(:,orderList(:,varHyp),varHyp)*w(varHyp))/denom;
    for(cur1=0;cur1<numTar;cur1++){
        size_t tarOffset=xDim*cur1;
        size_t tarOffsetX=xDim*orderList[cur1+curOrderOffset]+curHypOffset;
        
        //The inner loop copies the elements of the target.
        for(cur2=0;cur2<xDim;cur2++){
            xOptCur[tarOffset+cur2]=(MMOSPAEst[tarOffset+cur2]-x[cur2+tarOffsetX]*w[varHyp]);///denom;
        }
    }
    
    //Now, we must fill the cost matrix.
    for(cur2=0;cur2<numTar;cur2++){
        size_t curOffset=cur2*numTar;
        for(cur1=0;cur1<numTar;cur1++){
            workMem.C[cur1+curOffset]=dotProduct(xOptCur+xDim*cur1,x+xDim*cur2+curHypOffset,xDim);//c(curRow,curCol)=sum(xM(:,curRow).*x(:,curCol,curHyp));

            if(workMem.C[cur1+curOffset]>cMax){
                cMax=workMem.C[cur1+curOffset];
            }
        }
    }
        
    /* Now, evaluate c=cMax-c; so that we can use a 2D assignment
     * algorithm that performs minimization rather than one that
     * performs maximization.*/
    diffScalVec(cMax,workMem.C, CSize);

    /* Now, use the shortest path algorithm to find the optimal 2D
     * assignment. The result is placed in col4row in problemSol.*/
    shortestPathCPP(problemSol,workMem,numTar, numTar, 0);
    
    /*Record this ordering in orderList.*/
    copy(problemSol->col4row,problemSol->col4row+numTar,orderList+curOrderOffset);
    
    /*Update the MMOSPA estimate*/
    //xMMOSPAEst=xOptCur+w(varHyp)*x(:,orderList(:,varHyp),varHyp);
    for(cur1=0;cur1<numTar;cur1++){
        size_t tarOffset=xDim*cur1;
        size_t tarOffsetX=xDim*orderList[cur1+curOrderOffset]+curHypOffset;
        
        //The inner loop copies the elements of the target.
        for(cur2=0;cur2<xDim;cur2++){
            MMOSPAEst[tarOffset+cur2]=xOptCur[tarOffset+cur2]+w[varHyp]*x[cur2+tarOffsetX];
        }
    }
}


void MMOSPAApproxCPP(double *MMOSPAEst,
                     size_t *orderList,
                     const double *x,
                     const double *w,
                     const size_t xDim,
                     const size_t numTar,
                     const size_t numHyp,
                     size_t numScan) {
    MurtyHyp problemSol(numTar,numTar);
    ScratchSpace workMem(numTar,numTar);
    size_t curHyp;
    
    /*Run the forward algorithm*/
    MMOSPAApproxForward(MMOSPAEst,orderList,&problemSol,workMem,x,w,xDim,numTar,numHyp);
    
    if(numScan>1) {
        double *xOptCur = new double[xDim*numTar];//Scratch space for the update steps.
        
        //Do a reverse and then a forward scan numScans times.
        while(numScan>1) {
            //Re-evaluate the hypotheses going backwards.
            for(curHyp=numHyp-2;curHyp>0;curHyp--) {
                doUpdate4Col(MMOSPAEst,
                             orderList,
                             &problemSol,
                             workMem,
                             xOptCur,
                             x,
                             w,
                             xDim,
                             numTar,
                             curHyp);
            }

            //Re-evaluate the hypotheses going forwards.
            for(curHyp=2;curHyp<numHyp;curHyp++){
                doUpdate4Col(MMOSPAEst,
                             orderList,
                             &problemSol,
                             workMem,
                             xOptCur,
                             x,
                             w,
                             xDim,
                             numTar,
                             curHyp);
            }

            numScan--;
        }

        delete[] xOptCur;
    }
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
