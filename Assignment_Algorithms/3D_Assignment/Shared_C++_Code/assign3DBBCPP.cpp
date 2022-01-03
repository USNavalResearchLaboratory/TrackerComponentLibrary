/**ASSIGN3DBBCPP A C++ function to solve the axial operations research 3D
 *          assignment problem via a branch-and-bound algorithm. A routine
 *          for determinign the necessary memory as well as subroutines for
 *          some matrix operations are also included.
 * 
 *Better understanding of the algorithms can usually be obtained from
 *looking at the Matlab implementations.
 *
 *November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "assignAlgs3DCPP.hpp"

//For generating and counting tuples.
#include "combinatorialFuns.hpp"

//For std::numeric_limits<T>::infinity() and std::numeric_limits<T>::epsilon()
#include <limits>

//For max, min, and copy.
#include <algorithm>

//For isfinite
#include <cmath>

//For assign3DC, if used to get an initial estimate. 
#include "assignAlgs3D.h"

//For permuteMatrixColsCPP
#include "basicMatOpsCPP.hpp"

//For heapSortVecIdxCDouble
#include "mathFuncsC.h"

//For sub2Idx3D, sub2Idx2D
#include "basicMatOps.h"

size_t assign3DBBBufferSize(const size_t *nDims,const int boundType,const int initMethod) {
/**ASSIGN3DBBBUFFERSIZE Obtain the size (in bytes) of the temporary buffer
 *                      needed by the assign3DBBCPP function.
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    const size_t n1=nDims[0];
    const size_t n2=nDims[1];
    const size_t n3=nDims[2];
    const size_t numLevels=n1;
    const size_t maxTuplesPerLevel=n2*n3;
    size_t buffSize=0;
    
    size_t numCostMatEls=0;
    for(size_t curLevel=1;curLevel<numLevels+1;curLevel++) {
        const size_t offset=curLevel-1;
        numCostMatEls+=(n1-offset)*(n2-offset)*(n3-offset);
    }
    
    buffSize+=2*maxTuplesPerLevel*sizeof(size_t);
    buffSize+=2*maxTuplesPerLevel*numLevels*sizeof(size_t);
    buffSize+=maxTuplesPerLevel*numLevels*sizeof(size_t);
    buffSize+=numLevels*sizeof(size_t);
    buffSize+=2*maxTuplesPerLevel*numLevels*sizeof(size_t);
    buffSize+=maxTuplesPerLevel*numLevels*sizeof(size_t);
    buffSize+=numLevels*sizeof(size_t);
    buffSize+=2*(numLevels-1)*sizeof(size_t);
    buffSize+=numLevels*sizeof(size_t);
    buffSize+=(numLevels-1)*sizeof(double);
    buffSize+=maxTuplesPerLevel*(numLevels-1)*sizeof(double);
    buffSize+=(numLevels-1)*sizeof(double);
    if(boundType==0) {
        buffSize+=assign3DLBHungarianBufferSize(nDims);
    } else if(boundType!=1) {//BoundType==2
        buffSize+=assign3DLBDual0BufferSize(nDims);
    }
    buffSize+=(numLevels+1)*sizeof(size_t);
    buffSize+=numCostMatEls*sizeof(double);
    buffSize+=maxTuplesPerLevel*sizeof(size_t);

    if(initMethod==1) {
        const int subgradMethod=0;
        buffSize=n1*n2*n3*sizeof(double)+std::max(buffSize,n2*sizeof(double)+assign3DCBufferSize(nDims,subgradMethod)+n1*n2*n3*sizeof(double));
    } else {
        buffSize+=n1*n2*n3*sizeof(double);
    }
    
    return buffSize;
}

double assign3DBBCPP(ptrdiff_t *minCostTuples, const double *const COrig,const size_t *const nDims,const int maximize,const int boundType,const int initMethod,const size_t maxIter,const double epsVal, void *tempBuffer) {
/**ASSIGN3DBBCPP Solve the axial operations research 3D assign problem
 *               using a branch-and-bound algorithm.
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    
    double CDelta;
    const size_t n1=nDims[0];
    const size_t n2=nDims[1];
    const size_t n3=nDims[2];
    const size_t numEls=n1*n2*n3;
    uint8_t *curBuffPtr=reinterpret_cast<uint8_t*>(tempBuffer);
    
    if(n1==1&&n2==1&&n3==1) {
        //The special scalar case.
        minCostTuples[0]=0;
        minCostTuples[1]=0;
        minCostTuples[2]=0;
        
        return COrig[0];
    }
    
    //For some types of lower bounds, the cost matrix must have all non-
    //negative elements for the assignment algorithm to work. This forces
    //all of the elements to be positive. The delta is added back in when
    //computing the gain in the end.
    double *C=reinterpret_cast<double*>(curBuffPtr);
    curBuffPtr+=n1*n2*n3*sizeof(double);
           
    if(maximize==true) {
        //CDelta=max(C(:));
        CDelta=COrig[0];
        for(size_t i=1;i<numEls;i++) {
            CDelta=std::max(COrig[i],CDelta);
        }
        
        //C=-C+CDelta;
        for(size_t i=0;i<numEls;i++) {
            C[i]=-COrig[i]+CDelta;
        }
    } else{
        //CDelta=min(C(:));
        CDelta=COrig[0];
        for(size_t i=1;i<numEls;i++) {
            CDelta=std::min(COrig[i],CDelta);
        }
        //C=C-CDelta;
        for(size_t i=0;i<numEls;i++) {
            C[i]=COrig[i]-CDelta;
        }
    }
    
    double q;
    double gain;
    
    //Get an initial estimate.
    switch(initMethod) {
        case 0://Do not use an initial estimate.
        {
            q=-std::numeric_limits<double>::infinity();
            gain=std::numeric_limits<double>::infinity();
            break;
        }
        default://initMethod=1
        {
            //This section implements:
            //[minCostTuples,gain,q]=assign3D(C,false,[],[],maxIter,epsVal,epsVal);
            const int subgradMethod=0;
            double *u=reinterpret_cast<double*>(curBuffPtr);
            uint8_t *initBuffPtr=curBuffPtr+n2*sizeof(double);
            
            const size_t bufferSize=assign3DCBufferSize(nDims,subgradMethod);
            
            void *tempSpace=reinterpret_cast<void*>(initBuffPtr);
            initBuffPtr+=bufferSize;

            //assign3DC modified C, so we have to copy it.
            double *CCopy=reinterpret_cast<double*>(initBuffPtr);
            initBuffPtr+=n1*n2*n3*sizeof(double);
            
            std::copy(C,C+n1*n2*n3,CCopy);
            
            const bool maximize=false;
            const double param1=1;//gammaParam
            const double param2=0;
            const size_t param3=0;
 
            const ptrdiff_t exitCode=assign3DC(minCostTuples,
                                                       &gain,
                                                          &q,
                                                           u,
                                                   tempSpace,
                                                       nDims,
                                                       CCopy,
                                                    maximize,
                                               subgradMethod,
                                                     maxIter,
                                                      epsVal,
                                                      epsVal,
                                                      param1,
                                                      param2,
                                                      param3);
            
            if(exitCode==-3||exitCode==-1) {
                //If no valid assignment was found, then we just start the
                //algorithm without initialization. This most commonly
                //occurs when the problem is infeasible and no solution
                //will be found in the end anyway.
                q=-std::numeric_limits<double>::infinity();
                gain=std::numeric_limits<double>::infinity();
            }
            break;
        }
    }

    //If the initial approximation returned a suboptimal solution.
    if(fabs(gain-q)>=epsVal*gain) {
        const size_t numLevels=n1;//The number of things to assign.

        //There are n1*n2*n3 total tuples. However, we define a level by
        //the value of the first index of the tuple. In the initial level,
        //there are n2*n3 possible tuples (fix the value of the first
        //index of C and the other two are free). In subsequent levels,
        //there are (n2-curLevel+1)*(n3-curLevel+1) possibilities.
        const size_t maxTuplesPerLevel=n2*n3;
        //Generate all possible tuples of the second and third coordinates
        //that can be assigned and put them into origFreeTuples. The
        //assignment of the first coordinate is the level of the recursion
        //below.
        size_t *origFreeTuples;
        {
            const size_t maxVals[]={n2-1,n3-1};
            origFreeTuples=reinterpret_cast<size_t*>(curBuffPtr);
            curBuffPtr+=2*maxTuplesPerLevel*sizeof(size_t);
            
            genAllTuplesCPP(origFreeTuples, 2, maxVals, false);
        }

        //The list of tuples that are candidates in the current level.
        //freeTuples=zeros(2,maxTuplesPerLevel,numLevels);
        size_t *freeTuples=reinterpret_cast<size_t*>(curBuffPtr);
        curBuffPtr+=2*maxTuplesPerLevel*numLevels*sizeof(size_t);
        
        //Save the tuples that could possibly be assigned in the first
        //level.
        //freeTuples(:,:,1)=origFreeTuples;
        std::copy(origFreeTuples,origFreeTuples+2*maxTuplesPerLevel,freeTuples);
        
        //These values link the modified tuples at each level to the
        //original tuples. These are important, because the ordering of the
        //tuples in the lower levels will be changed.
        //freeTupleOrigIdx=zeros(maxTuplesPerLevel,numLevels);
        size_t *freeTupleOrigIdx=reinterpret_cast<size_t*>(curBuffPtr);
        curBuffPtr+=maxTuplesPerLevel*numLevels*sizeof(size_t);

        //freeTupleOrigIdx(:,1)=1:maxTuplesPerLevel in Matlab. Indices are
        //from 0 here.
        for(size_t i=0;i<maxTuplesPerLevel;i++) {
            freeTupleOrigIdx[i]=i;
        }
        
        //This keeps track of how many tuples are being considered at each
        //level after eliminating those that conflict with prior
        //assignments and eliminating those whose lower bounds are too
        //high.
        //numFreeTuples=zeros(numLevels,1);
        size_t *numFreeTuples=reinterpret_cast<size_t*>(curBuffPtr);
        curBuffPtr+=numLevels*sizeof(size_t);
        
        //In the top level, all possibilities are initially considered.
        numFreeTuples[0]=maxTuplesPerLevel;

        //The list of tuples that were candidates when entering the current
        //level. In the current level, some tuples might be removed due to
        //a low cost at that level. However, when going to a higher level,
        //all tuples that do not conflict with assignments at lower levels
        //must be considered. Thus, this stores the value of the tuples in
        //the current level before pruning.
        //freeTuplesEntering=zeros(2,maxTuplesPerLevel,numLevels);
        size_t *freeTuplesEntering=reinterpret_cast<size_t*>(curBuffPtr);
        curBuffPtr+=2*maxTuplesPerLevel*numLevels*sizeof(size_t);

        //freeTuplesEntering(:,:,1)=origFreeTuples;
        std::copy(origFreeTuples,origFreeTuples+2*maxTuplesPerLevel,freeTuplesEntering);

        //Like freeTupleOrigIdx, these link the values in
        //freeTuplesEntering to the original set of tuples.
        //freeTuplesEnteringOrigIdx=zeros(maxTuplesPerLevel,numLevels);
        size_t *freeTuplesEnteringOrigIdx=reinterpret_cast<size_t*>(curBuffPtr);
        curBuffPtr+=maxTuplesPerLevel*numLevels*sizeof(size_t);

        //freeTuplesEnteringOrigIdx(:,1)=1:maxTuplesPerLevel; -but the
        //indexation has been changed to 0.
        for(size_t i=0;i<maxTuplesPerLevel;i++) {
            freeTuplesEnteringOrigIdx[i]=i;
        }

        //This keeps track of how many tuples are being considered at each
        //level after eliminating those that conflict with prior
        //assignments but not eliminating anything where the lower bounds
        //are too high.
        //numFreeTuplesEntering=zeros(numLevels,1);
        size_t *numFreeTuplesEntering=reinterpret_cast<size_t*>(curBuffPtr);
        curBuffPtr+=numLevels*sizeof(size_t);

        //In the top level, all possibilities are considered.
        numFreeTuplesEntering[0]=maxTuplesPerLevel;

        //The tuple that is assigned at each level, after the indexation
        //has been changed to relate it to the current cost submatrix at
        //that level. We don't need to save the first index here, because
        //it is just the level number. We don't need to save the assigned
        //tuple in the maximum level, because we use it right away.
        //assignedTuples=zeros(2,numLevels-1);
        size_t *assignedTuples=reinterpret_cast<size_t*>(curBuffPtr);
        curBuffPtr+=2*(numLevels-1)*sizeof(size_t);

        //This links the assigned tuple to the index of the tuple in
        //origFreeTuples. It is simpler to use this than to try to
        //reconstruct the tuple from assignedTuples, where each tuple has
        //been modified based on assignments at lower levels. Unlike
        //assignedTuples, this needs to be numLevels in size.
        //assignedTupleOrigIdx=zeros(numLevels,1);
        size_t *assignedTupleOrigIdx=reinterpret_cast<size_t*>(curBuffPtr);
        curBuffPtr+=numLevels*sizeof(size_t);

        //The cumulative assigned cost at each level. The final cost is
        //computed at the top level and does not need to be stored in this.
        //cumAssignedCost=zeros(numLevels-1,1);
        double *cumAssignedCost=reinterpret_cast<double*>(curBuffPtr);
        curBuffPtr+=(numLevels-1)*sizeof(double);

        //Allocate space for the bounded total costs for assigning
        //something at each level. Each level is given by the value in the
        //first index that is assigned. For each first index, there are
        //n2*n3 values. There are no bounds for the deepest level, because
        //the cost matrix can be directly used.
        //boundVals=zeros(maxTuplesPerLevel,numLevels-1);
        double *boundVals=reinterpret_cast<double*>(curBuffPtr);
        curBuffPtr+=maxTuplesPerLevel*(numLevels-1)*sizeof(double);
        
        //This holds what was the minimum cost value before going to a
        //higher level of recursion. Thus, when returning to a particular
        //level, if the minimum cost value changed, one should check the
        //previous bounds and see if any of the candidate tuples should be
        //thrown out and never visited.
        //prevMinCost=zeros(numLevels-1,1);
        double *prevMinCost=reinterpret_cast<double*>(curBuffPtr);
        curBuffPtr+=(numLevels-1)*sizeof(double);

        //Some of the lower bound functions require temporary buffer space
        //to compute the results. Here, we allocate the space, if
        //necessary.
        void *LBTempBuffer=NULL;
        switch(boundType) {
            case 0:
            {
                const size_t bufferSize=assign3DLBHungarianBufferSize(nDims);
                LBTempBuffer=reinterpret_cast<void*>(curBuffPtr);
                curBuffPtr+=bufferSize;
                break;
            }
            case 1://This method does not require a buffer to compute the
                   //lower bounds.
                break;
            default://boundType=2
            {
                const size_t bufferSize=assign3DLBDual0BufferSize(nDims);
                LBTempBuffer=reinterpret_cast<void*>(curBuffPtr);
                curBuffPtr+=bufferSize;
            }
        }

        //Allocate space for the marginal costs at each level (the cost
        //submatrices). In Matlab one might this like
        //CostMatCur=zeros(n1,n2,n3,numLevels); and then easily address
        //each level with the final index. However, that wastes a lot of
        //space, because each level removes one element from every
        //dimension. The first matrix is n1Xn2Xn3 in size, the next is
        //(n1-1)*(n2-1)*(n3-1) in size and so on until the final level
        //which is 1X(n2-n1+1)X(n3-n1+1) in size. To save space, we shall
        //allocate the minimum amount of space required, as well as an
        //offset array, so that we know where the elements of the matrices
        //start and end so that we can later consider each matrix
        //individually.
        //costMatCurStartIdx=zeros(numLevels+1,1);
        size_t *costMatCurStartIdx=reinterpret_cast<size_t*>(curBuffPtr);
        curBuffPtr+=(numLevels+1)*sizeof(size_t);

        costMatCurStartIdx[0]=0;
        for(size_t curLevel=1;curLevel<numLevels+1;curLevel++) {
            const size_t offset=curLevel-1;
            costMatCurStartIdx[curLevel]=costMatCurStartIdx[curLevel-1]+(n1-offset)*(n2-offset)*(n3-offset);
        }
        //The final value in costMatCurStartIdx is the number of elements
        //in the matrix.
        //CostMatCur=zeros(costMatCurStartIdx(numLevels+1)-1,1);
        double *CostMatCur=reinterpret_cast<double*>(curBuffPtr);
        curBuffPtr+=costMatCurStartIdx[numLevels]*sizeof(double);

        //Assign the marginal cost matrix of the first level.
        //CostMatCur(costMatCurStartIdx(1):(costMatCurStartIdx(2)-1))=C(:);
        std::copy(C,C+n1*n2*n3,CostMatCur);
        
        //When sorting, we store the indices of the sorted items here so
        //that we can sort other things with the same ordering.
        size_t *idxVec=reinterpret_cast<size_t*>(curBuffPtr);
        curBuffPtr+=maxTuplesPerLevel*sizeof(size_t);
     
        //Enter the recursion (unrolled).
        ptrdiff_t curLevel=0;
        bool increasingLevelNumber=true;
        while(curLevel>=0) {            
            if(curLevel==static_cast<ptrdiff_t>(numLevels)-1) {
                //If we have reached the top level, the final assignment is
                //just choosing the most likely free tuple in the current
                //cost matrix.

                double minVal=std::numeric_limits<double>::infinity();
                size_t minIdx=0;
                //for curTuple=1:numFreeTuples(curLevel)
                for(size_t curTuple=0;curTuple<numFreeTuples[curLevel];curTuple++) {
                    //i=1;
                    const size_t i=0;
                    //j=freeTuples(1,curTuple,curLevel);
                    const size_t j=freeTuples[sub2Idx3D(0,curTuple,curLevel,2,maxTuplesPerLevel)];
                    //k=freeTuples(2,curTuple,curLevel);
                    const size_t k=freeTuples[sub2Idx3D(1,curTuple,curLevel,2,maxTuplesPerLevel)];
                    //The index within the cost matrix at the current level.
                    //idx=sub2ind(dims,i,j,k);
                    const size_t idx=sub2Idx3D(i,j,k,n1-curLevel,n2-curLevel);
                    const double curCost=CostMatCur[costMatCurStartIdx[curLevel]+idx];

                    if(curCost<minVal) {
                        minVal=curCost;
                        minIdx=curTuple;
                    }
                }
 
                if(numLevels>0) {
                    minVal+=cumAssignedCost[curLevel-1];
                }

                if(minVal<gain) {//If a new global minimum has been found.
                    gain=minVal;
                    
                    //assignedTupleOrigIdx(curLevel)=freeTupleOrigIdx(minIdx,curLevel);
                    assignedTupleOrigIdx[curLevel]=freeTupleOrigIdx[sub2Idx2D(minIdx,curLevel,maxTuplesPerLevel)];
                    
                    for(size_t i=0;i<numLevels;i++) {
                        ptrdiff_t * const curMinCostTuple=minCostTuples+i*3;
                        const size_t * const curOrigFreeTup=origFreeTuples+assignedTupleOrigIdx[i]*2;
                        //minCostTuples(:,i)=[i;origFreeTuples(:,assignedTupleOrigIdx(i))];
                        curMinCostTuple[0]=i;
                        curMinCostTuple[1]=curOrigFreeTup[0];
                        curMinCostTuple[2]=curOrigFreeTup[1];
                    }
                }
                
                increasingLevelNumber=false;
                curLevel--;
                continue;
            } else if(increasingLevelNumber==true) {
                //We have just entered this level, so we have to first
                //compute the lower bounds for the current level.

                size_t curTuple=0;
                while(curTuple<numFreeTuples[curLevel]) {
                    //i=1;
                    const size_t i=0;
                    //j=freeTuples(1,curTuple,curLevel);
                    const size_t j=freeTuples[sub2Idx3D(0,curTuple,curLevel,2,maxTuplesPerLevel)];
                    //k=freeTuples(2,curTuple,curLevel);
                    const size_t k=freeTuples[sub2Idx3D(1,curTuple,curLevel,2,maxTuplesPerLevel)];
                                        
                    //The i,j,k indices of the tuple are with respect to
                    //indices in the shrunken cost matrix
                    //CostMatCur(i,j,k,curLevel), not to the set of tuples
                    //with respect to the global cost matrix C. Since we
                    //are picking off rows one at a time, i will always be
                    //0.
                    
                    //The dimensionality of the cost matrix at the current
                    //level is (n1-curLevel)X(n2-curLevel)X(n3-curLevel)
                    //The index within the cost matrix at the current
                    //level.
                    const size_t idx=sub2Idx3D(i,j,k,n1-curLevel,n2-curLevel);
                    const double curCost=CostMatCur[costMatCurStartIdx[curLevel]+idx];
                    
                    //If the assignment is not allowed.
                    if(!std::isfinite(curCost)) {
                        curTuple++;
                        continue;
                    }

                    double newBoundVal;
                    //Create the cost matrix at the next level up by
                    //removing the appropriate entry from each dimension of
                    //the current cost matrix. 
                    //CCur=reshape(CostMatCur(costMatCurStartIdx(curLevel):(costMatCurStartIdx(curLevel+1)-1)),[n1-curLevel+1,n2-curLevel+1,n3-curLevel+1]);
                    //The indices of the cost matrix one level up.
                    //spanNext=costMatCurStartIdx(curLevel+1):(costMatCurStartIdx(curLevel+2)-1);
                    //CostMatCur(spanNext)=CCur(2:(n1-curLevel+1),[1:(j-1),(j+1):(n2-curLevel+1)],[1:(k-1),(k+1):(n3-curLevel+1)]);
                    {
                        const double * const sourcePtr=CostMatCur+costMatCurStartIdx[curLevel];
                        double *const destPtr=CostMatCur+costMatCurStartIdx[curLevel+1];
                        const size_t skipIdx[]={i,j,k};
                        const size_t sourceMatDims[]={n1-curLevel,n2-curLevel,n3-curLevel};
                        copyMat3DOmitDims(destPtr, sourceMatDims, sourcePtr,skipIdx);

                        //Evaluate the correct bound.
                        double curBoundVal;
                        const size_t destMatDims[]={n1-curLevel-1,n2-curLevel-1,n3-curLevel-1};
                        switch(boundType) {
                            case 0:
                                curBoundVal=assign3DLBHungarian(destMatDims,destPtr,LBTempBuffer);
                                break;
                            case 1:
                                curBoundVal=assign3DLBCPierskalla(destMatDims,destPtr);
                                break;
                            default://boundType=2
                                curBoundVal=assign3DLBDual0(destMatDims,destPtr,LBTempBuffer);
                                break;
                        }

                        if(curLevel>0) {
                            newBoundVal=cumAssignedCost[curLevel-1]+curCost+curBoundVal;
                        } else {
                            newBoundVal=curCost+curBoundVal;
                        }

                        boundVals[sub2Idx2D(curTuple,curLevel,maxTuplesPerLevel)]=newBoundVal;
                    }

                    //If the bound is no better than the current best
                    //solution, then remove that tuple as a contender in
                    //this level. We just overwrite it with the last tuple,
                    //decriment the number of tuples present and then
                    //revisit the tuple with the same index (which will now
                    //be the last tuple).
                    if(newBoundVal>=gain) {
                        CostMatCur[costMatCurStartIdx[curLevel]+idx]=std::numeric_limits<double>::infinity();
                        //freeTuples(:,curTuple,curLevel)=freeTuples(:,numFreeTuples(curLevel),curLevel);
                        const size_t *freeTupSource=freeTuples+sub2Idx3D(0,numFreeTuples[curLevel]-1,curLevel,2,maxTuplesPerLevel);
                        size_t *const freeTupDest=freeTuples+sub2Idx3D(0,curTuple,curLevel,2,maxTuplesPerLevel);
                        std::copy(freeTupSource,freeTupSource+2,freeTupDest);
                        
                        //freeTupleOrigIdx(curTuple,curLevel)=freeTupleOrigIdx(numFreeTuples(curLevel),curLevel);
                        freeTupleOrigIdx[sub2Idx2D(curTuple,curLevel,maxTuplesPerLevel)]=freeTupleOrigIdx[sub2Idx2D(numFreeTuples[curLevel]-1,curLevel,maxTuplesPerLevel)];
                        
                        numFreeTuples[curLevel]--;
                    } else {
                        curTuple++;
                    }
                }

                //If there are not enough tuples left to make a full
                //assignment. The assumption here is that numLevels>1, so
                //if we are here, there must be one tuple left to assign
                //after assigning the one in this level.
                if(numFreeTuples[curLevel]==0) {
                    increasingLevelNumber=false;
                    curLevel--;
                    continue;
                }
                
                //Sort the tuples by cost.
                //[boundVals(1:numFreeTuples(curLevel),curLevel),idxVec]=sort(boundVals(1:numFreeTuples(curLevel),curLevel),'descend');
                heapSortVecIdxCDouble(numFreeTuples[curLevel], boundVals+curLevel*maxTuplesPerLevel, idxVec,1);
                //freeTuples(:,1:numFreeTuples(curLevel),curLevel)=freeTuples(:,idxVec,curLevel);
                permuteMatrixColsCPP(2,numFreeTuples[curLevel],freeTuples+2*maxTuplesPerLevel*curLevel,idxVec);
                //freeTupleOrigIdx(1:numFreeTuples(curLevel),curLevel)=freeTupleOrigIdx(idxVec,curLevel);    
                permuteVectorCPP(numFreeTuples[curLevel], freeTupleOrigIdx+curLevel*maxTuplesPerLevel,idxVec);            
                
                //Record the current minimum cost value.
                prevMinCost[curLevel]=gain;
            }

            //If we are here, then either increasingLevelNumber=false and
            //we are backtracking, or increasingLevelNumber=true, and we
            //have just entered this level.
            if(increasingLevelNumber==false) {
                //If we are backtracking, then we remove the last tuple
                //visited in this level from consideration.
                //i=1;
                const size_t i=0;
                //j=assignedTuples(1,curLevel);
                const size_t j=assignedTuples[sub2Idx2D(0,curLevel,2)];
                //k=assignedTuples(2,curLevel);
                const size_t k=assignedTuples[sub2Idx2D(1,curLevel,2)];

                //The index within the cost matrix at the current
                //level.
                //idx=sub2ind(dims,i,j,k);
                const size_t idx=sub2Idx3D(i,j,k,n1-curLevel,n2-curLevel);
                //freeTuples(:,curLevel,numFreeTuples(curLevel),curLevel)
                //contains the last assigned tuple for this level. Here, we
                //eliminate it from further consideration by setting
                //CostMatCur(i,j,k,curLevel) to Inf and reducing the count
                //of the tuple by one.
                CostMatCur[costMatCurStartIdx[curLevel]+idx]=std::numeric_limits<double>::infinity();
                numFreeTuples[curLevel]--;
                
                //If the last tuple visited leaves too few tuples left for
                //a full assignment, then go up another level.
                if(numFreeTuples[curLevel]==0) {
                    increasingLevelNumber=false;
                    curLevel--;
                    continue;
                }

                //Also, if the minimum bound changed, remove all tuples
                //whose bounds are larger than the new minimum.
                if(prevMinCost[curLevel]!=gain) {
                    prevMinCost[curLevel]=gain;

                    //Find all of the tuples where the bound is no better
                    //than the current best solution and remove them as
                    //contenders in this level.

                    //Remember that when entering this level, the largest
                    //bounds are at the beginning of the boundVals array,
                    //because it has been sorted in decreasing order. Thus,
                    //we have to find the first index where the bound is
                    //<gain.
                    //idx=find(boundVals(1:numFreeTuples(curLevel),curLevel)<gain,1);
                    size_t idx=0;
                    double *curBoundVal=boundVals+curLevel*maxTuplesPerLevel;
                    for(idx=0;idx<numFreeTuples[curLevel];idx++) {
                        if(curBoundVal[idx]<gain) {
                            break;
                        }
                    }

                    //If this gets rid of all tuples.
                    //if(isempty(idx))
                    if(idx==numFreeTuples[curLevel]) {
                        curLevel--;
                        continue;
                    } else if(idx>0) {
                        //If here, remove the tuples with bounds that are
                        //too large.
                        //sel1=1:(numFreeTuples(curLevel)-idx+1);
                        //sel2=idx:numFreeTuples(curLevel);
                        //boundVals(sel1,curLevel)=boundVals(sel2,curLevel);
                        {
                            const double *startIdx=boundVals+sub2Idx2D(idx,curLevel,maxTuplesPerLevel);
                            const double *endIdx=boundVals+sub2Idx2D(numFreeTuples[curLevel],curLevel,maxTuplesPerLevel);
                            double *destIdx=boundVals+sub2Idx2D(0,curLevel,maxTuplesPerLevel);
                            std::copy(startIdx,endIdx,destIdx);
                        }

                        //Set the entries in the cost matrix associated
                        //with the tuples that are being removed to Inf, so
                        //they won't be assigned.
                        //for curTuple=1:(idx-1)
                        for(size_t curTuple=0;curTuple<idx-1;curTuple++) {
                            //i=1;
                            const size_t i=0;
                            //j=freeTuples(1,curTuple,curLevel);
                            const size_t j=freeTuples[sub2Idx3D(0,curTuple,curLevel,2,maxTuplesPerLevel)];
                            //k=freeTuples(2,curTuple,curLevel);
                            const size_t k=freeTuples[sub2Idx3D(1,curTuple,curLevel,2,maxTuplesPerLevel)];

                            //The dimensionality of the cost matrix at the
                            //current level is
                            //(n1-curLevel)X(n2-curLevel)X(n3-curLevel)
                            //The index within the cost matrix at the
                            //current level.
                            const size_t idx=sub2Idx3D(i,j,k,n1-curLevel,n2-curLevel);
                            CostMatCur[costMatCurStartIdx[curLevel]+idx]=std::numeric_limits<double>::infinity();
                        }

                        //sel1=1:(numFreeTuples(curLevel)-idx+1);
                        //sel2=idx:numFreeTuples(curLevel);
                        //freeTuples(:,sel1,curLevel)=freeTuples(:,sel2,curLevel);
                        {
                           const size_t *startIdx=freeTuples+sub2Idx3D(0,idx,curLevel,2,maxTuplesPerLevel);
                           const size_t *endIdx=freeTuples+sub2Idx3D(0,numFreeTuples[curLevel],curLevel,2,maxTuplesPerLevel);
                           size_t *destIdx=freeTuples+sub2Idx3D(0,0,curLevel,2,maxTuplesPerLevel);
                           std::copy(startIdx,endIdx,destIdx);
                        }
                        
                        //sel1=1:(numFreeTuples(curLevel)-idx+1);
                        //sel2=idx:numFreeTuples(curLevel);
                        //freeTupleOrigIdx(sel1,curLevel)=freeTupleOrigIdx(sel2,curLevel);
                        {
                            const size_t *startIdx=freeTupleOrigIdx+sub2Idx2D(idx,curLevel,maxTuplesPerLevel);
                            const size_t *endIdx=freeTupleOrigIdx+sub2Idx2D(numFreeTuples[curLevel],curLevel,maxTuplesPerLevel);
                            size_t *destIdx=freeTupleOrigIdx+sub2Idx2D(0,curLevel,maxTuplesPerLevel);
                            std::copy(startIdx,endIdx,destIdx);
                        }

                        numFreeTuples[curLevel]-=idx;
                    }
                }

                increasingLevelNumber=true;
            }
            
            //The tuple having the lowest cost is visited first.
            //Assign the tuple in this current level.
            //assignedTuples(:,curLevel)=freeTuples(:,numFreeTuples(curLevel),curLevel);
            {
                const size_t *const sourceIdx=freeTuples+sub2Idx3D(0,numFreeTuples[curLevel]-1,curLevel,2,maxTuplesPerLevel);
                size_t *const destIdx=assignedTuples+sub2Idx2D(0,curLevel,2);
                destIdx[0]=sourceIdx[0];
                destIdx[1]=sourceIdx[1];
            }
        
            //assignedTupleOrigIdx(curLevel)=freeTupleOrigIdx(numFreeTuples(curLevel),curLevel);
            assignedTupleOrigIdx[curLevel]=freeTupleOrigIdx[sub2Idx2D(numFreeTuples[curLevel]-1,curLevel,maxTuplesPerLevel)];
 
            //Add the cost of the assigned tuple to the cumulative assigned
            //cost.
            //i=1;
            const size_t i=0;
            //j=assignedTuples(1,curLevel);
            const size_t j=assignedTuples[sub2Idx2D(0,curLevel,2)];
            //k=assignedTuples(2,curLevel);
            const size_t k=assignedTuples[sub2Idx2D(1,curLevel,2)];
                       
            //The dimensionality of the cost matrix at the current
            //level is (n1-curLevel)X(n2-curLevel)X(n3-curLevel)
            //The index within the cost matrix at the current
            //level.
            const size_t idx=sub2Idx3D(i,j,k,n1-curLevel,n2-curLevel);
            const double curCost=CostMatCur[costMatCurStartIdx[curLevel]+idx];
        
            if(curLevel>0) {
                cumAssignedCost[curLevel]=cumAssignedCost[curLevel-1]+curCost;
            } else {
                cumAssignedCost[curLevel]=curCost;
            }

            //Construct the cost matrix for the next level. This means
            //removing the entire row, column, etc. in each dimension that
            //contains the assigned tuple.
            //CCur=reshape(CostMatCur(costMatCurStartIdx(curLevel):(costMatCurStartIdx(curLevel+1)-1)),[n1-curLevel+1,n2-curLevel+1,n3-curLevel+1]);
            //The indices of the cost matrix one level up.
            //spanNext=costMatCurStartIdx(curLevel+1):(costMatCurStartIdx(curLevel+2)-1);
            //CostMatCur(spanNext)=CCur(2:(n1-curLevel+1),[1:(j-1),(j+1):(n2-curLevel+1)],[1:(k-1),(k+1):(n3-curLevel+1)]);
            {
                const double * const sourcePtr=CostMatCur+costMatCurStartIdx[curLevel];
                double *const destPtr=CostMatCur+costMatCurStartIdx[curLevel+1];
                const size_t skipIdx[]={i,j,k};
                const size_t sourceMatDims[]={n1-curLevel,n2-curLevel,n3-curLevel};
                copyMat3DOmitDims(destPtr, sourceMatDims, sourcePtr,skipIdx);
            }

            //Copy the tuples that can be considered for the next and
            //subsequent levels into freeTuples(:,:,curLevel+1) and
            //freeTuplesEntering(:,:,curLevel+1).
            //This draws tuples for the next level down from
            //freeTuplesEntering(:,:,curLevel), because
            //freeTuples(:,:,curLevel) may have had tuples removed only for
            //curLevel due to bounding. At the same time as the tuples are
            //added, the coordinates of all valid tuples are adjusted to
            //deal with the removed elements of CostMatCur.

            numFreeTuples[curLevel+1]=0;
            //for curTuple=1:numFreeTuplesEntering(curLevel)
            for(size_t curTuple=0;curTuple<numFreeTuplesEntering[curLevel];curTuple++) {
                //if(any(freeTuplesEntering(:,curTuple,curLevel)==assignedTuples(:,curLevel)))
                if(anyArralElsEqual(freeTuplesEntering+sub2Idx3D(0,curTuple,curLevel,2,maxTuplesPerLevel), assignedTuples+curLevel*2,2)) {
                    continue;
                } else {
                    numFreeTuples[curLevel+1]++;
                    
                    //Free tuples are indexed with respect to the shrunken
                    //matrix of costs in CostMatCur(:,:,:,curLevel+1). This
                    //means that for index values > the assigned index, one
                    //must be subtracted.
                    //freeTuples(:,numFreeTuples(curLevel+1),curLevel+1)=freeTuplesEntering(:,curTuple,curLevel)-(freeTuplesEntering(:,curTuple,curLevel)>assignedTuples(:,curLevel));
                    {
                        const size_t *sourcePtr=freeTuplesEntering+sub2Idx3D(0,curTuple,curLevel,2,maxTuplesPerLevel);
                        size_t *destPtr=freeTuples+sub2Idx3D(0,numFreeTuples[curLevel+1]-1,curLevel+1,2,maxTuplesPerLevel);
                        const bool boolVal1=sourcePtr[0]>assignedTuples[sub2Idx2D(0,curLevel,2)];
                        const bool boolVal2=sourcePtr[1]>assignedTuples[sub2Idx2D(1,curLevel,2)];
   
                        destPtr[0]=sourcePtr[0]-boolVal1;
                        destPtr[1]=sourcePtr[1]-boolVal2;
                    }

                    //freeTupleOrigIdx(numFreeTuples(curLevel+1),curLevel+1)=freeTuplesEnteringOrigIdx(curTuple,curLevel);
                    freeTupleOrigIdx[sub2Idx2D(numFreeTuples[curLevel+1]-1,curLevel+1,maxTuplesPerLevel)]=freeTuplesEnteringOrigIdx[sub2Idx2D(curTuple,curLevel,maxTuplesPerLevel)];
                }
            }

            numFreeTuplesEntering[curLevel+1]=numFreeTuples[curLevel+1];

           //freeTuplesEntering(:,1:numFreeTuples(curLevel+1),curLevel+1)=freeTuples(:,1:numFreeTuples(curLevel+1),curLevel+1);
           //freeTuplesEnteringOrigIdx(1:numFreeTuples(curLevel+1),curLevel+1)=freeTupleOrigIdx(1:numFreeTuples(curLevel+1),curLevel+1);
            {
                const size_t *sourceTup=freeTuples+sub2Idx3D(0,0,curLevel+1,2,maxTuplesPerLevel);
                size_t *destTup=freeTuplesEntering+sub2Idx3D(0,0,curLevel+1,2,maxTuplesPerLevel);
                const size_t * const sourceIdx=freeTupleOrigIdx+sub2Idx2D(0,curLevel+1,maxTuplesPerLevel);
                size_t * const destIdx=freeTuplesEnteringOrigIdx+sub2Idx2D(0,curLevel+1,maxTuplesPerLevel);
                        
                for(size_t curFreeTuple=0;curFreeTuple<numFreeTuples[curLevel+1];curFreeTuple++ ) {
                    destTup[0]=sourceTup[0];
                    destTup[1]=sourceTup[1];
                    destIdx[curFreeTuple]=sourceIdx[curFreeTuple];
                    
                    //Go to the next tuple.
                    destTup+=2;
                    sourceTup+=2;
                }   
            }
            //Go to the next level.
            curLevel++;
        }
    }

    //Adjust the gain for the initial offset of the cost matrix.
    if(maximize==true) {
        gain=-gain+CDelta*n1;
    } else {
        gain=gain+CDelta*n1;
    }
    
    return gain;
}

/*LICENSE:
*
*The source code is in the public domain and not licensed or under
*copyright. The information and software may be used freely by the public.
*As required by 17 U.S.C. 403, third parties producing copyrighted works
*consisting predominantly of the material produced by U.S. government
*agencies must provide notice with such work(s) identifying the U.S.
*Government material incorporated and stating that such material is not
*subject to copyright protection.
*
*Derived works shall not identify themselves in a manner that implies an
*endorsement by or an affiliation with the Naval Research Laboratory.
*
*RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
*SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
*RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
*OF RECIPIENT IN THE USE OF THE SOFTWARE.
*/
