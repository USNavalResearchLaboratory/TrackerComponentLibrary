/*RUNNALLSGAUSSMIXRED Perform Gaussian mixture reduction using the greedy
*    merging algorithm by Runnals in [1]. See the comments to that Matlab
*   language implementation for more information.
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
* The algorithm is run in Matlab using the command format
* [w,mu,P,numPerCluster,clusterInfo,minCostClustPartition,wAll,muAll,PAll,clusterInfoAll,numPerClusterAll]=RunnalsGaussMixRed(w,mu,P,K,gammaBound,KMax);
*
*January 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4514 )
#endif

#include "mex.h"

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#include "MexValidation.h"

//The pragmas get rid of all warning in Visual Studio that arise in the
//Eigen library. There are tons of them.
#ifdef _MSC_VER
#pragma warning(push, 0)
#pragma warning( disable : 4710 4711 5045 )
#endif
#include "Eigen/Dense"
#ifdef _MSC_VER
#pragma warning(pop)
//Get rid of inlining warnings.
#pragma warning( disable : 4711 )
#endif

#include <limits.h>
#include <cmath>
#include <algorithm>

//Prototypes.
size_t BDistBufferSize(const size_t xDim) {
    //The size of the P12 matrix and the diff vector.
    return sizeof(double)*(xDim*xDim+xDim);
}
double BDist(const size_t N,const double &w1,const double &w2,const double * const mu1,const double * const mu2,const double * const P1,const double * const P2,const double wLogDetP1,const double wLogDetP2, double *buffer);
void RunnalsGaussMixRedCPP(const size_t xDim,const size_t N, double *w, double *mu, double *P, const size_t K, const double gammaBound, const size_t KMax, size_t &KCur, size_t *numPerCluster,size_t *clusterInfo,size_t &maxInClust,void *tempBuffer,double *wAll,double *muAll,double *PAll,size_t *clusterInfoAll,size_t *numPerClusterAll,size_t &curAllIdx);
size_t RunnalsGaussMixRedBufferSize(const size_t N, const size_t xDim) {
    size_t buffSize=BDistBufferSize(xDim);
    buffSize+=sizeof(double)*xDim;//For diff1
    buffSize+=sizeof(double)*xDim;//For diff2
    buffSize+=sizeof(double)*N*N;//For M
    buffSize+=sizeof(double)*N;//For wLogDetPVals
    buffSize+=sizeof(bool)*N;//For selIdxPresent

    return buffSize;
}

void getSizes4ReturningSteps(const size_t N,const size_t K,const size_t KMax, size_t &maxComp, size_t &totalRedHyps) {
    maxComp=std::min(KMax,N);
    if(K>maxComp) {
        totalRedHyps=1;
        maxComp=N;
    } else {
        totalRedHyps=maxComp-K+1;
    }
    //The matrices related to these are
    // wAllSize is maxCompXtotalRedHyps
    // muAllSize is xDimXmaxCompXtotalRedHyps
    // PAllSize is xDimXxDimXmaxCompXtotalRedHyps
    // clusterInfoAll is NXmaxCompXtotalRedHyps
    // numPerClusterAll is NXtotalRedHyps
}

void clusterInfo2MinCostClustPartition(const size_t KRed, const size_t maxInClust, const size_t *numPerCluster,const size_t * const clusterInfo,size_t *minCostClustPartition) {

    for(size_t curCluster=0;curCluster<KRed;curCluster++ ) {
        const size_t numCur=numPerCluster[curCluster];
        const size_t *const measIdxList=clusterInfo+maxInClust*curCluster;

        for(size_t k=0;k<numCur;k++) {
            const size_t curIdx=measIdxList[k];
            minCostClustPartition[curIdx]=curCluster;
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double gammaBound=std::numeric_limits<double>::infinity();
    size_t KMax=SIZE_MAX;
    mxArray *wMat=NULL;

    if(nrhs<4){
        mexErrMsgTxt("Not enough inputs.");
        return;
    }

    if(nrhs>6) {
        mexErrMsgTxt("Too many inputs.");
        return;
    }

    if(nlhs>11) {
        mexErrMsgTxt("Too many outputs.");
    }

    const size_t K=getSizeTFromMatlab(prhs[3]);
    if(nrhs>=5&&!mxIsEmpty(prhs[4])) {
        gammaBound=getDoubleFromMatlab(prhs[4]);
    }
    if(nrhs>=6&&!mxIsEmpty(prhs[5])) {
        KMax=getSizeTFromMatlab(prhs[5]);
    }

    const size_t xDim=mxGetM(prhs[1]);
    const size_t N=mxGetN(prhs[1]);
    const size_t N2=N*N;

    //Verify the sizes of w, mu and P
    if(mxIsEmpty(prhs[0])) {
        if(!mxIsEmpty(prhs[1])&&!mxIsEmpty(prhs[2])) {
            //Create a uniformly weighted w array.
            wMat=mxCreateDoubleMatrix(1,N,mxREAL);
            double *wTemp=mxGetDoubles(wMat);
            double NInvFloat=1.0/static_cast<double>(N);

            for(size_t k=0;k<N;k++) {
                wTemp[k]=NInvFloat;
            }
        } else {
            //All 3 are empty, so we will just return a bunch of empty matrices.
            plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
            for(size_t k=1;k<static_cast<size_t>(nlhs);k++) {
                plhs[k]=mxCreateDoubleMatrix(0,0,mxREAL);
            }
        }
    } else {
        checkRealDoubleArray(prhs[0]);//w
    }

    checkRealDoubleArray(prhs[1]);//mu
    checkRealDoubleHypermatrix(prhs[2]);//P

    if(wMat==NULL&&mxGetNumberOfDimensions(prhs[0])!=2) {
        mexErrMsgTxt("The dimensions of w are incorrect.");
    } else if(mxGetNumberOfDimensions(prhs[1])!=2) {
        mexErrMsgTxt("The dimensions of mu are incorrect.");
    } else if(mxGetNumberOfDimensions(prhs[2])!=3&&!(mxGetNumberOfDimensions(prhs[2])==2&&N==1)) {
        mexErrMsgTxt("The dimensions of P are incorrect.");
    }

    if(wMat==NULL&&(mxGetNumberOfElements(prhs[0])!=N||(mxGetM(prhs[0])!=N&&mxGetN(prhs[0])!=N))) {
        mexErrMsgTxt("The dimensions of w are inconsistent.");
    }

    {
        const mwSize *dims=mxGetDimensions(prhs[2]);
        if(dims[0]!=xDim||dims[1]!=xDim) {
            mexErrMsgTxt("The dimensions of P are inconsistent.");
        }
    }

    if(N>1) {
        const mwSize *dims=mxGetDimensions(prhs[2]);
        if(dims[2]!=N) {
            mexErrMsgTxt("The dimensions of P are inconsistent.");
        }
    }

    //The input arrays will be modified in place during the reduction, so
    //they must be copied.
    if(wMat==NULL) {
        wMat=mxDuplicateArray(prhs[0]);
    }
    mxArray *muMat=mxDuplicateArray(prhs[1]);
    mxArray *PMat=mxDuplicateArray(prhs[2]);
    double * w=mxGetDoubles(wMat);
    double * mu=mxGetDoubles(muMat);
    double * P=mxGetDoubles(PMat);
    const size_t bufferSize=RunnalsGaussMixRedBufferSize(N, xDim);

    uint8_t *tempBuffer=new uint8_t[bufferSize];
    if(tempBuffer==NULL) {
        mexErrMsgTxt("Error allocating temporary buffer.");
    }
   
    //If multiple steps should be saved. Initialize to NULL to get rid of
    //compiler warnings.
    mxArray *wAllMat=NULL;
    mxArray *muAllMat=NULL;
    mxArray *PAllMat=NULL;
    double *wAll, *muAll, *PAll;
    size_t *clusterInfoAll, *numPerClusterAll;
    size_t maxComp=0;//This is just initialized to get rid of a warning.
    size_t totalRedHyps;
    if(nlhs>6) {
        getSizes4ReturningSteps(N, K, KMax, maxComp, totalRedHyps);
        wAllMat=mxCreateDoubleMatrix(maxComp,totalRedHyps,mxREAL);

        {
            const mwSize nDim=3;
            const mwSize dims[]={xDim,maxComp,totalRedHyps};
            muAllMat=mxCreateNumericArray(nDim, dims,mxDOUBLE_CLASS, mxREAL);
        }
        {
            const mwSize nDim=4;
            const mwSize dims[]={xDim,xDim,maxComp,totalRedHyps};
            PAllMat=mxCreateNumericArray(nDim, dims,mxDOUBLE_CLASS, mxREAL);
        }

        wAll=mxGetDoubles(wAllMat);
        muAll=mxGetDoubles(muAllMat);
        PAll=mxGetDoubles(PAllMat);

        clusterInfoAll=new size_t[N*maxComp*totalRedHyps];
        numPerClusterAll=new size_t[N*totalRedHyps];

        //If the elements are not overwritten, the values are kept as
        //SIZE_MAX and then when returned at doubles for Matlab are set to
        //Inf (instead of having random values from unitialized memory).
        std::fill_n(clusterInfoAll,N*maxComp*totalRedHyps,SIZE_MAX);
        //If the elements are not changed, then the number in the cluster
        //will be left at 0.
        std::fill_n(numPerClusterAll,N*totalRedHyps,static_cast<size_t>(0));

        if(clusterInfoAll==NULL||numPerClusterAll==NULL) {
            mexErrMsgTxt("Error allocating temporary buffer.");
        }
    } else {
        wAll=NULL;
        muAll=NULL;
        PAll=NULL;
        clusterInfoAll=NULL;
        numPerClusterAll=NULL;
    }

    size_t *clusterInfo=new size_t[N2];
    size_t *numPerCluster=new size_t[N];
    size_t KCur, maxInClust, curAllIdx;
    RunnalsGaussMixRedCPP(xDim,N, w, mu, P, K, gammaBound, KMax, KCur, numPerCluster, clusterInfo, maxInClust, tempBuffer,
                          wAll, muAll, PAll, clusterInfoAll, numPerClusterAll,curAllIdx);
    delete[] tempBuffer;

    //Before returning the results, the arrays have to be resized.
    {
        mwSize nDim=2;
        mwSize dims[]={1,KCur,0};
        mxSetDimensions(wMat, dims, nDim);
        dims[0]=xDim;
        dims[1]=KCur;
        mxSetDimensions(muMat, dims, nDim);
        dims[0]=xDim;
        dims[1]=xDim;
        if(KCur>1) {
            nDim=3;
            dims[2]=KCur;
        }
        mxSetDimensions(PMat, dims, nDim);

        plhs[0]=wMat;
        if(nlhs>1) {
            plhs[1]=muMat;
            if(nlhs>2) {
                plhs[2]=PMat;
            } else {
                mxDestroyArray(PMat);
            }
        } else { 
            mxDestroyArray(muMat);
            mxDestroyArray(PMat);
        }
    }

    if(nlhs>=4) {
        mxArray *numPerClusterMat=mxCreateDoubleMatrix(KCur,1,mxREAL);
        double *numPerClusterData=mxGetDoubles(numPerClusterMat);
        
        for(size_t k=0;k<KCur;k++) {
            numPerClusterData[k]=static_cast<double>(numPerCluster[k]);
        }

        plhs[3]=numPerClusterMat;

        if(nlhs>=5) {
            mxArray *clusterInfoMat=mxCreateDoubleMatrix(maxInClust,KCur,mxREAL);
            double *clusterInfoData=mxGetDoubles(clusterInfoMat);

            for(size_t k=0;k<maxInClust*KCur;k++) {
                //Adding 1 changes it to Matlab's indexation from 1.
                //Unassigned values are marked with Inf.
                if(clusterInfo[k]==SIZE_MAX) {
                    clusterInfoData[k]=std::numeric_limits<double>::infinity();;
                } else {
                    clusterInfoData[k]=static_cast<double>(clusterInfo[k]+1);
                }
            }
            plhs[4]=clusterInfoMat;

            if(nlhs>=6) {
                //Compute a minimum cost partition to return.
                size_t *minCostClustPartition=new size_t[N];

                if(minCostClustPartition==NULL) {
                    mexErrMsgTxt("Error allocating temporary buffer.");
                }

                clusterInfo2MinCostClustPartition(KCur,maxInClust,numPerCluster,clusterInfo,minCostClustPartition);

                //Copy into a double Matlab matrix to return.
                mxArray *minCostClustPartitionMat=mxCreateDoubleMatrix(N,1,mxREAL);
                double *minCostClustPartitionData=mxGetDoubles(minCostClustPartitionMat);

                for(size_t k=0;k<N;k++) {
                    //The +1 changes the C indices to Matlab indices.
                    minCostClustPartitionData[k]=static_cast<double>(minCostClustPartition[k]+1);
                }

                plhs[5]=minCostClustPartitionMat;
                delete[] minCostClustPartition;

                if(nlhs>=7) {
                    //Size to fit.
                    mxSetN(wAllMat,curAllIdx);
                    plhs[6]=wAllMat;
                    if(nlhs>=8) {
                        {//Size to fit.
                            const mwSize nDim=3;
                            const mwSize dims[]={xDim,maxComp,curAllIdx};
                            mxSetDimensions(muAllMat, dims, nDim);
                        }

                        plhs[7]=muAllMat;

                        if(nlhs>=9) {
                            {//Size to fit.
                                const mwSize nDim=4;
                                const mwSize dims[]={xDim,xDim,maxComp,curAllIdx};
                                mxSetDimensions(PAllMat, dims, nDim);
                            }
                            plhs[8]=PAllMat;

                            if(nlhs>=10) {
                                mxArray *clusterInfoAllMat;
                                {
                                    const mwSize nDim=3;
                                    const mwSize dims[]={N,maxComp,curAllIdx};
                                    clusterInfoAllMat=mxCreateNumericArray(nDim, dims,mxDOUBLE_CLASS, mxREAL);
                                }
                                double *clusterInfoAllData=mxGetDoubles(clusterInfoAllMat);

                                for(size_t k=0;k<N*maxComp*curAllIdx;k++) {
                                    if(clusterInfoAll[k]==SIZE_MAX) {
                                        clusterInfoAllData[k]=std::numeric_limits<double>::infinity();
                                    } else {//+1 for Matlab indexation.
                                        clusterInfoAllData[k]=static_cast<double>(clusterInfoAll[k]+1);
                                    }
                                }

                                plhs[9]=clusterInfoAllMat;

                                if(nlhs>=11) {
                                    mxArray *numPerClusterAllMat=mxCreateDoubleMatrix(N,curAllIdx,mxREAL);
                                    double *numPerClusterAllData=mxGetDoubles(numPerClusterAllMat);

                                    for(size_t k=0;k<N*curAllIdx;k++) {
                                        numPerClusterAllData[k]=static_cast<double>(numPerClusterAll[k]);
                                    }

                                    plhs[10]=numPerClusterAllMat;
                                }
                            } 
                        } else {
                            mxDestroyArray(PAllMat);
                        }

                        delete[] clusterInfoAll;
                        delete[] numPerClusterAll;
                    } else {
                        mxDestroyArray(muAllMat);
                        mxDestroyArray(PAllMat);
                    }
                }
            }
        }
    }

    delete[] clusterInfo;
    delete[] numPerCluster;
}

void RunnalsGaussMixRedCPP(const size_t xDim,const size_t N, double *w , double *mu, double *P,
                           const size_t K, const double gammaBound, const size_t KMax, size_t &KCur,
                           size_t *numPerCluster,size_t *clusterInfo,size_t &maxInClust,void *tempBuffer,
                          double *wAll,double *muAll,double *PAll,size_t *clusterInfoAll,size_t *numPerClusterAll, size_t &curAllIdx) { 
    //Considering the inputs, numPerCluster is length N and clusterInfo is NXN.
    //If wAll is not NULL, then wAll,muAll,PAll,clusterInfoAll,numPerClusterAll
    //are all filled in. If wAll is NULL, then they are not filled in. Only
    //wAll is checked for being NULL.

    const size_t N2=N*N;//The number of elements in each matrix.
    //The number of elements in each covariance matrix.
    const size_t xDim2=xDim*xDim;
    //Break up the temporary buffer into parts.
    uint8_t *curBuffPtr=reinterpret_cast<uint8_t*>(tempBuffer);
    double *BDistBuffer=reinterpret_cast<double*>(curBuffPtr);
    curBuffPtr+=BDistBufferSize(xDim);
    double *diff1=reinterpret_cast<double*>(curBuffPtr);
    curBuffPtr+=sizeof(double)*xDim;
    double *diff2=reinterpret_cast<double*>(curBuffPtr);
    curBuffPtr+=sizeof(double)*xDim;
    double *M=reinterpret_cast<double*>(curBuffPtr);
    curBuffPtr+=sizeof(double)*N2;
    double *wLogDetPVals=reinterpret_cast<double*>(curBuffPtr);
    curBuffPtr+=sizeof(double)*N;
    bool *selIdxPresent=reinterpret_cast<bool*>(curBuffPtr);//A length N buffer.

    //If no reduction is necessary.
    if(N<=K) {
        std::fill_n(numPerCluster,N,static_cast<size_t>(1));//numPerCluster=ones(N,1);

        maxInClust=1;
        for(size_t k=0;k<N;k++) {
            clusterInfo[k]=k;
        }
        KCur=N;

        if(wAll!=NULL) {
            std::copy(w,w+N,wAll);
            std::copy(mu,mu+xDim*N,muAll);
            std::copy(P,P+xDim2*N2,PAll);
            for(size_t k=0;k<N;k++) {
                clusterInfoAll[k*N]=k;
                numPerClusterAll[k]=1;
            }
            curAllIdx=1;
        }
        return;
    } else {
        //These variables are only used if all steps of the reduction
        //process are to be saved.
        size_t totalRedHyps=0;
        size_t maxComp=0;

        if(wAll!=NULL) {
            //K is the minimum number of components.
            getSizes4ReturningSteps(N, K, KMax, maxComp, totalRedHyps);
            curAllIdx=0;

            //If the inital mixture should be saved.
            if(N==maxComp) {
                std::copy(w,w+N,wAll);
                std::copy(mu,mu+xDim*N,muAll);
                std::copy(P,P+xDim2*N,PAll);
                for(size_t k=0;k<N;k++) {
                    clusterInfoAll[k*N]=k;
                    numPerClusterAll[k]=1;
                }
                curAllIdx++;
            }
        }

        //Set the elements of N to Inf. Only one triangle of this matrix
        //will actually be used after this. This is the cost matrix.
        std::fill_n(M,N2,std::numeric_limits<double>::infinity());
        const auto xDimL=static_cast<Eigen::EigenBase<Eigen::MatrixXd>::Index>(xDim);
        for(size_t k=0;k<N;k++) {
            const Eigen::Map<Eigen::MatrixXd> PCurEigen(P+k*xDim2,xDimL,xDimL);
            const double detVal=PCurEigen.determinant();
          
            wLogDetPVals[k]=w[k]*std::log(std::fabs(detVal));
        }

        //Fill in the cost matrix with the cost of all pairs.
        for(size_t cur1=0;cur1<N-1;cur1++) {
            for(size_t cur2=cur1+1;cur2<N;cur2++) {
                const size_t idx=cur1+cur2*N;

                M[idx]=BDist(xDim,w[cur1],w[cur2],mu+xDim*cur1,mu+xDim*cur2,P+xDim2*cur1,P+xDim2*cur2,wLogDetPVals[cur1],wLogDetPVals[cur2],BDistBuffer);
            }
        }

        //Setting the elements all to SIZE_MAX leaves a quick way to see if
        //an element was never changed.
        std::fill_n(clusterInfo,N2,SIZE_MAX);
        std::fill_n(numPerCluster,N,static_cast<size_t>(1));
        for(size_t k=0;k<N;k++) {
            //The first element in each columns of the NXN matrix.
            clusterInfo[k*N]=k;
        }
    
        KCur=N;
        std::fill_n(selIdxPresent,N,true);
        for(size_t mergeRound=0;mergeRound<N-K;mergeRound++) {
            //Find the index and value of the minimum element in M
            size_t minIdx=0;
            double minVal=M[0];
            for(size_t k=1;k<N2;k++) {
                if(M[k]<minVal) {
                    minVal=M[k];
                    minIdx=k;
                }
            }

            //Convert the minimum linear index into a row and column
            //number. Get the minimum row and columns from minIdx and
            //knowing that there are N rows and N columns.
            const size_t minCol=minIdx/N;//Integer division.
            const size_t minRow=minIdx%N;//Remainder
       
            if(minVal>gammaBound&&KCur<=KMax) {
                //If the distribution has been sufficiently reduced in
                //terms of cost.
                break;
            }
       
            //Now we know which two hypotheses to merge: The ones with
            //indices minRow and minCol. We will merge those hypotheses
            //and put the results in minRow.
            const double *muRow=mu+xDim*minRow;
            const double *PMinRow=P+xDim2*minRow;
            const double *muCol=mu+xDim*minCol;
            const double *PMinCol=P+xDim2*minCol;

            const double wSum=w[minRow]+w[minCol];
            const double w1=w[minRow]/wSum;
            const double w2=w[minCol]/wSum;
            w[minRow]=wSum;
            //Where to start saving mu and P. We are overwriting the
            //original values.
            double *muMerged=mu+xDim*minRow;
            double *PMerged=P+xDim2*minRow;
            
            //The below lines are doing:
            //muMerged=w1*mu(:,minRow)+w2*mu(:,minCol);
            //diff1=mu(:,minRow)-muMerged;
            //diff2=mu(:,minCol)-muMerged;
            //However, the ordering is important, because we are
            //overwriting the original mu values with muMerged, so we have
            //to take the differences before overwriting it.
            for(size_t k=0;k<xDim;k++) {
                const double muMergedCur=w1*muRow[k]+w2*muCol[k];
                diff1[k]=muRow[k]-muMergedCur;
                diff2[k]=muCol[k]-muMergedCur;

                muMerged[k]=muMergedCur;
            }
          
            //PMerged(:,:,minRow)=w1*(P(:,:,minRow)+diff1*diff1')+w2*(P(:,:,minCol)+diff2*diff2');
            {
                size_t idx=0;//Linear indexation that will follow the row-column indexation.
                for(size_t col=0;col<xDim;col++) {
                    for(size_t row=0;row<xDim;row++) {
                        PMerged[idx]=w1*(PMinRow[idx]+diff1[col]*diff1[row])+w2*(PMinCol[idx]+diff2[col]*diff2[row]);
                        idx++;
                    }
                }
            }

            //wLogDetPVals(minRow)=wSum*log(det(PMerged));   
            const Eigen::Map<Eigen::MatrixXd> PMergedEigen(PMerged,xDimL,xDimL);
            const double detVal=PMergedEigen.determinant();
            wLogDetPVals[minRow]=wSum*std::log(detVal);

            //Update the cluster information.
            const size_t numInRow=numPerCluster[minRow];
            const size_t numInCol=numPerCluster[minCol];
            const size_t numTotal=numInRow+numInCol;
            //This is the starting index of
            //clusterInfo((numInRow+1):numTotal,minRow) if one were
            //using Matlab.
            size_t *clusterIdx1=clusterInfo+N*minRow+numInRow;
            //This is the starting index of
            //clusterInfo(1:numInCol,minCol); if one were using Matlab.
            const size_t *clusterIdx2=clusterInfo+N*minCol;
            //clusterInfo((numInRow+1):numTotal,minRow)=clusterInfo(1:numInCol,minCol);
            for(size_t k=0;k<numInCol;k++) {
                clusterIdx1[k]=clusterIdx2[k];
            }

            numPerCluster[minCol]=0;
            numPerCluster[minRow]=numTotal;

            //The column has been removed.
            selIdxPresent[minCol]=false;

            //Make all the costs infinite so that nothing will be
            //assigned to the removed column.
            //M(minCol,:)=Inf;   
            {
                size_t idx=minCol;
                for(size_t k=0;k<N;k++) {
                    M[idx]=std::numeric_limits<double>::infinity();
                    idx=idx+N;
                }
            }

            //M(:,minCol)=Inf;
            {
                size_t idx=N*minCol;
                for(size_t k=0;k<N;k++) {
                    M[idx]=std::numeric_limits<double>::infinity();
                    idx++;
                }
            }
        
            //We must now fill in the costs for the merged estimate, which
            //is in minRow. We shall fill the cost matrix with the cost of
            //all pairs.
            if(minRow>0) {
                for(size_t cur1=0;cur1<minRow;cur1++) {
                    if(selIdxPresent[cur1]) {
                        const size_t idx=cur1+minRow*N;
                        M[idx]=BDist(xDim,w[cur1],w[minRow],mu+xDim*cur1,mu+xDim*minRow,P+xDim2*cur1,P+xDim2*minRow,wLogDetPVals[cur1],wLogDetPVals[minRow],BDistBuffer);
                    }
                }
            }
        
            for(size_t cur2=minRow+1;cur2<N;cur2++) {
                if(selIdxPresent[cur2]) {
                    const size_t idx=minRow+cur2*N;
                    M[idx]=BDist(xDim,w[minRow],w[cur2],mu+xDim*cur2,mu+xDim*minRow,P+xDim2*minRow,P+xDim2*cur2,wLogDetPVals[minRow],wLogDetPVals[cur2],BDistBuffer);
                }
            }

            KCur--;
            //KCur is now the number of components left (The number of ones
            //in selIdxPresent).

            if(wAll!=NULL&&KCur<=maxComp) {
                //Copy the data on the remaining KCur coefficients.
                double *wStart=wAll+maxComp*curAllIdx;
                double *muAllStart=muAll+xDim*maxComp*curAllIdx;
                double *PAllStart=PAll+xDim2*maxComp*curAllIdx;
                size_t *clusterInfoAllStart=clusterInfoAll+N*maxComp*curAllIdx;
                size_t *numPerClusterAllStart=numPerClusterAll+maxComp*curAllIdx;

                size_t curIdx=0;
                for(size_t k=0;k<N;k++) {
                    if(selIdxPresent[k]) {
                        wStart[curIdx]=w[k];
                        std::copy(mu+k*xDim,mu+(k+1)*xDim,muAllStart+curIdx*xDim);
                        std::copy(P+k*xDim2,P+(k+1)*xDim2,PAllStart+curIdx*xDim2);
                        std::copy(clusterInfo+k*N,clusterInfo+(k+1)*N,clusterInfoAllStart+curIdx*N);
                        numPerClusterAllStart[curIdx]=numPerCluster[k];

                        curIdx++;
                    }
                }

                curAllIdx++;
            }
        }

        //The merged items are marked with selIdxPresent and can now be all
        //grouped together to be returned. This is doing
        // w=w(selIdxPresent);
        // w=w/sum(w);%Guarantee continued normalization.
        // mu=mu(:,selIdxPresent);
        // P=P(:,:,selIdxPresent);
        // numPerCluster=numPerCluster(selIdxPresent);
        // maxInClust=max(numPerCluster);
        // clusterInfo=clusterInfo(1:maxInClust,selIdxPresent);
        size_t idx2Copy2=0;
        double wSum=0;
        maxInClust=0;
        for(size_t curIdx=0;curIdx<N;curIdx++) {
            if(selIdxPresent[curIdx]) {
                w[idx2Copy2]=w[curIdx];
                wSum+=w[curIdx];
                numPerCluster[idx2Copy2]=numPerCluster[curIdx];
                maxInClust=std::max(maxInClust,numPerCluster[curIdx]);

                size_t startIdx2Copy=xDim*idx2Copy2;
                size_t startIdx=xDim*curIdx;
                for(size_t k=0;k<xDim;k++) {
                    mu[startIdx2Copy+k]=mu[startIdx+k];
                }

                startIdx2Copy=xDim2*idx2Copy2;
                startIdx=xDim2*curIdx;
                for(size_t k=0;k<xDim2;k++) {
                    P[startIdx2Copy+k]=P[startIdx+k];
                }

                idx2Copy2++;
            }
        }
    
        //Normalize the w values.
        for(size_t k=0;k<KCur;k++) {
            w[k]/=wSum;
        }

        //Shrink the clusterInfo matrix so that the information is packed
        //as tightly as possible.
        idx2Copy2=0;
        for(size_t curIdx=0;curIdx<N;curIdx++) {
            if(selIdxPresent[curIdx]){
                const size_t startIdx=N*curIdx;
                const size_t startIdx2Copy=maxInClust*idx2Copy2;

                for(size_t k=0;k<maxInClust;k++) {
                    clusterInfo[startIdx2Copy+k]=clusterInfo[startIdx+k];
                }
                idx2Copy2++;
            }
        }
    }
}

double BDist(const size_t xDim,const double &w1,const double &w2,const double * const mu1,const double * const mu2,const double * const P1,const double * const P2,const double wLogDetP1,const double wLogDetP2, double *buffer) {
    const double wSum=w1+w2;
    const double w1m=w1/wSum;
    const double w2m=w2/wSum;
    //Break the memory in the buffer up into parts.
    double *P12=buffer;
    double *diff=buffer+xDim*xDim;

    for(size_t k=0;k<xDim;k++) {
        diff[k]=mu1[k]-mu2[k];
    }

    for(size_t i2=0;i2<xDim;i2++) {//Going by column.
        const size_t colOffset=i2*xDim;
        
        for(size_t i1=i2;i1<xDim;i1++) {
            const size_t curIdx=i1+colOffset;
            const size_t flipIdx=i2+i1*xDim;

            P12[curIdx]=w1m*P1[curIdx]+w2m*P2[curIdx]+w1m*w2m*diff[i1]*diff[i2];
            //The matrix is symmetric.
            P12[flipIdx]=P12[curIdx];
        }
    }

    const auto xDimL=static_cast<Eigen::EigenBase<Eigen::MatrixXd>::Index>(xDim);
    const Eigen::Map<Eigen::MatrixXd> P12Eigen(P12,xDimL,xDimL);
    const double detP12=P12Eigen.determinant();
    double val=wSum*std::log(detP12)-wLogDetP1-wLogDetP2;

    //Deal with the case where w1 and w2 are both essentially zero.
    if(std::isfinite(val)==false) {
        val=0;
    }

    return val;
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
*OF RECIPIENT IN THE USE OF THE SOFTWARE.*/
