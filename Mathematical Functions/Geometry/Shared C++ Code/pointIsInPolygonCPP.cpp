/**POINTISINPOLYGONCPP  A direct C++ implementation of a function that
 *              given a polygon specified by a number of vertices,
 *              determines whether a point is in the polygon. See the
 *              Matlab implementation pointIsInPolygon.m and the C++ for
 *              Matlab function, which calls this one pointIsInPolygon.cpp
 *              for more details on the function.
 *
 *December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mathGeometricFuncs.hpp"

//Prototypes for subroutines used here.
bool pointsAreEqual(const double *P,const size_t idx1,const size_t idx2);
ptrdiff_t omegaInc4Edge(const double *P,const double *R, const size_t idxI, const size_t idxIPlus1, double *detVal);

bool pointIsInPolygonCPP(const double *P, const size_t numVertices, const double *R, const bool boundaryIsImportant,ptrdiff_t *omega) {
//The direct C++ implementation that solves the problem for a single point.
//P is the set of 2D vertices and R is the 2D point.
//omega is set by the function as an additional return value.
//There must be at least 3 vertices in 
    
//P=the vertices of the polygon
//R=the point
    *omega=0;
    if(boundaryIsImportant==false) {
        double unused;
    //This is Algorithm 6, which is faster, but which does not return
    //consistent results for points that are on the boundary.
        size_t i;
        for(i=0;i<numVertices-1;i++) {
            size_t idxI=i*2;
            size_t idxIPlus1=(i+1)*2;
            
            *omega+=omegaInc4Edge(P,R,idxI, idxIPlus1,&unused);
        }
        
        //Increment omega for the edge that goes from the end back to the
        //beginning.
        *omega+=omegaInc4Edge(P,R,2*(numVertices-1),0,&unused);
    } else {
        //This is algorithm 7 which is slower, but which correctly identifies
        //points on the boundary as being inside the polygon.
        if((P[0]==R[0])&&(P[1]==R[1])) {
            //If it is on the first vertex (and thus in the polygon)
            return true;
        } else {
            size_t i;
            double detVal;
            for(i=0;i<numVertices-1;i++) {
                size_t idxI=i*2;
                size_t idxIPlus1=(i+1)*2;
                
                //The first check for whether the point is on the edge.
                if(P[1+idxIPlus1]==R[1]) {
                    if(P[0+idxIPlus1]==R[0]) {
                        //If it is on a vertex
                        return true;
                    } else {
                        if((P[1+idxI]==R[1])&&((P[0+idxIPlus1]>R[0])==(P[0+idxI]<R[0]))) {
                            //If it is on an edge
                            return true;
                        }
                    }
                }

                *omega+=omegaInc4Edge(P,R,idxI, idxIPlus1,&detVal);
                //The second check for whether the point is on the edge.
                if(detVal==0&&!pointsAreEqual(P,idxI,idxIPlus1)) {
                    return true;//The point is on the edge.
                }
            }
            //Deal with the edge going from the last vertex to the first
            //vertex.
            //The first check for whether the point is on the edge.
            if(P[1+0]==R[1]) {
                if(P[0+0]==R[0]) {
                    //If it is on a vertex
                    return true;
                } else {
                    if((P[1+numVertices-1]==R[1])&&((P[0+0]>R[0])==(P[0+numVertices-1]<R[0]))) {
                        //If it is on an edge
                        return true;
                    }
                }
            }
            *omega+=omegaInc4Edge(P,R,2*(numVertices-1), 0,&detVal);
            //The second check for whether the point is on the edge.
            if(detVal==0&&!pointsAreEqual(P,numVertices-1,0)) {
                return true;//The point is on the edge.
            }
        }
    }
    return  *omega%2!=0;
}

bool pointsAreEqual(const double *P,const size_t idx1,const size_t idx2) {
//True is the two 2D points are equal, false otherwise.
    return (P[idx1]==P[idx2])&&P[idx1+1]==P[idx2+1];
}

ptrdiff_t omegaInc4Edge(const double *P,const double *R, const size_t idxI, const size_t idxIPlus1,double *detVal) {
//Get the increment to omega for a single edge pair when the boundary does
//not matter. By breaking it out into a subroutine, we can handle the last
//edge going between the first and last vertices by just calling the
//function with different values of vertexIIdx and vertexIPlus1Idx.
    ptrdiff_t omegaInc=0;
    *detVal=1;//Not zero, for the case where this is not set below.

    if((P[1+idxI]<R[1])!=(P[1+idxIPlus1]<R[1])) {
        //If crossing
        if(P[0+idxI]>=R[0]) {
            if(P[0+idxIPlus1]>R[0]) {
                //Modify_omega
                omegaInc=2*(P[1+idxIPlus1]>P[1+idxI])-1;
            } else {
                *detVal=(P[0+idxI]-R[0])*(P[1+idxIPlus1]-R[1])-(P[0+idxIPlus1]-R[0])*(P[1+idxI]-R[1]);
                if((*detVal>0)==(P[1+idxIPlus1]>P[1+idxI])) {
                //If right_crossing
                    //Modify_omega
                    omegaInc=2*(P[1+idxIPlus1]>P[1+idxI])-1;
                }
            }
        } else {
            if(P[0+idxIPlus1]>R[0]) {
                *detVal=(P[0+idxI]-R[0])*(P[1+idxIPlus1]-R[1])-(P[0+idxIPlus1]-R[0])*(P[1+idxI]-R[1]);
                if((*detVal>0)==(P[1+idxIPlus1]>P[1+idxI])) {
                //If right_crossing
                    //Modify_omega
                    omegaInc=2*(P[1+idxIPlus1]>P[1+idxI])-1;
                }
            }
        }
    }

    return omegaInc;
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
