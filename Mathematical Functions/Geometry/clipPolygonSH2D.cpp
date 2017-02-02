/*CLIPPOLYGONSH2D A C++ implementation of a function to clip a
 *                non-self-intersecting polygon to a region specified by a
 *                convex polygon using a 2D implementation of the
 *                Sutherland-Hodgman algorithm. If the polygon being
 *                clipped was concave at the vertices outside of the
 *                clipping region, the resulting clipped polygon will have
 *                overlapping edges on the boundary of the region rather
 *                than being slip into multiple separate polygons.
 *
 *INPUTS: polygon2Clip A 2XN list of N>=3 vertices making up the polygon of
 *                    the form [x;y]. It does not matter whether vertices
 *                    are repeated.
 *  convexClipPolygon A 2XNClip list of NClip>=3 vertices making up the
 *                    convex clipping polygon. Vertices should not be
 *                    repeated. The first vertex should not be repeated on
 *                    the end. The vertices should be in a counterclockwise
 *                    order.
 *
 *OUTPUTS: clippedPolygon The polygon polygon2Clip clipped to the convex
 *                       region specified by convexClipPolygon. If
 *                       polygon2Clip is completely does not intersect with
 *                       the clipping region, an empty matrix is returned.
 *                       If the clipping region is engulfed by the polygon,
 *                       then the polygon will be the clipping region.
 *
 *The algorithm is described in more detail in the comments to the Matlab
 *implementation.
 *
 *The Sutherland-Hodgman algorithm requires that the vertices in the
 *clipping polygon be in counterclockwise order.
 *
 *The maximum number of vertices on the clipped polygon is if every other
 *vertex of the polygon to clip is outside of the clipping region, as each
 *exit/ entry to the clipping region adds two vertices. Thus, the maximum
 *number of edges is 3/2 times the original number.
 *
 * The algorithm can be compiled for use in Matlab  using the 
 * CompileCLibraries function.
 *
 * The algorithm is run in Matlab using the command format
 * clippedPolygon=clipPolygonSH2D(polygon2Clip,convexClipPolygon);
 *
 *December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

//This header is required by Matlab.
#include "mex.h"
//This is for input validation
#include "MexValidation.h"
//To determine the sign of a sum with no finite precision errors.
#include "exactSignOfSumCPP.hpp"
//To determine the intersection point of two lines.
#include "mathGeometricFuncs.hpp"

//Prototype
bool vertexIsInsideClipEdge(const double *vertex,const double *clipVertex1,const double *clipEdge);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    size_t numPolyVertices, numClipVertices; 
    double *convexClipPolygon, *curPolygon, *clippedPolygonNew;
    double *prevClipVertex;
    mxArray *curPolygonMat, *clippedPolygonNewMat;
    size_t curClip,numVertices;
            
    if(nrhs!=2) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Incorrect number of outputs.");
        return;
    }
    
    //Return an empty matrix if an empty polygon is passed.
    if(mxIsEmpty(prhs[0])){
        plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
        return;
    }
    
    if(mxGetM(prhs[0])!=2||mxGetM(prhs[1])!=2) {
        mexErrMsgTxt("The polygons must be two-dimensional.");
        return;
    }
    
    checkRealDoubleArray(prhs[0]);
    numPolyVertices=mxGetN(prhs[0]);
    if(numPolyVertices<3) {
        mexErrMsgTxt("The polygon to be clipped must have at least three vertices.");
        return;
    }

    checkRealDoubleArray(prhs[1]);
    numClipVertices=mxGetN(prhs[1]);
    if(numClipVertices<3) {
        mexErrMsgTxt("The clipping polygon must have at least three vertices.");
        return;
    }
    convexClipPolygon=reinterpret_cast<double*>(mxGetData(prhs[1]));
    
    //The clipping polygon must have its vertices going in a
    //counterclockwise order. Check whether that is the case. If not, then
    //reverse the order.
    if(signedPolygonAreaCPP(convexClipPolygon,numClipVertices)<0) {
        mexErrMsgTxt("The vertices of the clipping polygon should be in a counterclockwise order. Reverse the order of the vertices and try again.");
        return;
    }
    
    //Allocate space for the new polygon and for temporary scratch space.
    {
        //The maximum possible number of vertices in the clipped polygon.
        const size_t maxClippedVertices=(numPolyVertices*3)/2+1;
        size_t i;
        double *inputEls=reinterpret_cast<double*>(mxGetData(prhs[0]));
                
        curPolygonMat=mxCreateDoubleMatrix(2,maxClippedVertices,mxREAL);
        clippedPolygonNewMat=mxCreateDoubleMatrix(2,maxClippedVertices,mxREAL);
        
        curPolygon=reinterpret_cast<double*>(mxGetData(curPolygonMat));
        clippedPolygonNew=reinterpret_cast<double*>(mxGetData(clippedPolygonNewMat));
        
        //Copy the values in prhs[0] into curPolygon.        
        for(i=0;i<2*numPolyVertices;i++) {
            curPolygon[i]=inputEls[i];
        }
        
        //The current number of vertices in curPolygon.
        numVertices=numPolyVertices;
    }
    
    //The first clipping edge will be the one from the end to the
    //beginning.
    prevClipVertex=convexClipPolygon+2*(numClipVertices-1);
    //For each edge, create the reduced polygon by clipping with that edge. 
    for(curClip=0;curClip<numClipVertices;curClip++) {
        double *curClipVertex=convexClipPolygon+2*curClip;
        double *prevVertex;
        //The polygon will be clipped to this edge this iteration.
        const double curClipEdge[2]={curClipVertex[0]-prevClipVertex[0],curClipVertex[1]-prevClipVertex[1]};
        size_t curV, numVerticesNew;
        double *curClippedPolygonNew=clippedPolygonNew;
        
        //Number of vertices that have been added to clippedPolygonNew
        numVerticesNew=0;
        prevVertex=curPolygon+2*(numVertices-1);
        
        for(curV=0;curV<numVertices;curV++) {
            double *curVertex=curPolygon+2*curV;
            
            //Note that the clip edges are infinitely extended so that
            //vertices outside of the clipping region are involved.
            if(vertexIsInsideClipEdge(curVertex,prevClipVertex,curClipEdge)) {
                //If the current vertex is inside of the clipping region,
                //then add it, but if the previous vertex was not in the
                //clipping region, then an extra vertex at the edge of the
                //boundary region needs to be added.
                if(!vertexIsInsideClipEdge(prevVertex,prevClipVertex,curClipEdge)) {
                    //[prevVertex,curVertex]
                    const double line1[4]={prevVertex[0],prevVertex[1],curVertex[0],curVertex[1]};
                    //[curClipVertex,prevClipVertex]
                    const double line2[4]={curClipVertex[0],curClipVertex[1],prevClipVertex[0],prevClipVertex[1]};
                    
                    twoLineIntersectionPoint2DCPP(line1,line2,curClippedPolygonNew);
                    curClippedPolygonNew+=2;
                    numVerticesNew++;
                }

                curClippedPolygonNew[0]=curVertex[0];
                curClippedPolygonNew[1]=curVertex[1];
                curClippedPolygonNew+=2;
                numVerticesNew++;
            }else if(vertexIsInsideClipEdge(prevVertex,prevClipVertex,curClipEdge)) {
                //If the previous vertex was inside of the clipping region
                //and this vertex is not, then add a line segment from the
                //previous vertex to the edge of the clipping region. 
                //[prevVertex,curVertex]
                const double line1[4]={prevVertex[0],prevVertex[1],curVertex[0],curVertex[1]};
                //[curClipVertex,prevClipVertex]
                const double line2[4]={curClipVertex[0],curClipVertex[1],prevClipVertex[0],prevClipVertex[1]};

                twoLineIntersectionPoint2DCPP(line1,line2,curClippedPolygonNew);
                curClippedPolygonNew+=2;
                numVerticesNew++;
            }
            
            prevVertex=curVertex;
        }

        {//Swap the buffers.
            double *temp=curPolygon;
            curPolygon=clippedPolygonNew;
            numVertices=numVerticesNew;
            clippedPolygonNew=temp;
        }
        
        //The object is not in the viewing area at all.
        if(numVertices==0) {
            mxFree(curPolygonMat);
            mxFree(clippedPolygonNewMat);
            plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
            return;
        }
        
        prevClipVertex=curClipVertex;
    }

    //curPolygon is what should be returned. This corresponds to either 
    //curPolygonMat or clippedPolygonNewMat, since the buffers were
    //exchanged during the loops. Determine which is curPolygon and return
    //it, freeing the other. 
    if(curPolygon==reinterpret_cast<double*>(mxGetData(curPolygonMat))) {
        //Resize to fit.
        mxSetN(curPolygonMat,numVertices);
        plhs[0]=curPolygonMat;
        mxDestroyArray(clippedPolygonNewMat);
    } else  {
        mxSetN(clippedPolygonNewMat,numVertices);
        plhs[0]=clippedPolygonNewMat;
        mxDestroyArray(curPolygonMat);
    }
}

bool vertexIsInsideClipEdge(const double *vertex,const double *clipVertex1,const double *clipEdge) {
//Determine whether a vertex is inside of the clipping region given an
//edge. This relies on the fact that the clipping region is convex and goes
//in a counterclockwise direction. Thus, the sign of the angle with respect
//to the clipping edge is used to determine whether it is inside or outside.
//This could be done by extending the vectors to 3D (padding a zero on the
//ends) and then using a cross product, whereby the z component would tell
//us the sign, but the value of the z-component is the same as the
//determinant value used here.
//
//INPUTS:  vertex The 2X1 vertex to test whether it is on the correct size
//                of the edge to be in the clipping region.
//    clipVertex1 The first 2X1 vertex of the boundary used to form the
//                edge for clipping.
//       clipEdge The edge vector, which would be clipVertex2-clipVertex1.
//
//The code below is equivalent to
//isInsideClipEdge=det([vertex-clipVertex1,clipEdge])<=0;
//but should be slightly less susceptible to finite precision errors.
    double S[4];
    
    S[0]=-clipEdge[0]*vertex[1];
    S[1]=clipEdge[1]*vertex[0];
    S[2]=-clipEdge[1]*clipVertex1[0];
    S[3]=clipEdge[0]*clipVertex1[1];
    
    return exactSignOfSumCPP(S,4)<=0;
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
