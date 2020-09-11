function [aVals,aJacob,aHess,papt]=aSpiralSimp(xPoints,numRetDims)
%%ASPIRALSIMP The drift function for a non-ballistic spiraling target
%         motion model in 3 dimensions, formulated in a manner such that it
%         is simple to make a spiraling target follow a nominal trajectory.
%         This model is probably bad to design into a target tracking
%         algorithm, but might be good for designing spiraling trajectories
%         to use to test target tracking algorithms.
%
%INPUTS: xPoints The target state at time t. It consists of position 
%                (3 elements), instantaneous velocity (3 elements), the
%                velocity vector (3 elements) of the overall direction of
%                motion of the spiraling model (ground speed) and the
%                spiral rate of the target. Thus xDim=10. If x is an
%                xDim X numStates matrix, then the spiraling model is
%                evaluated for all of the state vectors. The last 4
%                elements in xPoints are taken to be constants and thus
%                aVal entries will be 0 for them. If only aVals is
%                requested on the output, then xPoints can be a matrix of
%                points.
%     numRetDims If numRetDims=6, then the last 4 elements in xPoints (the
%                velocity vector of the overall spiraling model and the
%                spiral rate) are taken to be constants and the returned
%                aVals Vec is 6X1. Otherwise, numRetDims can be 10 and
%                aVals is 10X1. The default if omitted or an empty matrix
%                is passed is 6.
%
%OUTPUTS: aVals The 6X1 (or 10X1 depending on numRetDims) flat-Earth time-
%               derivative of the state. If xPoints was a matrix of N
%               points, then this will be 6XN (or 10XN).
%        aJacob The 6X6 (or 10X10) matrix of partial derivatives of aVals
%               such that aJacob(:,k) is the partial derivative of
%               aVals(:,k) with respect to
%               xPoints(k). xPoints can't be a matrix if this output is
%               requested.
%         aHess The 6X6X6 (or 10X10X10) matrix of second derivatives of
%               aVals such that aHess(:,k1,k2) is the second partial
%               derivative of aVals with respect to xPoints(k1) and
%               xPoints(k2).
%          papt The 6X1 or 10X1 derivative with resect to time of aVals.
%               This is all zeros, because the model is time invariant.
%
%A derivation of the simplified flat-Earth spiraling dynamic model is given
%in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Simulating aerial targets in 3D accounting for the
%    Earth's curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, Jun. 2015.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<3||isempty(numRetDims))
        numRetDims=6;
    end

    if(numRetDims<6||numRetDims>10)
       error('numRetDims must be between 6 and 10.')
    end

    numPoints=size(xPoints,2);
    aVals=zeros(numRetDims,numPoints);
    
    vl=xPoints(7:9,:);
    omega=xPoints(10,:);
    
    vlMag2=sum(vl.*vl,1);
    vlMag=sqrt(vlMag2);
    
    aVals(1:3,:)=xPoints(4:6,:);
    vs=xPoints(4:6,:)-vl;
    ul=bsxfun(@rdivide,vl,vlMag);
    vsDot=omega.*bsxfun(@cross,ul,vs);
    aVals(4:6,:)=vsDot;
    
    if(nargout>1)
        rdotx=xPoints(4);
        rdoty=xPoints(5);
        rdotz=xPoints(6);

        vlx=vl(1);
        vly=vl(2);
        vlz=vl(3);
        
        vlx2=vlx*vlx;
        vly2=vly*vly;
        vlz2=vlz*vlz;

        vlMag2=vlMag*vlMag;
        vlMag3=vlMag2*vlMag;
        
        pXrdy=-(omega*vlz)/vlMag;
        pXrdz=(omega*vly)/vlMag;

        pYrdx=(omega*vlz)/vlMag;
        pYrdz=-((omega*vlx)/vlMag);

        pZrdx=-((omega*vly)/vlMag);
        pZrdy=(omega*vlx)/vlMag;

        if(numRetDims==6)
            aJacob=[0,0,0,1,        0,      0;
                    0,0,0,0,        1,      0;
                    0,0,0,0,        0,      1;
                    0,0,0,0,        pXrdy,  pXrdz;
                    0,0,0,pYrdx,    0,      pYrdz;
                    0,0,0,pZrdx,    pZrdy,  0];
        else
            pXvx=(omega*vlx*(-rdotz*vly+rdoty*vlz))/vlMag3;
            pXvy=(omega*(rdoty*vly*vlz+rdotz*(vlx2+vlz2)))/vlMag3;
            pXvz=-((omega*(rdoty*(vlx2+vly2)+rdotz*vly*vlz))/vlMag3);

            pXomega=(rdotz*vly-rdoty*vlz)/vlMag;

            pYvx=-((omega*(rdotx*vlx*vlz+rdotz*(vly2+vlz2)))/vlMag3);
            pYvy=(omega*vly*(rdotz*vlx-rdotx*vlz))/vlMag3;
            pYvz=(omega*(rdotx*(vlx2+vly2)+rdotz*vlx*vlz))/vlMag3;
            pYomega=(-rdotz*vlx+rdotx*vlz)/vlMag;
            
            pZvx=(omega*(rdotx*vlx*vly+rdoty*(vly2+vlz2)))/vlMag3;
            pZvy=-((omega*(rdoty*vlx*vly+rdotx*(vlx2+vlz2)))/vlMag3);
            pZvz=(omega*(-rdoty*vlx+rdotx*vly)*vlz)/vlMag3;
            pZomega=(rdoty*vlx-rdotx*vly)/vlMag;

            aJacob=[0,0,0,1,        0,      0,      0,      0,      0,      0;
                    0,0,0,0,        1,      0,      0,      0,      0,      0;
                    0,0,0,0,        0,      1,      0,      0,      0,      0;
                    0,0,0,0,        pXrdy,  pXrdz,  pXvx,   pXvy,   pXvz,   pXomega;
                    0,0,0,pYrdx,    0,      pYrdz,  pYvx,   pYvy,   pYvz,   pYomega;
                    0,0,0,pZrdx,    pZrdy,  0,      pZvx,   pZvy,   pZvz,   pZomega;
                    zeros(4,10)];
        end

        if(nargout>2)
            aHess=zeros(numRetDims,numRetDims,numRetDims);
 
            vlMag5=vlMag3*vlMag2;
            
            %Second deriavtives are all 0 if only 
            if(numRetDims>6)
                pXrdyvx=(omega*vlx*vlz)/vlMag3;
                pXrdyvy=(omega*vly*vlz)/vlMag3;
                pXrdyvz=-((omega*(vlx2+vly2))/vlMag3);
                pXrdyOmega=-(vlz/vlMag);

                pXrdzvx=-((omega*vlx*vly)/vlMag3);
                pXrdzvy=(omega*(vlx2+vlz2))/vlMag3;
                pXrdzvz=-((omega*vly*vlz)/vlMag3);
                pXrdzOmega=vly/vlMag;
                
                pXvxrdx=0;
                pXvxrdy=pXrdyvx;
                pXvxrdz=pXrdzvx;
                pXvxvx=-((omega*(rdotz*vly-rdoty*vlz)*(-2*vlx2+vly2+vlz2))/vlMag5);
                pXvxvy=-((omega*vlx*(3*rdoty*vly*vlz+rdotz*(vlx2-2*vly2+vlz2)))/vlMag5);
                pXvxvz=(omega*vlx*(3*rdotz*vly*vlz+rdoty*(vlx2+vly2-2*vlz2)))/vlMag5;
                pXvxOmega=(-rdotz*vlx*vly+rdoty*vlx*vlz)/vlMag3;

                pXvyrdx=0;
                pXvyrdy=pXrdyvy;
                pXvyrdz=pXrdzvy;
                pXvyvx=pXvxvy;
                pXvyvy=(omega*(-3*rdotz*vly*(vlx2+vlz2)+rdoty*vlz*(vlx2-2*vly2+vlz2)))/vlMag5;
                pXvyvz=(omega*(rdoty*vly*(vlx2+vly2-2*vlz2)-rdotz*vlz*(vlx2-2*vly2+vlz2)))/vlMag5;
                pXvyOmega=(rdoty*vly*vlz+rdotz*(vlx2+vlz2))/vlMag3;

                pXvzrdx=0;
                pXvzrdy=pXrdyvz;
                pXvzrdz=pXrdzvz;
                pXvzvx=pXvxvz;
                pXvzvy=pXvyvz;
                pXvzvz=(omega*(3*rdoty*(vlx2+vly2)*vlz-rdotz*vly*(vlx2+vly2-2*vlz2)))/vlMag5;
                pXvzOmega=(-rdoty*(vlx2+vly2)-rdotz*vly*vlz)/vlMag3;

                pXOmegardx=0;
                pXOmegardy=pXrdyOmega;
                pXOmegardz=pXrdzOmega;
                pXOmegavx=pXvxOmega;
                pXOmegavy=pXvyOmega;
                pXOmegavz=pXvzOmega;
                pXOmegaOmega=0;

                pYrdxvx=-((omega*vlx*vlz)/vlMag3);
                pYrdxvy=-((omega*vly*vlz)/vlMag3);
                pYrdxvz=(omega*(vlx2+vly2))/vlMag3;
                pYrdxOmega=vlz/vlMag;

                pYrdzvx=-((omega*(vly2+vlz2))/vlMag3);
                pYrdzvy=(omega*vlx*vly)/vlMag3;
                pYrdzvz=(omega*vlx*vlz)/vlMag3;
                pYrdzOmega=-(vlx/vlMag);
                                
                pYvxrdx=pYrdxvx;
                pYvxrdy=0;
                pYvxrdz=pYrdzvx;
                pYvxvx=(omega*(3*rdotz*vlx*(vly2+vlz2)-rdotx*vlz*(-2*vlx2+vly2+vlz2)))/vlMag5;
                pYvxvy=(omega*vly*(3*rdotx*vlx*vlz+rdotz*(-2*vlx2+vly2+vlz2)))/vlMag5;
                pYvxvz=(omega*(-rdotx*vlx*(vlx2+vly2-2*vlz2)+rdotz*vlz*(-2*vlx2+vly2+vlz2)))/vlMag5;
                pYvxOmega=(-rdotx*vlx*vlz-rdotz*(vly2+vlz2))/vlMag3;

                pYvyrdx=pYrdxvy;
                pYvyrdy=0;
                pYvyrdz=pYrdzvy;
                pYvyvx=pYvxvy;
                pYvyvy=(omega*(rdotz*vlx-rdotx*vlz)*(vlx2-2*vly2+vlz2))/vlMag5;
                pYvyvz=-((omega*vly*(3*rdotz*vlx*vlz+rdotx*(vlx2+vly2-2*vlz2)))/vlMag5);
                pYvyOmega=(vly*(rdotz*vlx-rdotx*vlz))/vlMag3;

                pYvzrdx=pYrdxvz;
                pYvzrdy=0;
                pYvzrdz=pYrdzvz;
                pYvzvx=pYvxvz;
                pYvzvy=pYvyvz;
                pYvzvz=(omega*(-3*rdotx*(vlx2+vly2)*vlz+rdotz*vlx*(vlx2+vly2-2*vlz2)))/vlMag5;
                pYvzOmega=(rdotx*(vlx2+vly2)+rdotz*vlx*vlz)/vlMag3;
            
                pYOmegardx=pYrdxOmega;
                pYOmegardy=0;
                pYOmegardz=pYrdzOmega;
                pYOmegavx=pYvxOmega;
                pYOmegavy=pYvyOmega;
                pYOmegavz=pYvzOmega;
                pYOmegaOmega=0;
                
                pZrdxvx=(omega*vlx*vly)/vlMag3;
                pZrdxvy=-((omega*(vlx2+vlz2))/vlMag3);
                pZrdxvz=(omega*vly*vlz)/vlMag3;
                pZrdxOmega=-(vly/vlMag);
                                
                pZrdyvx=(omega*(vly2+vlz2))/vlMag3;
                pZrdyvy=-((omega*vlx*vly)/vlMag3);
                pZrdyvz=-((omega*vlx*vlz)/vlMag3);
                pZrdyOmega=vlx/vlMag;

                pZvxrdx=pZrdxvx;
                pZvxrdy=pZrdyvx;
                pZvxrdz=0;
                pZvxvx=(omega*(-3*rdoty*vlx*(vly2+vlz2)+rdotx*vly*(-2*vlx2+vly2+vlz2)))/vlMag5;
                pZvxvy=(omega*(rdotx*vlx*(vlx2-2*vly2+vlz2)-rdoty*vly*(-2*vlx2+vly2+vlz2)))/vlMag5;
                pZvxvz=-((omega*vlz*(3*rdotx*vlx*vly+rdoty*(-2*vlx2+vly2+vlz2)))/vlMag5);
                pZvxOmega=(rdotx*vlx*vly+rdoty*(vly2+vlz2))/vlMag3;

                pZvyrdx=pZrdxvy;
                pZvyrdy=pZrdyvy;
                pZvyrdz=0;
                pZvyvx=pZvxvy;
                pZvyvy=(omega*(3*rdotx*vly*(vlx2+vlz2)-rdoty*vlx*(vlx2-2*vly2+vlz2)))/vlMag5;
                pZvyvz=(omega*vlz*(3*rdoty*vlx*vly+rdotx*(vlx2-2*vly2+vlz2)))/vlMag5;
                pZvyOmega=(-rdoty*vlx*vly-rdotx*(vlx2+vlz2))/vlMag3;
                
                pZvzrdx=pZrdxvz;
                pZvzrdy=pZrdyvz;
                pZvzrdz=0;
                pZvzvx=pZvxvz;
                pZvzvy=pZvyvz;
                pZvzvz=-((omega*(rdoty*vlx-rdotx*vly)*(vlx2+vly2-2*vlz2))/vlMag5);
                pZvzOmega=(-rdoty*vlx*vlz+rdotx*vly*vlz)/vlMag3;

                pZOmegardx=pZrdxOmega;
                pZOmegardy=pZrdyOmega;
                pZOmegardz=0;
                pZOmegavx=pZvxOmega;
                pZOmegavy=pZvyOmega;
                pZOmegavz=pZvzOmega;
                pZOmegaOmega=0;

                aHess(:,:,4)=[zeros(3,10);
                              0,0,0,0,0,0,  pXvxrdx,   pXvyrdx,   pXvzrdx,   pXOmegardx;
                              0,0,0,0,0,0,  pYvxrdx,   pYvyrdx,   pYvzrdx,   pYOmegardx;
                              0,0,0,0,0,0,  pZvxrdx,   pZvyrdx,   pZvzrdx,   pZOmegardx;
                              zeros(4,10)];
                aHess(:,:,5)=[zeros(3,10);
                              0,0,0,0,0,0,  pXvxrdy,   pXvyrdy,   pXvzrdy,   pXOmegardy;
                              0,0,0,0,0,0,  pYvxrdy,   pYvyrdy,   pYvzrdy,   pYOmegardy;
                              0,0,0,0,0,0,  pZvxrdy,   pZvyrdy,   pZvzrdy,   pZOmegardy;
                              zeros(4,10)];
                aHess(:,:,6)=[zeros(3,10);
                              0,0,0,0,0,0,  pXvxrdz,   pXvyrdz,   pXvzrdz,   pXOmegardz;
                              0,0,0,0,0,0,  pYvxrdz,   pYvyrdz,   pYvzrdz,   pYOmegardz;
                              0,0,0,0,0,0,  pZvxrdz,   pZvyrdz,   pZvzrdz,   pZOmegardz;
                              zeros(4,10)];     
                aHess(:,:,7)=[zeros(3,10);
                              0,0,0,0,          pXrdyvx,  pXrdzvx,  pXvxvx,   pXvyvx,   pXvzvx,   pXOmegavx;
                              0,0,0,pYrdxvx,    0,        pYrdzvx,  pYvxvx,   pYvyvx,   pYvzvx,   pYOmegavx;
                              0,0,0,pZrdxvx,    pZrdyvx,  0,        pZvxvx,   pZvyvx,   pZvzvx,   pZOmegavx;
                              zeros(4,10)];
                aHess(:,:,8)=[zeros(3,10);
                              0,0,0,0,          pXrdyvy,  pXrdzvy,  pXvxvy,   pXvyvy,   pXvzvy,   pXOmegavy;
                              0,0,0,pYrdxvy,    0,        pYrdzvy,  pYvxvy,   pYvyvy,   pYvzvy,   pYOmegavy;
                              0,0,0,pZrdxvy,    pZrdyvy,  0,        pZvxvy,   pZvyvy,   pZvzvy,   pZOmegavy;
                              zeros(4,10)];
                aHess(:,:,9)=[zeros(3,10);
                              0,0,0,0,          pXrdyvz,  pXrdzvz,  pXvxvz,   pXvyvz,   pXvzvz,   pXOmegavz;
                              0,0,0,pYrdxvz,    0,        pYrdzvz,  pYvxvz,   pYvyvz,   pYvzvz,   pYOmegavz;
                              0,0,0,pZrdxvz,    pZrdyvz,  0,        pZvxvz,   pZvyvz,   pZvzvz,   pZOmegavz;
                              zeros(4,10)];
               aHess(:,:,10)=[zeros(3,10);
                              0,0,0,0,             pXrdyOmega,  pXrdzOmega,  pXvxOmega,   pXvyOmega,   pXvzOmega,   pXOmegaOmega;
                              0,0,0,pYrdxOmega,    0,           pYrdzOmega,  pYvxOmega,   pYvyOmega,   pYvzOmega,   pYOmegaOmega;
                              0,0,0,pZrdxOmega,    pZrdyOmega,  0,           pZvxOmega,   pZvyOmega,   pZvzOmega,   pZOmegaOmega;
                              zeros(4,10)];
            end
            if(nargout>3)
                papt=zeros(10,1);
            end
        end
    end
end

%LICENSE:
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
%OF RECIPIENT IN THE USE OF THE SOFTWARE.
