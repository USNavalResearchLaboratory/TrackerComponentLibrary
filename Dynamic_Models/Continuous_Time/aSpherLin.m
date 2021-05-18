function [aVal,aJacob,aHess,papt]=aSpherLin(x)
%%ASPHERLIN The drift function for a 3D continous-time constant direction
%           and speed dynamic model where the target state is given as
%           position, and a velocity that is given as azimuth (phi),
%           elevation (theta) and speed (s). Thus, the state is
%           [x;y;z;phi;theta;s].
%
%INPUTS: x The 6X1 set of target state vectors in 3D space in the order of
%          of [x;y;z;phi;theta;s] where the angles are in radians and the
%          coordinate system for azimuth and elevation corresponds to
%          system type 0 in the spher2Cart: Azimuth is measured 
%          counterclockwise from the x-axis in the x-y plane. Elevation
%          is measured up from the x-y plane (towards the z-axis).
%
%OUTPUTS: aVal The 6X1 time-derivative of the state under the constant
%              heading motion model.
%       aJacob This is the 6X6 matrix of partial derivatives of aVal such
%              that aJacob(:,i) is the partial derivative of aVal with
%              respect to x(i).
%        aHess The 6X6X6 matrix of second derivatives of aVal such that
%              aHess(:,k1,k2) is the second partial derivative of aVal with
%              respect to x(k1) and x(k2).
%         papt The 6X1 partial derivative with resect to time of aVal. This
%              is all zeros, because the model is time invariant.
%
%The model here is model 1 in [1]. The drift function corresponding to that
%used in [1] just adds the noise component to the last 3 terms and thus one
%could use DPoly(1,q0,1,3).
%
%REFERENCES:
%[1] D. Laneuville, "New models for 3D maneuvering target tracking," in New
%    Models for 3D Maneuvering Target Tracking, Big Sky, MT, 1-8 Mar. 2014.
%
%January 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

phi=x(4);
theta=x(5);
s=x(6);

sinPhi=sin(phi);
cosPhi=cos(phi);
sinTheta=sin(theta);
cosTheta=cos(theta);

sinPhiCosTheta=sinPhi*cosTheta;
cosPhiCosTheta=cosPhi*cosTheta;

%The coordinate system used in [1] is different from coordinate system 0 in
%spher2Cart: The sine and cosine terms in x and y are switched, meaning
%that effectively x and y are switched (or they define their elevation
%angle differently). Thus, to be consistent with coordinate system 0,  we
%switch what is in X and Y. The variable names below are consistent with
%that is used as X and Y in the paper, but the elements into which the
%values are placed have been reversed, to be consistent with coordinate
%system 0 in spher2Cart.
aVal=[s*cosPhiCosTheta;
      s*sinPhiCosTheta;
      s*sinTheta;
      0;
      0;
      0];
  
if(nargout>1)
    sinPhiSinTheta=sinPhi*sinTheta;
    cosPhiSinTheta=cosPhi*sinTheta;

    dXdPhi=s*cosPhiCosTheta;
    dXdTheta=-s*sinPhiSinTheta;
    dXds=sinPhiCosTheta;

    dYdPhi=-s*sinPhiCosTheta;
    dYdTheta=-s*cosPhiSinTheta;
    dYds=cosPhiCosTheta;

    dZdPhi=0;
    dZdTheta=s*cosTheta;
    dZds=sinTheta;
    
    aJacob=[0,0,0,dYdPhi,dYdTheta,   dYds;
            0,0,0,dXdPhi,dXdTheta,   dXds;
            0,0,0,dZdPhi,dZdTheta,   dZds;
            zeros(3,6)];

    if(nargout>2)
        aHess=zeros(6,6,6);
        
        dXdPhidPhi=-s*sinPhiCosTheta;
        dXdThetadPhi=-s*cosPhiSinTheta;
        dXdsdPhi=cosPhiCosTheta;
        
        dXdPhidTheta=dXdThetadPhi;
        dXdThetadTheta=-s*sinPhiCosTheta;
        dXdsdTheta=-sinPhiSinTheta;
        
        dXdPhids=dXdsdPhi;
        dXdThetads=dXdsdTheta;
        dXdsds=0;
        
        %%%
        dYdPhidPhi=-s*cosPhiCosTheta;
        dYdThetadPhi=s*sinPhiSinTheta;
        dYdsdPhi=-sinPhiCosTheta;
        
        dYdPhidTheta=dYdThetadPhi;
        dYdThetadTheta=-s*cosPhiCosTheta;
        dYdsdTheta=-cosPhiSinTheta;
        
        dYdPhids=dYdsdPhi;
        dYdThetads=dYdsdTheta;
        dYdsds=0;
        
        %%%
        dZdPhidPhi=0;
        dZdThetadPhi=0;
        dZdsdPhi=0;
        
        dZdPhidTheta=dZdThetadPhi;
        dZdThetadTheta=-s*sinTheta;
        dZdsdTheta=cosTheta;
        
        dZdPhids=dZdsdPhi;
        dZdThetads=dZdsdTheta;
        dZdsds=0;

        aHess(:,:,4)=[0,0,0,dYdPhidPhi,dYdThetadPhi,dYdsdPhi;
                      0,0,0,dXdPhidPhi,dXdThetadPhi,dXdsdPhi;
                      0,0,0,dZdPhidPhi,dZdThetadPhi,dZdsdPhi;
                      zeros(3,6)];
        aHess(:,:,5)=[0,0,0,dYdPhidTheta,dYdThetadTheta,dYdsdTheta;
                      0,0,0,dXdPhidTheta,dXdThetadTheta,dXdsdTheta;
                      0,0,0,dZdPhidTheta,dZdThetadTheta,dZdsdTheta;
                      zeros(3,6)];
        aHess(:,:,6)=[0,0,0,dYdPhids,dYdThetads,dYdsds;
                      0,0,0,dXdPhids,dXdThetads,dXdsds;
                      0,0,0,dZdPhids,dZdThetads,dZdsds;
                      zeros(3,6)];
        
        if(nargout>3)
            papt=zeros(6,1);
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
