function z=genUniformGridOnSphere(numPts)
%%GENUNIFORMGRIDONSPHERE Generate a uniform grid of Cartesian points on the
%                        surface of a unit sphere. 
%
%INPUTS: numPts This is the desired number of points to obtain on the
%               sphere. The actual number of points obtained can be
%               different.
%
%OUTPUTS: z A 3XnumGenPoints set of 3D Cartesian points approximately
%           uniformly distributed on the 3D unit sphere.
%
%This function implements the deterministic algorithm of [1]. 
%
%EXAMPLE:
%This example generates 500 points on the unit sphere and then display
%them, so that one can see that the points are approximately uniform on the
%sphere.
% pts=genUniformGridOnSphere(500);
% scatter3(pts(1,:),pts(2,:),pts(3,:))
%
%REFERENCES:
%[1] M. Deserno, "How to generate equidistributed points on the surface of
%    a sphere," Max-Planck-INntitut fuer Polymerforschung, Tech. Rep., 28
%    Sep. 2004. [Online]. Available:
%    https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
%
%October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<1||isempty(numPts))
    numPts=1;
end

z=zeros(3,ceil(1.2*numPts));

NCount=0;
a=4*pi/numPts;
d=sqrt(a);
MTheta=round(pi/d);
dTheta=pi/MTheta;
dPhi=a/dTheta;

for m=0:(MTheta-1)
    theta=pi*(m+(1/2))/MTheta;
    MPhi=round(2*pi*sin(theta)/dPhi);  
   
    sinTheta=sin(theta);
    cosTheta=cos(theta);
    for n=0:(MPhi-1)
        phi=2*pi*n/MPhi;
        NCount=NCount+1;
        
        z(:,NCount)=[sinTheta*cos(phi);
                     sinTheta*sin(phi);
                     cosTheta];  
    end
end

%Size to fit.
z=z(:,1:NCount);

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

