function rphiList=Cart2Hypersphere(x)
%%HYPERSPHERE2CART Given a point in n-dimensional Cartesian coordinates
%                  convert the point to n-dimensional spherical coordinates
%                  (hyperspherical coordinates.
%
%INPUTS: x An nXN matrix of the values converted into n-dimensional
%          Cartesian coordinates. n>=2.
%
%OUTPUTS: rPhiList The nXN values convert from Cartesian coordinates to 
%                 hyperspherical coordinates. Each value is a range r
%                 (0<=r<Inf) followed by (n-1) angular values. The first
%                 (n-2) values are bound from 0 to pi and the final one can
%                 span the entire circle (-pi to pi).
%
%Hyperspherical coordinates are derived in [1]. The definition here is
%nearly identical, except the final two coordinates are switched (the final
%term is all sines). Hyperspherical coordinates arise in tracking, for
%example, when generating false measurements in gates in simulations, such
%as in [2]. Note that in [2], yet another definition of the system is used.
%
%The definition here is 
% x(1)=cos(phi(1));
% x(2)=sin(phi(1))*cos(phi(2));
% x(3)=sin(phi(1))*sin(phi(2))*cos(phi(3));
% ...
% x(n-1)=sin(phi(1))*...*sin(phi(n-2))*cos(phi(n-1));
% x(n)=sin(phi(1))*...*sin(phi(n-2))*sin(phi(n-1));
%
%The definition of the angles differs in 3D from those use in the
%Cart2Sphere function.
%
%REFERENCES:
%[1] L. E. Blumenson, "A derivation of n-dimensional spherical
%    coordinates," The American Mathematical Monthly, no. 1, pp. 63-66,
%    Jan. 1960.
%[2] H. Sun and M. Farooq, "Note on the generation of random points
%    uniformly distributed in hyper-ellipsoids," in Proceedings of the
%    Fifth International Conference on Information Fusion, Annapolis, MD,
%    8-11 Jul. 2002, pp. 489-496.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(x);

phiList=zeros(n-1,1);
r=norm(x);

for curDim=1:(n-1)
    normVal=norm(x(curDim:n));
    
    if(normVal==0)
        if(x(curDim)>=0)
            phiList(curDim)=0;
        else
            phiList(curDim)=pi;
        end
    else
        phiList(curDim)=acos(x(curDim)/normVal);
    end
end

if(x(n)<0)
    phiList(n-1)=2*pi-phiList(n-1);
end
phiList(n-1)=wrapRange(phiList(n-1),-pi,pi);

rphiList=[r;phiList];
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
