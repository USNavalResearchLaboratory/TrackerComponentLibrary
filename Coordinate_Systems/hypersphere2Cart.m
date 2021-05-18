function x=hypersphere2Cart(rPhiList)
%%HYPERSPHERE2CART Given a point in n-dimensional spherical coordinates
%                  (hyperspherical coordinates), convert the point to
%                  Cartesian coordinates.
%
%INPUTS: rPhiList An nXN list of N values to convert from hyperspherical to
%                 Cartesian coordinates. The values are in n-dimensional
%                 space (n>=2). Each value is a range r (0<=r<Inf) followed
%                 by (n-1) angular values. The first (n-2) values are bound
%                 from 0 to pi and the final one can span the entire circle
%                 (-pi to pi). In 2D, the angle is just from (-pi to pi).
%
%OUTPUTS: x An nXN matrix of the values converted into n-dimensional
%           Cartesian coordinates.
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
%spher2Cart function.
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
    
    N=size(rPhiList,2);
    n=size(rPhiList,1);
    
    x=zeros(n,N);
    for curEl=1:N
        r=rPhiList(1,curEl);
        phiList=rPhiList(2:end,curEl);

        cosList=cos(phiList);
        sinProdList=cumprod(sin(phiList));

        x(1,curEl)=r*cosList(1);
        for curDim=2:(n-1)
            x(curDim,curEl)=r*sinProdList(curDim-1)*cosList(curDim);
        end
        x(n,curEl)=r*sinProdList(end);
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
