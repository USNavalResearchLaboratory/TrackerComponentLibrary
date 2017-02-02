function R=rotZAxis2Vec(u)
%%ROTZAXIS2VEC  Get a rotation matrix that rotates the z-axis in the
%               direction of a given vector u. That is u=R*z. As the
%               solution to the problem is not unique, the
%               shortest-distance rotation is chosen.
%
%INPUTS: u1     A 3X1 unit vector representing the direction into which the
%               z-axis is to be rotated.
%
%OUTPUTS: R     A rotation matrix such that u=R*z, where z is the z-axis
%               [0;0;1].
%
%In Appendix C of [1], the minimal rotation to point the z-axis in a
%desired direction is derived.
%
%REFERENCES:
%[1] D. F. Crouse, "On measurement-based light-time corrections for
%    bistatic orbital debris tracking," IEEE Transactions on Aerospace and
%    Electronics Systems, vol. 51, no. 3, pp. 2502-2518,  Jul. 2015.
%
%August 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    s=-sqrt(0.5*(1-u(3)));
    dist=sqrt(u(1)^2+u(2)^2);
    q=[sqrt(1-s^2),-s*u(2)/dist,s*u(1)/dist,0];
    R=quat2RotMat(q,'left');
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
