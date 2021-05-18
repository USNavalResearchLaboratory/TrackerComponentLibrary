function S=tria(A)
%%TRIA Square root matrix triangularization. Given a rectangular square
%      root matrix, obtain a lower-triangular square root matrix that is
%      square.
%
%INPUTS: A A numRowXnumCol matrix that is generally not square.
%
%OUTPUTS: S A lower-triangular matrix such that S*S'=A*A'. If
%           numCol>=numRow, then S is a square numRowXnumRow matrix.
%           Otherwise, S is a numRowXnumCol matrix.
%
%This is the tria function needed for various steps in the cubature Kalman
%filter and the square root Kalman filter. It is described in [1]. It has
%been slightly modified from the paper so that the diagonal elements remain
%positive.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%July 2012 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    [~,R]=qr(A',0);
    S=R';
    
    %Make the diagonal elements all positive.
    sel=diag(S)<0;
    S(:,sel)=-S(:,sel);
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
