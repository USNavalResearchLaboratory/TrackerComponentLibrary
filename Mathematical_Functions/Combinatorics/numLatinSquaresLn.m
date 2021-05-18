function logNumSquares=numLatinSquaresLn(n,isReduced)
%%NUMLATINSQUARESLN Return the natural logarithm of the number of Latin
%           squares when the number is known exactly (for 1<=n<=11), or
%           return bounds on the number, when the number is not known
%           exactly (n>=12). This number can either be the total number of
%           latin squares, or the reduced number, as described below. The
%           number of Latin squares is related to the number of feasible
%           solution to a planar 3D assignment problem and it increases
%           very rapidly with n.
%
%INPUTS: n The scalar size of the square, n>=1.
% isReduced A boolean variable indicating whether or not the number of
%          reduced square or all squares is desired. A reduced latin
%          square is one where the values in the first row and the first
%          column have been permuted so as to be in increasing order. A
%          latin square in reduced form is also said to be in normalized or
%          standard form. The default if this parameter is omitted or an
%          empty matrix is passed is false.
%
%OUTPUTS: logNumSquares If n<=11, then this is the natural logarithm of the
%                   exact scalar number of nXn latin squares. If n>11,
%                   then this is a 2X1 vector where logNumSquares(1) is a
%                   lower bound and logNumSquares(2) is an upper bound on
%                   the natural logarithm of the number of possible latin
%                   squares.
%
%A Latin square is an nXn matrix where each row and each column contain the
%values 1:n exactly once. For example,
% [1 2 3;
%  3 1 2;
%  2 3 1]
%is a valid latin square.
%
%The exact number of latin squares for n=1 to 11 are tabulated in [1] and
%the tabulated values are used here. An exact algorithm is given in [2]
%(and a different one in [1]), but the computational complexity of the
%algorithms increases so rapidly with n that they are not really practical.
%
%For n>11, the upper and low bounds given in Ch. 17 of [3] are used. Note
%that the upper bound is based on van der Waerden's conjecture, which is
%proven in [4].
%
%REFERENCES:
%[1] B. D. McKay and I. M. Wanless, "On the number of latin squares,"
%    Annals of Combinatorics, vol. 9, no. 3, pp. 335-344, Oct. 2005.
%[2] J.-y. Shao and W.-d. Wei, "A formula for the number of latin squares,"
%    Discrete Mathematics, vol. 110, no. 1-3, pp. 293-296, 11 Dec. 1992.
%[3] J. H. van Lint and R. M. Wilson, A Course in Combinatorics. York:
%    Cambridge University Press, 1992.
%[4] B. Gyires, "Elementary proof for a van der Waerden's conjecture and
%    related theorems," Computational Mathematics with Applications, vol.
%    31, no. 10, pp. 7-21, May 1996.
%
%October 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(isReduced))
    isReduced=false; 
end

%Exact R and L values go from n=1 to 11.
R= [1;
    1;
    1;
    4;
    56;
    9408;
    16942080;
    535281401856;
    377597570964258816;
    7580721483160132811489280;
    5363937773277371298119673540771840];

L= [1;
    2;
    12;
    576;
    161280;
    812851200;
    61479419904000;
    108776032459082956800;
    5524751496156892842531225600;
    9982437658213039871725064756920320000;
    776966836171770144107444346734230682311065600000];

if(isempty(n))
    logNumSquares=[];
    return
end

if(~isscalar(n)||~isreal(n)||n<1)
    error('n is invalid.')
end

if(n<=11)
    if(isReduced)
        logNumSquares=log(R(n));
    else
        logNumSquares=log(L(n));
    end
else
    if(isReduced)
        deltaVal=gammaln(n+1)+gammaln(n);
        LB=-deltaVal;
        UB=-deltaVal;
    else
        LB=0;
        UB=0;
    end

    LB=LB+(2*n)*gammaln(n+1)-n^2*log(n);
    for k=1:n
        UB=UB+(n/k)*gammaln(k+1);
        if(~isfinite(UB))
            break;
        end
    end

    logNumSquares=[LB;UB];
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
