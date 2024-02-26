function val=numberOfHyperplanePartitions(N,k)
%%NUMBEROFHYPERPLANEPARTITIONS Given a set of N points in k-dimensional
%       space, this is the number of unique partitions that can be formed
%       from a (k-1) dimensional hyperplane splitting the data into two
%       partitions. Thus, in 2D, this is the number of ways that the points
%       can be split using a line and in 3D, the number of ways that the
%       points can be divided using a plane.
%
%INPUTS: N The total number of points.
%        k The number of dimensions of space (e.g. 2D, 3D).
%
%This implements the formula given in [1], with just simplified expressions
%for N<=6.
%
%REFERENCES:
%[1] E. F. Harding, "The number of partitions of a set of n points in k
%    dimensions induced by hyperplanes," Proceedings of the Edinburgh
%    Mathematical Society, vol. 15, no. 4, pp. 285-289, Dec. 1967.
%
%June 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(N-1<=k)
    val=2^(N-1);
elseif(k<=6)
    switch(k)
        case 1
            val=N;
        case 2
            val=1+1/2*(N-1)*N;
        case 3
            val=1/6*N*(8+(N-3)*N);
        case 4
            val=1+1/24*(-1+N)*N*(18+(-5+N)*N);
        case 5
            val=1/120*N*(184+N*(-110+N*(55+(-10+N)*N)));
        case 6
            val=1/120*N*(184+N*(-110+N*(55+(-10+N)*N)))+binomial(N-1,6);
        otherwise
            error('Invalid k.')
    end
else
    val=0;
    for i=1:k
        val=val+binomial(N-1,i);
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
