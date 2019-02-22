function points=vanDerCorputSeq(N,skipFirst,stepSize,base,sigma)
%%VANDERCORPUTSEQ Obtain values of the van der Corput sequence. This is a
%           low-discrepance 1D sequency with values ranging from 0 to 1.
%           Such sequence values can be used for 1D quasi-Monte Carlo
%           integration, which can have better performance than Monte Carlo
%           integration. For multivariate quasi-Monte Carlo integration
%           problems, other sequences, such as the Halton sequence, should
%           be used as this will perform poorly.
%
%INPUTS: N A vector indicating the size of the point set to generate. N(i)
%          is the size of dimension i of points. For example, for a 10X1
%          vector, one should set N=[10;1]. If a scalar is passed, then a
%          matrix of size [N,N] is generated.
% skipFirst The number of values at the start of the van der Corput
%          sequence to skip before recording values in points. The default
%          if omitted or an empty matrix is passed is 0.
% stepSize Rather than taking every point in the sequence, every stepSize
%          point is taken. The default if omitted or an empty matrix is
%          passed is 1, meaning that every point is taken.
%     base Different sequences are defined for different bases (integer
%          values >=2). This is the base to use. If this parameter and
%          sigma are both omitted, then base=36 and the good permutation
%          sigma for base 36 that is given in Theorem 1.2 of [1] is used.
%          If base is omitted but sigma is given, then the value implied by
%          sigma is used.
%    sigma Improved performance can be obtained by permuting the order of
%          the digits creating the terms of the sequence, as given in the
%          equation below. This is a baseX1 or 1Xbase permutation of the
%          digits 0:(base-1) that starts with 0. If this parameter is
%          omitted and base is provided, then the good permutation that is
%          given in Section 3.1 of [1] is used. If base is omitted, then
%          base and sigma are both set as described for the base input.
%
%OUTPUTS: points The set of values of the van der Corput sequence for the
%                given parameters. These are values between 0 and 1.
%
%The value of the ith entry in a van der Corput sequence S(n) is
%S(n)=sum_{j=0}^Inf sigma(a(j))*base^{-j-1}
%where the vector a is the digits of the number n converted into the
%specified base and sigma is a permutation of the values of 0:(base-1) that
%begins with 0.
%
%EXAMPLE:
%A plot of a histogram of the points indicates that they can be relatively
%well distributed between 0 and 1 even though they are not pseudo-random.
% vals=vanDerCorputSeq([1000,1]);
% figure(1)
% clf
% histogram(vals,'normalization','probability')
% h1=xlabel('Value');
% h2=ylabel('Fraction of Samples');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] H. Faure, "Good permutations for extreme discrepancy," Journal of
%    Number Theory, vol. 42, no. 1, pp. 47?56, Sep. 1992.
%
%January 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(skipFirst))
    skipFirst=0;
end

if(nargin<3||isempty(stepSize))
    stepSize=1;
end

if((nargin<4||isempty(base))&&(nargin<5||isempty(sigma)))
    base=36;
    sigma=[0;25;17;7;31;11;20;3;27;13;34;22;5;15;29;9;23;1;18;32;8;...
              28;14;4;21;33;12;26;2;21;10;30;6;16;24;35];
elseif((nargin<4||isempty(base))&&~(nargin<5||isempty(sigma)))
    base=length(sigma);
end

if(nargin<5||isempty(sigma))
    sigma=makeFaurePerm(base);
else
    sigma=sigma(:);
end

if(isscalar(N))
    numPoints=N^2;
else
    numPoints=prod(N(:));
end

%Get the powers of the base.
[~,baseValsInt]=intVal2Base(0,base,false,0);
baseVals=double(baseValsInt)*base;
numBase=length(baseVals);

points=zeros(N);
for k=1:numPoints
    rank=skipFirst+(k-1)*stepSize+1;

    digits=intVal2Base(rank,base,false,numBase,baseValsInt);
    digits=sigma(digits+1);
    points(k)=sum(digits./baseVals);
end
end

function sigma=makeFaurePerm(b)
%%MAKEFAUREPERM This generates the permutation suggested by Faure for
%               improving the van der Corput sequence, given in Section 3.1 
%               of [1]. The expressions in [1] can be a bit difficult to
%               understand. the formulation given in Section 4.3 of [2] is
%               clearer.
%
%REFERENCES:
%[1] H. Faure, "Good permutations for extreme discrepancy," Journal of
%    Number Theory, vol. 42, no. 1, pp. 47?56, Sep. 1992.
%[2] B. Vandewoestyne and R. Cools, "Good permutations for deterministic
%    scrambled Halton sequences in terms of L2-discrepancy," Journal of
%    Computational and Applied Mathematics, vol. 189, no. 1?2, pp. 341-361,
%    May 2006.
%
%January 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    %We have to generate subsequent permutations. These are stored in a
    %table. The permutations have increasing sizes, but it is faster in
    %matlab to use a table with a bunch of unused elements than to use cell
    %arrays with entries of different sizes.
    sigmaTable=zeros(b,b);
    sigmaTable(1:2,2)=[0;1];
    for bCur=3:b
        if(mod(bCur,2)==0)
            k=bCur/2;
            sigmaTable(1:bCur,bCur)=[2*sigmaTable(1:k,k);2*sigmaTable(1:k,k)+1];
        else
            eta=sigmaTable(1:(bCur-1),bCur-1);
            k=(bCur-1)/2;
            sel=(eta>=k);
            eta(sel)=eta(sel)+1;
            sigmaTable(1:bCur,bCur)=[eta(1:k);k;eta((k+1):(bCur-1))];
        end
    end
    
    sigma=sigmaTable(1:b,b);
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
