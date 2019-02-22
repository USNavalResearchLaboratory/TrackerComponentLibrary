function points=HaltonSequence(numDim,numPoints,algorithm,skipFirst,stepSize)
%%HALTONSEQUENCE Obtain values of the Halton sequence. This is a low-
%           discrepance multivariate sequency with values ranging from 0 to
%           1. Such sequence values can be used for multivariate quasi-
%           Monte Carlo integration, which can have better performance than
%           Monte Carlo integration. The Halton sequence has van der Corput
%           sequences across the dimensions, where each is a different
%           prime base and where each can be differently scrambled.
%
%INPUTS: numDim The number of dimensions of the point set to generate.
%               numDim>=1.
%     numPoints The number of quasi-random ponts to generate.
%     algorithm This selects how the elements of the different van der
%               Corput sequences that make up the Halton sequence are
%               permuted (scrambled). See below for what what means. Note
%               that each sequences has a different prime base. The default
%               algorithm chosen if this parameter is omitted or an empty
%               matrix is passed is shown below. Possible values are:
%               0 Do not scramble the sequences.
%               1 Use the Faure scrambling of Section 3.1 of [1] (The
%                 default if numDim>50).
%               2 Use the "optimal" scrambling of [2]. This is only valid
%                 for numDim<=50 (The default if 16<numDim<=50).
%               3 Use the Atanassov's scrambling with the parameters given
%                 in Section 4.6 (Table 3) of [3].
%               4 Use the Braaten-Weller scrambling of [4]. This is only
%                 valid for numDim<=16 (The default if numDim<=16).
%               5 Use the backwards scrambling of [3]. Despite the claims
%                 in [3], this is a rather bad scrambling.
%     skipFirst The number of values at the start of the sequence to skip
%               before recording values in points (>=0). The default if
%               omittedor an empty matrix is passed is 0. This can be a
%               scalar if it is the same for all dimensions, or it can be a
%               length numDim vector, if one wishes it to be different in
%               different dimensions.
%      stepSize Rather than taking every point in the sequence, every
%               stepSize point is taken. The default if omitted or an empty
%               matrix is passed is 1, meaning that every point is taken.
%               This can be a scalar if it is the same for all dimensions,
%               or it can be a length numDim vector, if one wishes it to be
%               different in different dimensions.
%
%OUTPUTS: points The numDimXnumPoints set of values of the Halton sequence
%                for the given parameters. Each entry is between 0 and 1.
%
%The scalar value of the ith entry in a van der Corput sequence S(n) is
%S(n)=sum_{j=0}^Inf sigma(a(j))*base^{-j-1}
%where the vector a is the digits of the number n converted into the
%specified base and sigma is a permutation of the values of 0:(base-1) that
%begins with 0. The Halton sequence is composed of multiple van der Corput
%sequences across the different dimensions of a multidimensional problem.
%Each sequence has a base that is a different consecuritive prime number,
%here starting with 2. The permutation sigma for each van der Corput
%sequence in the Halton sequence is determined by the selected algorithm.
%
%EXAMPLE:
%In higher dimensions, the unscrambled Halton sequence can perform poorly.
%Here, we compare the values of algorithms 0 and 4 in dimensions 7 and 8 of
%an 8D set of points.
% points0=HaltonSequence(8,1000,0);
% points4=HaltonSequence(8,1000,4);
% figure(1)
% clf
% hold on
% scatter(points0(7,:),points0(8,:),'.r')
% scatter(points4(7,:),points4(8,:),'.b')
% h1=xlabel('Dimension 7');
% h2=ylabel('Dimension 8');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% legend('Algorithm 0','Algorithm 4')
%One can see that the unscrambled Halton points (red) are arranged in a
%linear pattern to a greater extent than those scrambled with algorithm 4.
%
%REFERENCES:
%[1] H. Faure, "Good permutations for extreme discrepancy," Journal of
%    Number Theory, vol. 42, no. 1, pp. 47?56, Sep. 1992.
%[2] H. Chi, M. Mascagni, and T. Warnock, "On the optimal Halton sequence,"
%    Mathematics and Computers in Science, vol. 70, no. 1, pp. 9-21, 1 Sep.
%   2005.
%[3] B. Vandewoestyne and R. Cools, "Good permutations for deterministic
%    scrambled Halton sequences in terms of L2-discrepancy," Journal of
%    Computational and Applied Mathematics, vol. 189, no. 1?2, pp. 341-361,
%    May 2006.
%[4] E. Braaten and G. Weller, "An improved low-discrepancy sequence for
%    multidimensional quasi-monte carlo integration," Journal of
%    Computational Physics, vol. 33, pp. 249?258, Nov. 1979.
%
%January 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(stepSize))
    stepSize=1;
end

if(nargin<5||isempty(skipFirst))
    skipFirst=1;
end

if(nargin<3||isempty(algorithm))
    if(numDim<=16)
        algorithm=4;
    elseif(numDim<=50)
        algorithm=2;
    else
        algorithm=1;
    end
end

if(isscalar(stepSize))
    stepSize=stepSize*ones(numDim,1);
end

if(isscalar(skipFirst))
    skipFirst=skipFirst*ones(numDim,1);
end

b=firstNPrimes(numDim);

points=zeros(numDim,numPoints);
switch(algorithm)
    case 0%Do not scramble.
        for curDim=1:numDim
            sigma=0:(b(curDim)-1);
            points(curDim,:)=vanDerCorputSeq([1,numPoints],skipFirst(curDim),stepSize(curDim),b(curDim),sigma);
        end
    case 1%Use the Faure scrambling.
        for curDim=1:numDim
            points(curDim,:)=vanDerCorputSeq([1,numPoints],skipFirst(curDim),stepSize(curDim),b(curDim));
        end
    case 2%Use the "optimal" scrambling from [2].
        if(numDim>50)
            error('This algorithm is not available for numDim>50.')
        end
        %For the first 50 primes.
        w=[1,2,2,5,3,7,3,10,18,11,17,5,17,26,40,14,40,44,12,31,45,70,8,38,82,8,12,...
            38,47,70,29,57,97,110,32,48,84,124,155,26,69,83,157,171,8,32,112,205,15,31];

        for curDim=1:numDim
            sigma=mod(w(curDim)*(0:(b(curDim)-1)),b(curDim));
            points(curDim,:)=vanDerCorputSeq([1,numPoints],skipFirst(curDim),stepSize(curDim),b(curDim),sigma);
        end
    case 3%Use Atanassov's scrambling as given in [3]. 
        if(numDim>30)
            error('This algorithm is not available for numDim>30.')
        end
        %For the first 30 primes.
        kVals=[1 1 4 2 9 9 2 1 13 6 22 7 37 36 36 39,4 26 13 12 35 66 60 68 63 47 15 104 4 64];

        kPows=(0:(b(numDim)-1)).^(1:(b(numDim)));
        for curDim=1:numDim
            sigma=mod(kVals(curDim)*kPows(1:b(numDim)),b(curDim));

            points(curDim,:)=vanDerCorputSeq([1,numPoints],skipFirst(curDim),stepSize(curDim),b(curDim),sigma);
        end
    case 4%The Braaten-
        if(numDim>30)
            error('This algorithm is not available for numDim>16.')
        end
        sigmaVals=cell(16,1);
        sigmaVals{1}=[0,1];
        sigmaVals{2}=[0,2,1];
        sigmaVals{3}=[0,3,1,4,2];
        sigmaVals{4}=[0,4,2,6,1,5,3];
        sigmaVals{5}=[0,5,8,2,10,3,6,1,9,7,4];
        sigmaVals{6}=[0,6,10,2,8,4,12,1,9,5,11,3,7];
        sigmaVals{7}=[0,8,13,3,11,5,16,1,10,7,14,4,12,2,15,6,9];
        sigmaVals{8}=[0,9,14,3,17,6,11,1,15,7,12,4,18,8,2,16,10,5,13];
        sigmaVals{9}=[0,11,17,4,20,7,13,2,22,9,15,5,18,1,14,10,21,6,16,3,19,8,12];
        sigmaVals{10}=[0,15,7,24,11,20,2,27,9,18,4,22,13,26,5,16,10,23,1,19,28,6,14,17,3,25,12,8,21];
        sigmaVals{11}=[0,15,23,5,27,9,18,2,29,12,20,7,25,11,17,3,30,14,22,1,21,8,26,10,16,28,4,19,6,24,13];
        sigmaVals{12}=[0,18,28,6,23,11,34,3,25,14,31,8,20,36,1,16,27,10,22,13,32,4,29,17,7,35,19,2,26,12,30,9,24,15,33,5,21];
        sigmaVals{13}=[0,20,31,7,26,12,38,3,23,34,14,17,29,5,40,10,24,1,35,18,28,9,33,15,21,4,37,13,30,8,39,22,2,27,16,32,11,25,6,36,19];
        sigmaVals{14}=[0,21,32,7,38,13,25,3,35,17,28,10,41,5,23,30,15,37,1,19,33,11,26,42,8,18,29,4,39,14,22,34,6,24,12,40,2,31,20,27,9,36,16];
        sigmaVals{15}=[0,24,12,39,6,33,20,44,3,29,16,36,10,42,22,8,31,26,14,46,1,35,18,28,5,40,19,37,11,25,43,4,30,15,34,9,45,21,2,32,17,41,13,27,7,38,23];
        sigmaVals{16}=[0,26,40,9,33,16,49,4,36,21,45,12,29,6,51,23,38,14,43,1,30,19,47,10,34,24,42,3,27,52,15,18,39,7,46,31,11,35,20,48,2,28,41,8,22,50,13,32,17,44,5,37,25];
        
        for curDim=1:numDim
            points(curDim,:)=vanDerCorputSeq([1,numPoints],skipFirst(curDim),stepSize(curDim),b(curDim),sigmaVals{curDim});
        end
    case 5%Use the backward scrambling of [1] (as bad as not scrambling).
        for curDim=1:numDim
            sigma=[0,b(curDim)-(1:(b(curDim)-1))];
            points(curDim,:)=vanDerCorputSeq([1,numPoints],skipFirst(curDim),stepSize(curDim),b(curDim),sigma);
        end
    otherwise
        error('Unknown algorithm selected.')
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
