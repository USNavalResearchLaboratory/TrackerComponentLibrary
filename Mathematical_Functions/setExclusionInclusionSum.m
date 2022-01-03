function sumVal=setExclusionInclusionSum(A,algorithm)
%%SETEXCLUSIONINCLUSIONSUM This evaluates the set inclusion-exclusion sum:
%              sumVal=sum_1 A_{i_1}-sum_2 A_{i_1}*A_{i_2}+sum_3
%                     A_{i_1}*A_{i_2}*A_{i_3}-...sum_{n}A_{i_1}*...*A_{i_n}
%              where sum_k is the sum over the set of i_1,...,i_k length-k
%              combinations of n things and the sign of each summation
%              group alternates.
%
%INPUTS: A An nX1 or 1XN vector. It can be real or complex.
% algorithm An optional parameter specifying the algorithm to use. Possible
%          values are:
%          0 Use the algorithm of [1], which involves combinatorics.
%          1 (The default if omitted or an empty matrix is passed) Just
%            compute 1-prod(1-A), which is equivalent to the above sum.    
%
%OUTPUTS: sumVal The desired sum value.
%
%The algorithm of [1] is implemented here. Note that if all of the values
%in A are probabilities of independent events, then sumVal is the
%probability of any event occurring which is equivalent to 1-prod(1-A). One
%notices that this identity applies regardless of the values of A (real or
%somplex), which is how we come to algorithm 1 (the default), which can be
%much faster than algorithm 0 on small problems.
%
%EXAMPLE 1:
%Here, we verify that the value for 4 random values equals the explicit
%solution given complex values of A.
% A=randn(4,1)+1j*randn(4,1);
% sol0=setExclusionInclusionSum(A,0);
% sol1=setExclusionInclusionSum(A,1);
% solExp=sum(A)-A(1)*A(2)-A(1)*A(3)-A(1)*A(4)-A(2)*A(3)-A(2)*A(4)-A(3)*A(4)+A(1)*A(2)*A(3)+A(1)*A(2)*A(4)+A(1)*A(3)*A(4)+A(2)*A(3)*A(4)-A(1)*A(2)*A(3)*A(4);
% RelErr0=abs(solExp-sol0)/abs(solExp)
% RelErr1=abs(solExp-sol1)/abs(solExp)
%The relative errors of both solutions should be within finite precision
%limits.
%
%EXAMPLE 2:
%We verify that both algorithms produce the same result, within finite
%precision limits.
% A=randn(4,1)+1j*randn(4,1);
% sol0=setExclusionInclusionSum(A,0);
% sol1=setExclusionInclusionSum(A,1);
% RelErr=abs(sol0-sol1)/abs(sol0)
%
%REFERENCES:
%[1] C. J. Mifsud, "Algorithm 156: Algebra of sets," Communications of the
%    ACM, vol. 6, no. 3, pg. 103, Mar. 1963.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
   algorithm=1; 
end

switch(algorithm)
    case 0
        n=length(A);

        sumVal=0;
        j=-1;

        %%%%B
        for r=1:n
            part=0;
            I=[];

            while(1)%%%A
                if(isempty(I))
                    I=1:r;
                else
                    I=getNextCombo(I,n,r);
                end

                if(isempty(I))
                    j=-j;
                    part=j*part;
                    sumVal=sumVal+part;
                    break;
                else%I is not empty
                    T=1;
                    for i=1:r
                        T=A(I(i))*T;
                    end
                    part=part+T;
                end
            end
        end
    case 1
        sumVal=1-prod(1-A);
    otherwise
        error('Unknown algorithm specified.')
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
