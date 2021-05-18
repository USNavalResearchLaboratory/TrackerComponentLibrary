function aHiMed=weightedMedianHi(a,w,algorithm)
%WEIGHTEDMEDIANHI Compute the weighted high (upper) median. This finds the
%                 smallest value a(k) such that the sum of the weights
%                 associated with all a(i)<=a(k) is greater than half the
%                 total weight.
%
%INPUTS: a An nX1 or 1Xn real vector.
%        w A 1Xn or nX1 real vector of weights associated with the elements
%          of a. All w(i)>=0.
% algorithm An optional parameter specifying the algorithm to use. Possible
%          values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            efficient algorithm given in [1]. This is an O(n) algorithm.
%          1 Sort the values and just go through w until the cumulative
%            value of w is > wHalf. The sorting means that this will not be
%            better than O(n*log(n)).
%
%OUTPUTS: aHiMed The value of the high weighted median in a. If a contains
%                a NaN, then an empty matrix will be returned.
%
%EXAMPLE:
%Here, we demonstrate on random problems that both approaches give the same
%result. When run, the assertion should always be true (no error occur).
% numEls=50;
% numRuns=10000;
% for curRun=1:numRuns
%     a=fix(10*rand(numEls,1));
%     w=fix(10*rand(numEls,1));
%     aHiMed0=weightedMedianHi(a,w,0);
%     aHiMed1=weightedMedianHi(a,w,1);
%     assert(aHiMed0==aHiMed1)
% end
%
%REFERENCES:
%[1] C. Croux and P. J. Rousseeuw, "Time-efficient algorithms for two-
%    highly robust estimators of scale," in Computational Statistics.
%    Berlin: Springer-Verlag, Aug. 1992, vol. 1: Proceedings of the 10th
%    Symposium on Computational Statistics, pp. 411-428, Conference
%    Location: Neuchâtel, Switzerland.
%
%August 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
    algorithm=0;
end

if(any(isnan(a)))
   aHiMed=[];
   return;
end

switch(algorithm)
    case 0
        nOrig=length(a);

        %Allocate space.
        aCand=zeros(nOrig,1);
        wCand=zeros(nOrig,1);

        wHalf=sum(w)/2;
        wRest=0;
        n=nOrig;
        while(1)
            %This is something of a pivot point.
            trial=kthOrderStat(a(1:n),fix(n/2)+1);

            %These sums could all be done in one loop, but separating them is
            %probably faster in Matlab.
            wLeft=sum(w(1:n).*(a(1:n)<trial));
            wMid=sum(w(1:n).*(a(1:n)==trial));
            %Now that we could also find:
            %wRight=sum(w(1:n).*(a(1:n)>trial));

            if(wRest+wLeft>wHalf)
                %If there is too much weight in the left side.
                kCand=0;
                for i=1:n
                    if(a(i)<trial)
                        kCand=kCand+1;
                        aCand(kCand)=a(i);
                        wCand(kCand)=w(i);
                    end
                end
                n=kCand;
            elseif(wRest+wLeft+wMid>wHalf)
                aHiMed=trial;
                return;
            else
                kCand=0;
                for i=1:n
                    if(a(i)>trial)
                       kCand=kCand+1;
                       aCand(kCand)=a(i);
                       wCand(kCand)=w(i);
                    end
                end
                n=kCand;
                %We are essentially taking the higher-branch, so the points
                %in the lower part that will no longer be considered are in
                %wRest.
                wRest=wRest+wLeft+wMid;
            end

            a(1:n)=aCand(1:n);
            w(1:n)=wCand(1:n);
        end
    case 1
        n=length(a);
        [a,idx]=sort(a);
        w=w(idx);
        wHalf=sum(w)/2;

        wCum=0;
        for i=1:n
            wCum=wCum+w(i);
            if(wCum>wHalf)
                aHiMed=a(i);
                return;
            end
        end
        
        %The high median does not exist. This should not occur.
        aHiMed=[];
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
