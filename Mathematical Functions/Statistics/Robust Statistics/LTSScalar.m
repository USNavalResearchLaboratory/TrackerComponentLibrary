function [thetaMin,minCost]=LTSScalar(y)
%%LTSSCALAR Find the least trimmed squares estimate of a set of scalar
%           samples. This is the value
%           theta=min_theta \sum_{i=1}^h r_i^2
%           where r.^2=sort(abs(y-theta)^2,'ascend). That is, the cost
%           function is the sum of the smallest h residuals. h is taken to
%           be floor(n/2)+1.
%
%INPUTS: y An NX1 or 1XN vector of N real values.
%
%OUTPUTS: theta The least trimmed squares estimate of the given scalar
%               samples.
%
%The recursive formula described on page 172 of Chapter 4-2 of [1] is
%implemented here. 
%
%EXAMPLE:
%This is the example from page 165 of [1] (Chapter 4-2). The point 40 is an
%outlier.
% y=[90,93,86,92,95,83,75,40,88,80];
% theta=LTSScalar(y)
% meanVal=mean(y)
%One will see that theta is about 90.6667. The mean (the least squares
%solution) is drawn towards the outlier and is 82.2.
%
%REFERENCES:
%[1] P. J. Rousseeuw and A. M. Leroy, Robust Regression and Outlier
%    Detection. New York: John Wiley and Sons, 1987.
%
%August 2018 David F. Crouse, Naval Research Laboratory, Washington D.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(y);
h=floor(n/2)+1;

y=sort(y);

thetaMin=mean(y(1:h));
minCost=sum((y(1:h)-thetaMin).^2);

thetaPrev=thetaMin;
costPrev=minCost;
for j=2:(n-h+1)
    thetaCur=(h*thetaPrev-y(j-1)+y(j+h-1))/h;
    costCur=costPrev-y(j-1)^2+y(j+h-1)^2-h*thetaCur^2+h*thetaPrev;
    
    if(costCur<minCost)
        minCost=costCur;
        thetaMin=thetaCur;
    end
    
    thetaPrev=thetaCur;
    costPrev=costCur;
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
