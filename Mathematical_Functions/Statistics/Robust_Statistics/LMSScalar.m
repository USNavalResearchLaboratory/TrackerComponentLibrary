function theta=LMSScalar(y)
%%LMSSCALAR Find the least median of squares estimate of a set of scalar
%           samples. This is the value
%           theta=min_theta \median_i (y(i)-theta)^2
%           Compared to the mean, this estimate is robust to outliers.
%
%INPUTS: y An NX1 or 1XN vector of N real values.
%
%OUTPUTS: theta The least median of squares estimate of the given scalar
%               samples.
%
%The algorithm described in page 169 of [1] (Chapter 4-2) based on Theorem
%1 on page 166 is used.
%
%Theorem 2 on page 170 of [1] notes that each solution to 
% theta=min_theta \median_i (y(i)-theta)^2
%is a solution to 
%theta=min_theta \median_i abs(y(i)-theta)
%However, the latter problem can have many solutions.
%
%EXAMPLE:
%This is the example from page 165 of [1] (Chapter 4-2). The point 40 is an
%outlier.
% y=[90,93,86,92,95,83,75,40,88,80];
% theta=LMSScalar(y)
% meanVal=mean(y)
%One will see that theta is 90.5. The mean (the least squares solution) is
%drawn towards the outlier and is 82.2.
%
%REFERENCES:
%[1] P. J. Rousseeuw and A. M. Leroy, Robust Regression and Outlier
%    Detection. New York: John Wiley and Sons, 1987.
%
%August 2018 David F. Crouse, Naval Research Laboratory, Washington D.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(y);
h=floor(n/2)+1;

y=sort(y,'ascend');

minCost=Inf;
theta=[];
for offset=1:(n-h+1)
    cost=(y(h+offset-1)-y(offset))^2;
    
    if(cost<minCost)
        minCost=cost;
        theta=(y(h+offset-1)+y(offset))/2;
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
