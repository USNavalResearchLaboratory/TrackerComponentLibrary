function probVal=bivarGaussRectangleCDF(rectMin,rectMax,mu,R)
%%BIVARGAUSSRECTANGLECDF Evaluate the an integral of the probability
%      density function of the bivariate normal distribution with a
%      specified mean and covariance matrix over a specified rectangular
%      region. In other words, it is the cumulative distribution function
%      (CDF) over a rectangular region.
%
%INPUTS: rectMin The 2X1 minimum bounds of the rectangular region.
%   rectMax The 2X1 maximum bounds of the rectangular region. Note that
%           rectMax>=rectMin.
%        mu The 2X1 mean of the distribution. If this is omitted or an
%           empty matrix is passed, then the default of [0;0] is used.
%         R The 2X2 covariance matrix of the distribution. If this is
%           omitted or an empty matrix is passed, then R=eye(2,2) is used.
%
%OUTPUTS: val The value of the integral.
%
%We can use the function bivarNormCDF to evaluate the probability of any
%region of the form Pr{x1<b(1), x2<b(2)). Consider 4 points of a rectangle
%dividing 2D space into 9 region.
%   a | b | c
% ----1---2---
%   d | e | f
% ----3---4---
%   g | h | i
%We would like to find the probability in region e. Well, if Pr{number} is
%the probability that everything is below and left of point specified by
%number (i.e. what one gets with the bivarNormCDF function), then one can
%write 4 equations in 4 unknowns:
% Pr{1}=Pr{d}+Pr{g}
% Pr{2}=Pr{d}+Pr{e}+Pr{g}+Pr{h}
% Pr{3}=Pr{g}
% Pr{4}=Pr{g}+Pr{h}
%where Pr{region name} is the probability of of that region. This is a
%system of 4 equation and 4 unknowns, and we solve it to find Pr{e}:
% Pr{e}=-Pr{1}+Pr{2}+Pr{3}-Pr{4}
%
%EXAMPLE 1:
%As a quick check, we know that the integral over the lower-left quadrant
%fof the bivarate normal distribution equals 0.25. That result is obtained
%here.
% mu=[0,;0];
% R=eye(2,2);
% rectMin=[-Inf;-Inf];
% rectMax=[0;0];
% P=bivarGaussRectangleCDF(rectMin,rectMax,mu,R)
%
%EXAMPLE 2:
%Here, we compare the result of this function with the lower-fidelity
%approximation from GaussianD.integralOverRegion. The values are in the
%general ballpark.
% mu=[1,;2];
% R=[3,-2;
%   -2, 4];
% rectMin=[-3;-2];
% rectMax=[2;1];
% P=bivarGaussRectangleCDF(rectMin,rectMax,mu,R)
% PApprox=GaussianD.integralOverRegion(mu,R,rectMin,rectMax)
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public releas

if(nargin<4||isempty(R))
    R=eye(2,2);
end

if(nargin<3||isempty(mu))
    mu=[0;0];
end

p1=[rectMin(1);rectMax(2)];
p2=[rectMax(1);rectMax(2)];
p3=[rectMin(1);rectMin(2)];
p4=[rectMax(1);rectMin(2)];

Pr1=bivarNormCDF(p1,mu,R);
Pr2=bivarNormCDF(p2,mu,R);
Pr3=bivarNormCDF(p3,mu,R);
Pr4=bivarNormCDF(p4,mu,R);

%The min and max is just to handle finite precision effects making the
%result slightly invalid.
probVal=min(1,max(0,-Pr1+Pr2+Pr3-Pr4));

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
