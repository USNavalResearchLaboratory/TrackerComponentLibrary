function y=evalSinCosSeries(c,cosX,sinX)
%%EVALSINCOSSERIES Evaluate the cosine series
%                  y=sum_{i=1}^Nc(i)*cos((i-1)*x)
%                  or the sine series
%                  y=sum_{i=1}^(N+1)c(i)*sin(i*x)
%                  efficiently using Clenshaw summation.
%
%INPUTS: c An NX1 or a 1XN vector of the coefficients in the sum.
%     cosX The value cos(x). This value is always required. To evaluate
%          multiple sums at once, this can be a vector or a matrix.
%     sinX If this value is provided, it is sin(x) and the sum evaluated by
%          this function will be a sine series. Otherwise, if this
%          parameter is omitted or an empty matrix is passed, this function
%          will evaluate a cosine series. If this parameter is provided, it
%          must be the same size as cosX.
%
%OUTPUTS: y The value of the sum. If cosX was a matrix, then this will be a
%           matrix.
%
%This function chooses the correct parameters and calls the function
%evalClenshawRecurSeries.
%
%EXAMPLE:
% N=200;
% c=rand(N,1);
% x=2*pi*rand(4,4);
% cosX=cos(x);
% yF=evalSinCosSeries(c,cosX);
% ySum=zeros(4,4);
% for k=1:N
%    ySum=ySum+c(k)*cos((k-1)*x); 
% end
% relativeError=max(max(abs(ySum-yF)./yF))
% %This will typically be on the order of 1e-12 or less, due simply to
% %finite precision differences.
% %Similarly, if we wanted the sine series.
% sinX=sin(x);
% yF=evalSinCosSeries(c,cosX,sinX);
% ySum=zeros(4,4);
% for k=1:N
%     ySum=ySum+c(k)*sin(k*x); 
% end
% relativeError=max(max(abs(ySum-yF)./yF))
% %This will again be on the order of 1e-12 or less, due simply to finite
% %precision differences.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

alphaVal=2*cosX;
betaVal=-1;

if(nargin<3||isempty(sinX))
    %A cosine series.
    F0=ones(size(cosX));%=cos(0)
    F1=cosX;
    y=evalClenshawRecurSeries(c,alphaVal,betaVal,F0,F1);
else%A sine series
    F0=sinX;
    F1=2*sinX.*cosX;
    y=evalClenshawRecurSeries(c,alphaVal,betaVal,F0,F1);
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
