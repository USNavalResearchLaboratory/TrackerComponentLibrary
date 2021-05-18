function W=powerSeriesReversion(V)
%%POWERSERIESREVERSION Given a power series of the form 
%           z=V(1)*t+V(2)*t^2+V(3)*t^3... invert the series to get a series
%           of the form t=W(1)*z+W(2)*z^2+W(3)*z^3+... The resulting series
%           will be correct in a region around the origin. It cannot,
%           however, be correct everywhere, even though a strict inversion
%           might have been performed, because polynomials are usually not
%           bijective.
%                      
%INPUTS: V A 1XN or NX1 set of coefficients of the original power series.
%          They are in the order z=V(1)*t+V(2)*t^2+V(3)*t^3... Zero padding
%          the higher powers will result in a higher-order result and can
%          improve the fit of the inverted polynomial. The ordering of the
%          coefficients differs from that used in polyval and it lacks a
%          constant term.
%
%OUTPUTS: W The 1XN power series of the inverse series.
%
%This implements Algorithm L (Lagrangian power series reversion) in Chapter
%4.7 of [1].
%
%EXAMPLE:
%We try to invert z=t-t^2. We wil plot the original value as well as the
%value of the inverted series to compare. It will be seen that the inverted
%series is valid over a certain region
% numPoints=5000;
% t=linspace(-1,1,numPoints);
% z=t-t.^2;
% figure(1)
% clf
% hold on
% plot(t,z)
% 
% z=linspace(min(z),max(z),numPoints);
% %The inverse polynomial order, which is achieved by zero-padding.
% numCoeffs=28;
% W=powerSeriesReversion([1;-1;zeros(numCoeffs-2,1)]);
% 
% tVals=zeros(1,numPoints);
% for curPoint=1:numPoints
%     tVals(curPoint)=sum(W.*z(curPoint).^(1:numCoeffs));
% end
% plot(tVals,z)
% axis([min(t) max(t) -2 0.5])
% %One sees a good fit in a region near zero, and outside of that region, one
% %sees that the result is bad. nonetheless, the results can be inverted to
% %get back the original polynomial.
% powerSeriesReversion(W)
%
%REFERENCES:
%[1] D. Knuth, The Art of Computer Programming: Seminumerical Algorithms,
%    2nd ed. Reading, MA: Addison-Wesley, 1998, vol. 2.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Normalize to make the first coefficient 1 as in the book. After the
%conversion, V1 will have to be used to adjust the final output
%coefficients.
V1=V(1);
V=V/V1;

N=length(V);

W=zeros(1,N);
U=zeros(1,N);

V=V(:).';

W(1)=1;

%Step L1
U(0+1)=1;
for n=2:N
    %Step L3
    for k=1:(n-2)
        U(k+1)=U(k+1)-sum(U(((k-1):-1:0)+1).*V(2:(k+1)));
    end
    U(n-1+1)=-sum((2:n).*U(((n-2):-1:0)+1).*V(2:n));
    
    %Step L4
    W(n)=U(n-1+1)/n;
end

%Now, put back in the effects of the first coefficient
W=W./V1.^(1:N);

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
