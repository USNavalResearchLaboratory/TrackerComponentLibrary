function evalPts=evalCubSpline(x,C,tau,numDerivs)
%%EVALCUBSPLINE Given the output of the fitCubSpline function, evaluate the
%               value of the fitted spline and also a number of derivatives
%               if desired.
%
%INPUTS: x A length numPts set of points at which interpolation should be
%          performed.
%        C The 4XN matrix of coefficients returned by the fitCubSpline
%          function.
%      tau The 1XN or NX1 vector of independent variable points fitted in
%          the fitCubSpline function that produced C (breakpoints).
%       pp The structure returned by the fitCubSpline function, which has
%          members tau (the array of breakpoints) and C (the matrix of
%          coefficients).
% numDerivs A number from 0 to 3 specifying the number of derivatives
%          desired. The default if omitted or an empty matrix is passed is
%          0.
%
%OUTPUTS: evalPts A (numDerivs+1)X1 matrix of the interpolated points and
%                 derivatives.
%
%The interpolation equation given the coefficients is given in Chapter 4 of
%[1]. Having determined that a point x(i) lies between knots tau(k) and
%tau(k+1) (or any point after tau(N-1) when dealing with a set of N knots),
%the value is interpolated as
%val(i)=C(1,k)+C(2,k)*(x(i)-tau(k))+C(3,k)*(x(i)-tau(k))^2/2+C(4,k)*(x(i)-tau(k))^3/6
%The derivatives just come from the derivative of that equation with
%respect to x(i).
%
%REFERENCES:
%[1] C. de Boor, A Practical Guide to Splines. New York: Springer-Verlag,
%    1978.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(numDerivs))
    numDerivs=0;
end

numTau=length(tau);
N=length(x);

evalPts=zeros(numDerivs+1,N);

[~,idxPrev]=binSearch(tau,x(1),1);
idxPrev=idxPrev-(idxPrev==numTau);

for k=1:N
    if((x(k)>=tau(idxPrev)&&x(k)<=tau(idxPrev+1)))
        %No need to search for the correct region.
        idxCur=idxPrev;
    elseif(idxPrev<numTau&&(x(k)>=tau(idxPrev+1)&&x(k)<=tau(idxPrev+2)))
        %If it is just the next region. This helps avoid calling binsearch
        %when given a large number of sorted points to interpolate.
        idxCur=idxPrev+1;
    else
        [~,idxCur]=binSearch(tau,x(k),1);
      	idxCur=idxCur-(idxCur==numTau);
    end
    
    diff=x(k)-tau(idxCur);
    diff2=diff*diff;
    diff3=diff*diff2;
    evalPts(1,k)=C(1,idxCur)+C(2,idxCur)*diff+C(3,idxCur)*diff2/2+C(4,idxCur)*diff3/6;
    
    if(numDerivs>0)
        evalPts(2,k)=C(2,idxCur)+C(3,idxCur)*diff+C(4,idxCur)*diff2/2;
        
        if(numDerivs>1)
            evalPts(3,k)=C(3,idxCur)+C(4,idxCur)*diff;
            
            if(numDerivs>2)
                evalPts(4,k)=C(4,idxCur);
            end
        end
    end

    idxPrev=idxCur;
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
