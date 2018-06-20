function [a,c]=convertPolynomialForm(formOrig,formDes,a,cIn,cNew)
%%CONVERTPOLYNOMIALFORM Convert a polynomial given as a power series, in
%                       Newton's form, or in Taylor form to a different
%                       form. Note that finite precision limitations can
%                       cause problems when given coefficients over a large
%                       range.
%
%INPUTS: formOrig,formDes Strings specifying the original and desired forms
%                 of the polynomials. Possible values and the corresponding
%                 power series in terms of coefficients a and possible
%                 point or points c are
%                 'PowerSeries' for
%                        y(z)=a(end)+a(end-1)*z+a(end-2)*z^2+...
%                        (The same as Matlab's polyval function uses).
%                 'Taylor' 
%                        y(z)=a(1)+a(2)*(z-c)+a(3)*(z-c)^2+...
%                 'Newton'
%                        y(z)=a(1)+a(2)*(z-c(1))+a(3)*(z-c(2))^2+...
%               a The nX1 original polynomial coefficient vector.
%             cIn If the input series is 'Taylor', then this is the scalar
%                 point of the expansion. If the input series is 'Newton',
%                 then this is an (n-1)X1 vector of points. If the input
%                 series is 'PowerSeries', then this parameter is not
%                 required and an empty matrix can be passed.
%            cNew If the output is 'Newton' or 'Taylor', then this is
%                 either the (n-1)X1 vector of coefficients or the scalar
%                 point of expansion.
%
%OUTPUTS: a The coefficients in the new format given by formDes.
%         c The points associated with the new format. If formDes is
%           'PowerSeries', then c is the empty matrix, because no
%           coefficients are associated with a power series. If formDes is
%           'Newton', then c is the set of points associated with the
%           Newton form of the polynomial. If formDes is 'Taylor, then c is
%           a vector of repeated values of cNew, because the Taylor form is
%           the same as a Newton form where all of the points are the same.
%
%The conversions come from Chapter 19 of [1].
%
%REFERENCES:
%[1] A. Nijenhuis and H. S. Wilf, Combinatorial Algorithms for Computers
%    and Calculators, 2nd ed. New York: Academic press, 1978.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The number of coefficients.
n=length(a);

%Using the switch statement to change formOrig to numbers checks the
%validity of the input and avoids the use of lots of strcmp for further
%comparisons.
switch(formOrig)
    case 'PowerSeries'
        formOrig=0;
        %Convert to the form y(z)=a(1)+a(2)*z+a(3)*z^2+...
        a=flipud(a(:));
    case 'Taylor'
        %The Taylor series is the same as Newton's form with all of the
        %points being the same.
        formOrig=1;
        cIn=repmat(cIn(1),[n-1,1]);
    case 'Newton'
        formOrig=1;
    otherwise
        error('Invalid original form given.')
end

switch(formDes)
    case 'PowerSeries'
        formDes=0;
    case 'Taylor'
        %The Taylor series is the same as Newton's form with all of the
        %points being the same.
        formDes=1;
        cNew=repmat(cNew(1),[n-1,1]);
    case 'Newton'
        formDes=1;
    otherwise
        error('Invalid destination form given.')
end

%Power Series -> Power Series Do not do anything.
if(formOrig==0&&formDes==0)
    c=[];
    %Convert a back to the original ordering of
    %y(z)=a(end)+a(end-1)*z+a(end-2)*z^2+...
    a=flipud(a(:));
    return;
end

%Power Series -> Newton
if(formOrig==0&&formDes==1)
    for m=1:(n-1)    
        for i=(n-1):-1:m
            a(i)=a(i)+cNew(m)*a(i+1);
        end
    end
    c=cNew;
    return;
end

%Newton/ Taylor -> Power Series
if(formOrig==1&& formDes==0)
    for m=1:(n-1)
        for i=(n-1):-1:1
            a(i)=a(i)+(0-cIn(i))*a(i+1);
        end
        cIn=circshift(cIn(:),[1,0]);
        cIn(1)=0;
    end
    c=[];
    %Convert to the ordering of y(z)=a(end)+a(end-1)*z+a(end-2)*z^2+...
    %as used in the polyval function.
    a=flipud(a(:));
    return;
end

%If we are here, then the conversion is Newton/ Taylor -> Newton/ Taylor
%and we can assume that the basis points are probably changing.
for m=1:(n-1)
    cCur=cNew(end-m+1);
    for i=(n-1):-1:1
        a(i)=a(i)+(cCur-cIn(i))*a(i+1);
    end
    cIn=circshift(cIn(:),[1,0]);
    cIn(1)=cCur;
end
c=cNew;

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
