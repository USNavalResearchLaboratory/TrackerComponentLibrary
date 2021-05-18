function [C,S]=fitSpherHarmonics(M,points,vals,a,c)
%%FITSPHERHARMONICS Given a number of field values at a diversity
%               of points in spherical coordinates (either range, azimuth
%               and elevation, or just azimuth and elevation), fit fully
%               normalized spherical harmonic coefficients up to a
%               particular degree to the data. The data being fitted can be
%               real or complex. One can fit multiple sets of data at once
%               as long as they share the same set of points.
%
%INPUTS: M The maximum degree of the coefficients to use.
%   points The 2XN or 3XN set of points in spherical coordinates in the
%          form [r;azimuth;elevation] or just [azimuth;elevation]. Azimuth
%          is measured counterclockwise from the x-axis in the x-y plane.
%          Elevation is measured up from the x-y plane (towards the
%          z-axis). Angles must be given in radians. Note the system to be
%          observable, N>=(M+1)^2.
%     vals The NXnumSets sets of real or complex values at the given
%          points.
%        a The real or complex numerator in the (a/r)^n term in the
%          spherical harmonic sum. If this parameter is omitted or an empty
%          matrix is passed, then a=1 is used.
%        c The real or complex constant value by which the spherical
%          harmonic series is multiplied. If this parameter is omitted or
%          an empty matrix is passed, then c=1 is used.
%
%OUTPUTS: C A numCoeffsXnumSets collections of arrays, one for each set of
%           values passed, holding the coefficient terms that are
%           multiplied by cosines in the spherical harmonic expansion. If
%           given to a CountingClusterSet class, C(n+1,m+1) is the
%           coefficient of degree n and order m. This and S can be given to
%           the function spherHarmonicEval.
%         S A numCoeffsXnumSets collection of arrays, one for each set of
%           values passed, holding the coefficient terms that are
%           multiplied by sines in the spherical harmonic expansion. S has
%           the same structure as C. 
%
%The spherical harmonic series being fit is assumed to be of the form
%V=(c/r)*sum_{n=0}^M\sum_{m=0}^n(a/r)^n*(C(n+1,m+1)*cos(m*lambda)+S(n+1,m+1)*sin(m*lambda))*\bar{P}_{nm}(sin(theta))
%where \bar{P}_{nm} is a fully-normalized associated Legendre polynomial.
%This function uses the LegendreCos function to evaluate the Legendre
%polynomials and it them constructs a linear system to solve for the values
%of C and S. Unobservable values of S are where m=0, so S(n+1,0+1)=0 and
%those components are not part of the estimation problem.
%
%Real spherical harmonic models are often used for gravitational, magnetic
%field and some terrain modeling, as discussed in part of [1], whereas
%complex spherical harmonic models are often used for electromagnetic
%modeling.
%
%The complex spherical harmonic model that is usually used in
%electromagnetics is 
%V=(c/r)*sum_{n=0}^M\sum_{m=-n}^n(a/r)^n*A(n+1,m+1)*\bar{P}_{nm}(sin(theta))*exp(1j*m*lambda)
%Initially, it just looks like one has combined the S and C values into one
%complex coefficient A and then multiplied the sin term by 1j. However, if
%that were all that had been done, one wouldn't be able to match many
%models, because if A(n+1,m+1)=AR+1j*AI, expanding the term
%A(n+1,m+1)*exp(1j*m*lambda), one gets 
%aR*cos(m*lambda)-aI*sin(m*lambda)+1j*(aI*cos(m*lambda)+aR*sin(m*lambda))
%Where one can see that it is not possible to independently set the real
%and the imaginary terms. However, note that the formulation of the inner
%sum used for complex values starts at -m and not at 0. This extension of
%the sum allows one to uniquely form real and complex terms. Though such a
%symmetry does not occur at m=0, that does not actually matter, since the
%sine terms are zero, so one can still independently set the real and the
%imaginary coefficients parts. Indeed, it can be shown that the resulting
%expression allows one to fit the real and complex parts independently if a
%and c are real.
%
%EXAMPLE:
%We take the coefficients for the World Magnetic Model, evaluate the
%magnetic field values at a number of points around the globe, use this
%function to reconstruct the spherical harmonic coefficients from the field
%values and compare the computed coefficients with the coefficients used to
%generate the field values.
% [C,S,a,c]=getWMMCoeffs();
% M=(1/2)*(sqrt(1+8*length(C))-1);
% 
% %Evaluate the model on a grid.
% numPoints=21;
% totalGridPoints=numPoints*numPoints;
% deltaLat=180/numPoints;
% lat=-90+(0:(numPoints-1))*deltaLat;
% lat=lat*(pi/180);%Convert latitude to radians.
% 
% deltaLon=360/numPoints;
% lon=-180+(0:(numPoints-1))*deltaLon;
% lon=lon*(pi/180);%Convert longitude to radians.
% [latGrid,lonGrid]=ndgrid(lat,lon);
% latLonEllipse=[latGrid(:)';lonGrid(:)'];
% %Convert from ellipsoidal latitudes to spherical latitudes (This assumes
% %that the points are on the surface of the reference ellipsoid).
% points=ellips2Sphere([latLonEllipse;zeros(1,totalGridPoints)]);
% 
% %Evaluate the magnetic field values at the given points.
% V=spherHarmonicEval(C,S,points,a,c);
% 
% %Reconstruct the coefficients from the magnetic field values.
% [CFit,SFit]=fitSpherHarmonics(M,points,V,a,c);
% 
% %Check how good the fit is to the original set of coefficients.
% sel1=C(:)~=0;
% max(abs(CFit(sel1)-C(sel1))./C(sel1))
% sel1=S(:)~=0;
% max(abs(SFit(sel1)-S(sel1))./S(sel1))
%In both instances, the maximum relative error should have been less than
%8e-10. The fact that there is any error stems from finite precision
%effects.
%
%REFERENCES:
%[1] D. F. Crouse, "An overview of major terrestrial, celestial, and
%    temporal coordinate systems for target tracking," Naval Research
%    Laboratory, Washington, DC, Tech. Rep. NRL/FR/5344-16-10,279, 10 Aug.
%    2016.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(c))
    c=1; 
end

if(nargin<4||isempty(a))
    a=1; 
end

pointDims=size(points,1);

numEq=size(points,2);

%The total number of associated Legendre polynomials.
numPBar=(M+1)*(M+2)/2;

%The M+1 being subtracted is the number of S_{n,0} values,
%because these can just be set to zero as they are never used; they are
%always multiplied by sin(0).
numCCoeffs=numPBar;
numSCoeffs=numPBar-(M+1);

%The total number of coefficients if all values for C and S must be
%computed. 
coeffsTotal=numCCoeffs+numSCoeffs;

if(coeffsTotal>numEq)
   warning('The value of M is too large for the number of values provided. The exact solution is unobservable.') 
end

%The columns of Y corresponds to the C coefficients followed by the S
%coefficients.
Y=zeros(numEq,coeffsTotal);
for curEq=1:numEq
    curPoint=points(:,curEq);
    
    if(pointDims==3)
        r=curPoint(1);
        Az=curPoint(2);
        coLat=pi/2-curPoint(3);
    else
        r=1;
        Az=curPoint(1);
        coLat=pi/2-curPoint(2);
    end
    %We now have r, azimuth, and colatitude.

    %Get all of the associated Legendre polynomials.
    PBarVals=LegendreCos(coLat,M);
    PBarVector=PBarVals(:);
    
    outerCoeff=c/r;
    innerCoeff=1;
    for n=0:M
        for m=0:n
            %The index of PBarVals(n+1,m+1) in the vector PBarVector.
            idx=n*(n+1)/2+m+1;
            
            coeffVal=outerCoeff*innerCoeff*PBarVector(idx);

            Y(curEq,idx)=coeffVal*cos(m*Az);
            if(m>0)
                %The -(n+1) deals with the fact that there are no values of
                %S for m=0.
                Y(curEq,numCCoeffs+idx-(n+1))=coeffVal*sin(m*Az);
            end
        end
        innerCoeff=innerCoeff*(a/r);
    end
end

%If one were to use linsolve(Y,vals(:)), many systems would complain about
%Y having poor conditioning.
coeffs=lsqminnorm(Y,vals);

C=coeffs(1:numPBar,:);

%Fill in S, noting that the unused parts are left set to zero.
numSets=size(vals,2);
SData=zeros(numPBar,numSets);
curSIdx=numCCoeffs+1;
for n=1:M
    for m=1:n
        idx=n*(n+1)/2+m+1;
        SData(idx,:)=coeffs(curSIdx,:);
        curSIdx=curSIdx+1;
    end
end
S=SData;

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
