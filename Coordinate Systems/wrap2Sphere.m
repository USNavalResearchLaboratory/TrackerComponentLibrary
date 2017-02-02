function wrapVals=wrap2Sphere(azEl)
%%WRAP2SPHERE  Give az-muth-elevation directions in spherical coordinates
%              that might be outside of the standard range of -pi to pi for
%              azimuth and -pi/2 to pi/2 for evevation, wrap the points to
%              the sphere. This function is useful when offsets (like
%              simulated noise) are added to spherical coordinates. For
%              example, if elevation exceeds pi/2 by epsilon, then it must
%              be reduced to pi/2-epsilon and azimuth offset by pi.
%
%INPUTS: azEl  A 2XN set of N spherical coordinates in radians that should
%              be wrapped to the normal sphere. azEl(1,:) are azimuth
%              values measured counterclockwise from the x-axis. azEl(2,:)
%              are elevation values measured up from the equator.
%
%OUTPUTS: wrapVals Values of azEl that have been wrapped to the sphere in
%                  radians. azEl(1,:) is from -pi to pi and azEl(2,:) is
%                  from -pi/2 to pi/2.
%
%August 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numVals=size(azEl,2);

wrapVals=zeros(2,numVals);
for curVal=1:numVals
    lambda=azEl(1,curVal);
    theta=azEl(2,curVal);
    
    %Wrap the azimuth value, dependent on the elevation value.
    if(mod((abs(theta)-pi/2)/pi,2)>1)
        wrapVals(1,curVal)=wrapRange(lambda,-pi,pi);
        wrapVals(2,curVal)=asin(sin(theta));
    else
        wrapVals(1,curVal)=wrapRange(lambda+pi,-pi,pi);
        wrapVals(2,curVal)=asin(sin(theta));
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
