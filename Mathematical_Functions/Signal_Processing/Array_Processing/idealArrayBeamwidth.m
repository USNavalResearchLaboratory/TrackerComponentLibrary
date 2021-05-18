function [BW,HBWL,HBWR,Rsp,U]=idealArrayBeamwidth(T,xyPoints,BWType,lineParams,bounds,numPoints)
%%IDEALARRAYBEAMWIDTH Determine the beamwidth of a tapered linear or
%               planar array that is composed of ideal isotropic antenna
%               elements. As there is no exact expression for generic
%               tapered arrays, this function find the beamwidth via a
%               brute-force search. The results are given in terms of
%               direction cosines, because in an ideal array, such
%               beamwidth are idependent of the steering direction, unlike
%               when the results are given in terms of angles.          
%
%INPUTS: T The numSubarraysXnumElements tapering matrix of the array. The
%          elements can be complex. If an empty matrix is passed, then it
%          is assumed that no tapering is used so T will be an identity
%          matrix. numElements>1.
% xyPoints A 1XNumDim or 2XnumDim matrix of the [x;y] locations of the
%          points in the linear or planar array. The units of the distances
%          are in terms of wavelengths for the narrowband model.
%   BWType The type of beamwidth to compute. Possible options are
%          'HalfPower' The beamwidth is determined by the distance from the
%                      peak to the half-power point on either side of the
%                      main beam. This is the default if this parameter is
%                      omitted or an empty matrix is passed. 
%          'NullToNull' This is the beamwidth from one null to another on
%                      either side of the main beam. The algorithm
%                      evaluates the  response on a gird. Rather than
%                      directly finding the null, this just finds the
%                      points on either side of the mean beam after which
%                      the response value increases. This will not work
%                      with two-element arrays.
% lineParams This is an optional parameter. The beamwidth in u might not be
%          the same as the beamwidth in v. These optional parameter
%          specifies which cut of the u-v plane (direction cosines, not
%          angles) to take when determining the beamwidth. This is only
%          used for 2D arrays. The default if omitted or an empty matrix is
%          passed is to take a cut along the U axis. The equation for the
%          line along with the beam pattern will be evaluated is
%          v=lineParams.intercept+lineParams.slope*u if the value
%          lineParams.vIndep=false or the vIndep component of
%          lineParams is omitted. Otherwise, the line is
%          u=lineParams.intercept+lineParams.slope*v .
%   bounds A 2X1 or 1X2 vector specifying the bounds of the independent
%          variable to use for the line search for the peak and beamwidth.
%          bounds(1) is the minimum; bounds(2) is the maximum. If this
%          parameter is omitted or an empty matrix is passed, the default
%          of bounds=[-1;1] is used.
% numPoints The function just evaluates the array response on a grid. This
%           is the number of points to use for the grid. If this parameter
%           is omitted or an empty matrix is passed, the default of 10000
%           is used. The grid is on the idependent parameter. If bounds do
%           not sufficiently limit the valid independent parameter values
%           for a particular set of lineParams, then fewer than numPoints
%           might be generated. The limit arises because values of
%           u^2+v^2>1 are not valid points.
%
%OUTPUTS: BW The full beamwidth. If either the left or right halfbeamwidths
%            could not be found, then this will be an empty matrix.
%       HBWL The distance from the peak to the left-side (lower independent
%            parameter values) of the beam until the beamwidth condition
%            in BWType is satisfied.
%       HBWR The distance from the peak to the right-side (higher 
%            independent parameter values) of the beam until the beamwidth
%            condition in BWType is satisfied.
%        Rsp The cut of the beampattern that was evaluated to find the
%            beamwidth. This corresponds to the absolute value of the
%            'RawOutput' option in the standardUVBeamPattern function.
%          U The values of the independent parameter at which the response
%            in Rsp is evaluated.
%
%Chapter 2.4.1 of [1] defines the beamwidth of an array.
%
%This function uses standardUVBeamPattern to evaluate the tapered array
%beam pattern on a linear grid of points. The peak is found and then the
%array is scanned in either direction until the beamwidth condition is
%satisfied. For both BWType values, one point after the required
%condition is satisfied is required. Using too coarse a grid or
%inappropriate bounds can lead to incorrect results.
%
%EXAMPLE 1:
%We consider a uniform linear array consisting of N elements space a
%half wavelenegth apart. All elements are uniformly weighted 1/N. In
%Chapter 2.4 of [1], it is shown that the array response is real and as a
%function of the direction cosine u is B=(1/N)*(sin(pi*u*N/2)/sin(pi*u/2))
%over all valid u values (-1<=u<=1). From that, in Section 2.4.1, the
%halfpower beamwidth is shown to be approximately 0.891*(2/N) and the
%null-to-null beamwidth (also known as the Rayleigh resolution) is 4/N.
%Here, we verify that this function gives the correct values.
% xPoints=0:(1/2):10;
% N=length(xPoints);
% T=(1/N)*eye(N,N);
% [BW,HBWL,HBWR]=idealArrayBeamwidth(T,xPoints,'HalfPower')
% BWApprox=0.891*(2/N)
% [BWNN,HBWLNN,HBWRNN]=idealArrayBeamwidth(T,xPoints,'NullToNull')
% BWNNExact=4/N
%One will see that the values are close to where they should be. Also, the
%two half-beamwidths are essentially the same as they 
%
%EXAMPLE 2:
%Here we use a circular planar array and rather than producing the width
%of a sum beam, we show how a difference beamwidth can be determined. This
%is the width of one  half of the difference beam. 
% xyPoints=getShaped2DLattice([30;30],'circular');
% sidelobedB=-30;
% N=17;
% T=BaylissTapering(sidelobedB,N,xyPoints).';
% [BW,HBWL,HBWR]=idealArrayBeamwidth(T,xyPoints)
%
%REFERENCES:
%[1] H. L. Van Trees, Optimum Array Processing. New York: Wiley-
%    Interscience, 2002.
%
%September 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(xyPoints,1);

if(nargin<3||isempty(BWType))
    BWType='HalfPower';
end

if((nargin<4||isempty(lineParams))&&numDim~=1)
    lineParams=[];
    lineParams.intercept=0;
    lineParams.slope=0;
    lineParams.vIndep=false;
else
    lineParams=[];
end

if(nargin<5||isempty(bounds))
    bounds=[-1;1];
end

if(nargin<6||isempty(numPoints))
    numPoints=10000;
end

[Rsp,U,V]=standardUVBeamPattern(T,xyPoints,'RawOutput',[],lineParams,bounds,numPoints);

numPoints=length(Rsp);%Actaul number of points generated.
Rsp=abs(Rsp);

if(~isempty(lineParams)&&lineParams.vIndep)
    U=V; 
end

%Find the maximum point.
[maxVal,maxIdx]=max(Rsp);

%We will then just scan element by element until we hit the desired value.
switch(BWType)
    case 'HalfPower'
        halfPowVal=maxVal/sqrt(2);

        %Scan to find the right half-beamwidth.
        curIdx=maxIdx+1;
        HBWR=[];
        while(curIdx<=numPoints)
            curVal=Rsp(curIdx);
            
            if(curVal<=halfPowVal)
                %Use linear interpolation to find the exact point. 
                slope=(curVal-Rsp(curIdx-1))/(U(curIdx)-U(curIdx-1));
                yIntercept=curVal-slope*U(curIdx);
                
                UHBW=(halfPowVal-yIntercept)/slope;
                HBWR=UHBW-U(maxIdx);
                break;
            end
            curIdx=curIdx+1;
        end
        
        %Scan to find the left half-beamwidth.
        curIdx=maxIdx-1;
        HBWL=[];
        while(curIdx>=1)
            curVal=Rsp(curIdx);
            
            if(curVal<=halfPowVal)
                %Use linear interpolation to find the exact point. 
                slope=(curVal-Rsp(curIdx+1))/(U(curIdx)-U(curIdx+1));
                yIntercept=curVal-slope*U(curIdx);
                
                UHBW=(halfPowVal-yIntercept)/slope;
                HBWL=U(maxIdx)-UHBW;
                break;
            end
            
            curIdx=curIdx-1;
        end

        BW=HBWL+HBWR;
    case 'NullToNull'
        %Scan to find the right half-beamwidth.
        curIdx=maxIdx+1;
        HBWR=[];
        valPrev=maxVal;
        while(curIdx<=numPoints)
            curVal=Rsp(curIdx);
            
            if(curVal>valPrev)
                HBWR=U(curIdx-1)-U(maxIdx);
                break;
            end
            
            valPrev=curVal;
            curIdx=curIdx+1;
        end
        
        curIdx=maxIdx-1;
        HBWL=[];
        valPrev=maxVal;
        while(curIdx>=1)
            curVal=Rsp(curIdx);
            
            if(curVal>valPrev)
                HBWL=U(maxIdx)-U(curIdx+1);
                break;
            end
            
            valPrev=curVal;
            curIdx=curIdx-1;
        end
        
        BW=HBWL+HBWR;
    otherwise
        error('Unknown beamwidth type specified')
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
