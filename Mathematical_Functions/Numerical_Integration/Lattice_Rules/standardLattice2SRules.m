function [xi,w,errParams]=standardLattice2SRules(numDim,methodIdx,periodizeAlg,upperBounds,lowerBounds)
%%STANDARD2SLATTICERULES Obtain lattice points for multidimensional
%               numerical integration based on the 2^s lattice rules that
%               are Given in Appendix A of [1], wherbey the method of
%               obtaining the points is described in Chapter 10 of [1].
%               These points and be used for typical numerical integration
%               like those from the KorobovLatticeRules function, or the
%               integrateNsLatticeRule function can be used to perform
%               integration AND obtain an error estimate. This is possible
%               because of the embedded nature of the points. The standard
%               points are meant for integrating a function that is
%               periodic in all dimensions with the period being on the
%               region 0-1 in all dimensions. Inputs to this function can
%               perform a transformation to periodize the points (for use
%               with non-periodic functions, though periodization is not
%               always necessary), and can change the integration region to
%               other finite or infinite bounds.
%
%INPUTS: numDim The number of dimensions of the lattice points that are to
%               be generated. Rules for 2<=numDim<=12 are currently
%               avaiable.
%     methodIdx The index of the formula to use. For each number of
%               dimensions, 7 formula are available. They are those taken
%               from Appendix A of [1]. methodIdx thus ranges from 1 to 7
%               with an increasing number corresponding to an increasing
%               number of lattice points used. If omitted or an empty
%               matrix is passed, then method 1 is used.
%  periodizeAlg If this parameter is omitted or an empty matrix is passed,
%               then no periodization of the points will be performed.
%               Periodization of the points can improve performance when
%               used for integration over non-periodic functions. If this
%               parameter is provided, then it is a sturcture with element
%               periodizeAlg.method and optionally periodizeAlg.param1 and
%               periodizeAlg.param2. The meaning of the parameters is the
%               same as in the function periodizeLatticePoints.
% lowerBounds,upperBounds If these parameters are provided (and are not
%               empty matrices), then the lattice points will be
%               transformed for integration over regions other than just
%               +/-1 in all dimensions. These are numDimX1 or 1XnumDim
%               vectors where lowerBounds are the lower integration bounds
%               and upperBounds are the upper integration bounds. If
%               transformed, the following substitutions are used:
%               Considering only the current dimensions of xi where a is
%               the lower bound and b is the upper bound:
%               * If a and b are finite:
%                 xi(curDim,:)=(b-a)*xi(curDim,:)+a, w=(b-a)*w
%               * If a=-Inf and b is finite:
%                 xi(curDim,:)=-1./xi(curDim,:)+b+1; w=w./xi(curDim,:).^2;
%                 --Note this may cause problems if any elements of xi are
%                 zero.
%               * If a is finite and b=Inf:
%                 xi(curDim,:)=1./(1-xi(curDim,:))+a-1;
%                 w=w./(1-xi(curDim,:)).^2;
%                 --Note this may cause problems if any elements of xi are
%                 1.
%               * If a=-Inf and b=Inf:
%                 xi(curDim,:)=1./(1-xi(curDim,:))-1./xi(curDim,:);
%                 w=w.*(1./(1-xi(curDim,:)).^2+1./xi(curDim,:).^2);
%                 --Note this may cause problems if any elements of xi are
%                 zero or 1.
%               * Same as any above, but the lower bound is greater than
%                 the upper bound:
%                 The bounds are flipped, but the same substitution is
%                 made, just w is multiplied by -1.
%
%OUTPUTS: xi A numDimXnumPoints matrix of the lattice points generated
%            using the selected method.
%          w A numPointsX1 set of weights associated with the lattice
%            points. If no periodization or transformation of the
%            integration region took place, then these are all just
%            1/numPoints.
%  errParams This is a structure containing parameters that can be passed
%            to the integrateNsLatticeRule function so that these points
%            can be used to perform numerical integration with an estimate
%            of the error of the integral.
%
%REFERENCES:
%[1] I. H. Sloan and S. Joe, Lattice Methods for Multiple Integration.
%    Oxford: Clarendon Press, 1994.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%These parameters are from tables A.1 to A.12 from Appendix A of [1]. The
%first column is m, the second l.
switch(numDim)
    case 2
        paramTable=[1249,   512;
                    2503,   672;
                    5003,   1850;
                    10007,  3822;
                    20011,  6103;
                    40009,  15152;
                    80021,  30954];
    case 3
        paramTable=[619,    178;
                    1249,   136;
                    2503,   652;
                    5003,   1476;
                    10007,  2325;
                    20011,  2894;
                    40009,  8789];
    case 4
        paramTable=[313,    51;
                    619,    73;
                    1249,   197;
                    2503,   792;
                    5003,   792;
                    10007,  1206;
                    20011,  8455];
    case 5
        paramTable=[157,    18;
                    313,    51;
                    619,    104;
                    1249,   165;
                    2503,   792;
                    5003,   380;
                    10007,  1927];
    case 6
        paramTable=[79,     27;
                    157,    18;
                    313,    80;
                    619,    102;
                    1249,   175;
                    2503,   253;
                    5003,   162;
                    10007,  2286];
    case 7
        paramTable=[41,     13;
                    79,     27;
                    157,    11;
                    313,    70;
                    619,    161;
                    1249,   303;
                    2503,   306;
                    5003,   2037;
                    10007,  2286];
    case 8
        paramTable=[19,     2;
                    41,     13;
                    79,     27;
                    157,    11;
                    313,    70;
                    619,    161;
                    1249,   155;
                    2503,   153;
                    5003,   137;
                    10007,  373];
    case 9
        paramTable=[11,     2;
                    19,     2;
                    41,     13;
                    79,     27;
                    157,    11;
                    313,    93;
                    619,    161;
                    1249,   18;
                    2503,   153;
                    5003,   657;
                    10007,  563];
    case 10
        paramTable=[5,      2;
                    11,     2;
                    19,     2;
                    41,     13;
                    79,     27;
                    157,    36;
                    313,    62;
                    619,    106;
                    1249,   184;
                    2503,   1019;
                    5003,   189;
                    10007,  1047];
    case 11
        paramTable=[3,      1;
                    5,      2;
                    11,     2;
                    19,     4;
                    41,     13;
                    79,     8;
                    157,    36;
                    313,    15;
                    619,    106;
                    1249,   113;
                    2503,   29;
                    5003,   189;
                    10007,  2499];
    case 12
       paramTable=[3,       1;
                   5,       2;
                   11,      3;
                   19,      3;
                   41,      13;
                   79,      8;
                   157,     36;
                   313,     15;
                   619,     80;
                   1249,    467;
                   2503,    29;
                   5003,    331;
                   10007,   745];
    otherwise
        error('Unsupported dimensionality chosen')
end

if(nargin<2||isempty(methodIdx))
    methodIdx=1;
end

m=paramTable(methodIdx,1);
l=paramTable(methodIdx,2);

z=mod(l.^(0:(numDim-1)),m)';

[xi,errParams]=convertNSLatticeRule2Points(z,m);

%If the points should be periodized.
if(nargin>2&&~isempty(periodizeAlg))
    method=periodizeAlg.method;
    if(isfield(periodizeAlg,'param1'))
        param1=periodizeAlg.param1;
    else
        param1=[];
    end
    if(isfield(periodizeAlg,'param2'))
        param1=periodizeAlg.param2;
    else
        param2=[];
    end
    
    [xi,w]=periodizeLatticePoints(method,xi,[],param1,param2);
else%Otherwise, the xi do not change and the weights are just uniform.
    numPoints=size(xi,2);
    w=(1/numPoints)*ones(numPoints,1);
end

%If lower and upper bounds are provided.
if(nargin>3&&~isempty(lowerBounds))
    for curDim=1:numDim
        a=lowerBounds(curDim);
        b=upperBounds(curDim);
        
        if(a<b)
            flipSign=1;
        else%Make a<b and flip the sign of the result.
            temp=a;
            a=b;
            b=temp;
            flipSign=-1;
        end
        %Now, a<b.

        if(isfinite(a)&&isfinite(b))%Both bounds finite
            xi(curDim,:)=a+(b-a)*xi(curDim,:);
            w=flipSign*w*(b-a);
        elseif(~isfinite(a)&&isfinite(b))%a=-Inf
            xi(curDim,:)=-1./xi(curDim,:)+b+1;
            w=flipSign*w./xi(curDim,:).^2;
        elseif(isfinite(a)&&~isfinite(b))%b=Inf
            xi(curDim,:)=1./(1-xi(curDim,:))+a-1;
            w=flipSign*w./(1-xi(curDim,:)).^2;
        else%a=-Inf and b=Inf
            xi(curDim,:)=1./(1-xi(curDim,:))-1./xi(curDim,:);
            w=flipSign*w.*(1./(1-xi(curDim,:)).^2+1./xi(curDim,:).^2);
        end
    end
end
end

function [xi,errParams]=convertNSLatticeRule2Points(z,m,n)
%%CONVERTNSLATTICERULES2POINTS This algorithm is essentially an
%              implementation of the algorithm in Chapter 10.4 of [1] that
%              has been modified to collect the points rather than feed
%              them into the function and evaluate Q and Q^(i). To be able
%              to give the points to Q and Q^(i) at a later time to
%              evaluate an integral and estimate its accuracy, the
%              integrateNsLatticeRule function can be used.
%
%REFERENCES:
%[1] I. H. Sloan and S. Joe, Lattice Methods for Multiple integration.
%    Oxford: Clarendon Press, 1994.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Assume a 2^s rule if not specified.
if(nargin<3||isempty(n))
    n=2;
end

s=size(z,1);
numPoints=2^s*m;

xi=zeros(s,numPoints);

errParams.omegaStart=zeros(numPoints,1);
errParams.includeQParenth=false(s,numPoints);

k=zeros(s,1);
k(1)=-1;%Summation variables in Equation 10.1, pg. 165.
l=1;
%Omega labels the outermost sum in 10.1 up to which a summation index
%has eveer been incremented.
omega=0;

curPoint=1;
while(l<=s)
   if(k(l)<n-1)
       k(l)=k(l)+1;
       l=1;
       for j=0:(m-1)
           xi(:,curPoint)=mod((j/m)*z+k/n,1);
           errParams.omegaStart(curPoint)=omega;
           errParams.includeQParenth(k==0,curPoint)=true;
           curPoint=curPoint+1;
       end

       if(omega==0)
           omega=1;
       end
   else
      k(l)=0;
      l=l+1;
      if(omega<l)
          omega=l;
      end
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
