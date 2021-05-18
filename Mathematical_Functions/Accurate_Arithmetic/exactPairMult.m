function [z,e]=exactPairMult(xIn,yIn)
%%EXACTPAIRMULT Calculate the product of the floating point double values
%               x and y such that the result z+e (if added with infinite
%               precision) is exactly the product x*y, where z is the
%               floating point product and e is the part of the result that
%               cannot be expressed exactly. This is implemented for x and
%               y being floating point doubles.
%
%INPUTS xIn, yIn Two real arrays of sizes that can be added.
%
%OUTPUTS: z The real product x.*y as limited by floating point precision.
%         e The real error such that if done with inifite precision
%           (z+e)==x*y. This does not necessarily hold when numbers become
%           denormalized.
%
%This function implements Dekker's algorithm, which is given as Mul12 in
%[1]. Note that additional scaling is performed for the instance where
%Dekker's algorithm might result in an overflow. That is when x or y is
%above 2^(1023 (maximum exponent) -53 (mantissa length))=2^970.
%
%The output pair (z,e) can be considered to be a "doublelength" floating
%point number, as defined in [1]. This means that
%abs(e)<=abs(z+e)*2^(-t)/(1+2^(-t)) (considering exact arithmetic), where t
%is the number of bits in the mantissa of the floating point number. Here,
%that is 53.
%
%Note that this assumes that the processor rounding mode has been set to
%round to "nearest," which is the default in Matlab, but which can be
%changed using the function setProcRoundingMode.
%
%EXAMPLE:
%Consider the product
% [z,e]=exactPairMult(1+2^50,1+2^(-15))
%One will get z=1.125934266580993e+15 and e=3.051757812500000e-05.
%The exact value of the product can be found to be
%1125934266580993.000030517578125
%If one adds z and e exactly, one gets precisely that result.
%
%REFERENCES:
%[1] T. J. Dekker, "A Floating Point Technique for Extending the Available
%    Precision," Numerische Mathematik, vol. 18, no. 3, Jun. 1971, pp.
%    224-242.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(~isa(xIn,'double')||~isa(yIn,'double'))
    error('This function is only implemented for x and y being double floating point values.')
end

if(isempty(xIn)||isempty(yIn))
    %Returning empty matrices if one input is an empty matrix is consistent
    %with what Matlab does when trying to multiply an empty matrix with
    %something.
    z=[];
    e=[];
    return
end

%This is 2^(1023 (maximum exponent) -53)=2^970
maxUnscaled=9.9792015476735990582818635651842e291;
two53=9007199254740992;%2^53;
twom53=1.1102230246251565404236316680908e-16;%2^(-53)

isScalX=isscalar(xIn);
isScalY=isscalar(yIn);

if(xor(isScalX,isScalY))
    %If one is scalar and the other is not.
    if(isScalX)
        dims=size(yIn);
        xIn=repmat(xIn,dims);
    else
        dims=size(xIn);
        yIn=repmat(yIn,dims);
    end
else
    %If neither is scalar or both are scalar, assume that they are the same
    %size.
    dims=size(xIn);
end

z=zeros(dims);
e=zeros(dims);
numEl=numel(z);

for curEl=1:numEl
    x=xIn(curEl);
    y=yIn(curEl);
    if(x>maxUnscaled)
        x=x*twom53;
        invScalFactX=two53;
    else
        invScalFactX=1;
    end

    if(y>maxUnscaled)
        y=y*twom53;
        invScalFactY=two53;
    else
        invScalFactY=1;
    end
    invScalFact=invScalFactX*invScalFactY;

    [zCur,eCur]=exactPairMultPrescaled(x,y);

    z(curEl)=zCur*invScalFact;
    e(curEl)=eCur*invScalFact;
end
end

function [z,e]=exactPairMultPrescaled(x,y)
    %The double mantissa is 53 bits (including the implicit one).
    const=134217729;%This is 2^ceil(53/2)+1=2^27+1
    
    p=x*const;
    hx=(x-p)+p;
    tx=x-hx;
    p=y*const;
    hy=(y-p)+p;
    ty=y-hy;
    p=hx*hy;
    q=(hx*ty)+(tx*hy);
    z=p+q;
    e=p-z+q+tx*ty;

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
