function [s,e]=exactPairSum(a,b)
%%EXACTPAIRSUM This evaluates a+b such that s is the floating point value
%              given by a+b with fixed precision and e is a value such
%              that, with infinite precision, (s+e)=(a+b).
%
%INPUTS: a,b Two real arrays of sizes that can be added.
%
%OUTPUTS: s The real sum a+b as limited by floating point precision.
%         e The real error such that if done with infinite precision
%           (s+e)==(a+b). This does not necessarily hold when numbers
%           become denormalized.
%
%The basic notion is that when assing a+b, then the value with the smaller
%magnitude will have bits of precision tuncated off in the sum. Thus, if
%abs(b)<abs(a), then the error is e=b-(s-a), where the parentheses are
%important, and just the same swapping a and b if abs(a)>abs(b). This is
%described in [1].
%
%The output pair (s,e) can be considered to be a "doublelength" floating
%point number, as defined in [1]. This means that
%abs(e)<=abs(s+e)*2^(-t)/(1+2^(-t)) (considering exact arithmetic), where t
%is the number of bits in the mantissa of the floating point number. Here,
%that is 53 for double floating point values.
%
%EXAMPLE:
% Here we have a sum where it is easy to verify the result.
% [s,e]=exactPairSum(2e12,1+2e-15)
%Here, s=2.000000000001000e+12, which is 2e12+1. The extra 2e-15 is cut off
%and thus ends up in e, so e=1.998401444325282e-15.
%
%REFERENCES:
%[1] T. J. Dekker, "A Floating Point Technique for Extending the Available
%    Precision," Numerische Mathematik, vol. 18, no. 3, Jun. 1971, pp.
%    224-242.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

s=a+b;

aScal=isscalar(a);
bScal=isscalar(b);
if(aScal&&bScal)
    if(abs(a)>abs(b))
        e=b-(s-a);
    else
        e=a-(s-b);
    end
elseif(~aScal&&~bScal)
    e=zeros(size(s));
    numE=numel(e);
    
    for k=1:numE
        if(abs(a(k))>abs(b(k)))
            temp=s(k)-a(k);
            e(k)=b(k)-temp;
        else
            temp=s(k)-b(k);
            e(k)=a(k)-temp;
        end
    end
elseif(aScal)
    e=zeros(size(s));
    numE=numel(e);
    absA=abs(a);

    for k=1:numE
        if(absA>abs(b(k)))
            temp=s(k)-a;
            e(k)=b(k)-temp;
        else
            temp=s(k)-b(k);
            e(k)=a-temp;
        end
    end
else%bScal is true
    e=zeros(size(s));
    numE=numel(e);
    absB=abs(b);
    
    for k=1:numE
        if(abs(a(k))>absB)
            temp=s(k)-a(k);
            e(k)=b-temp;
        else
            temp=s(k)-b;
            e(k)=a(k)-temp;
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
