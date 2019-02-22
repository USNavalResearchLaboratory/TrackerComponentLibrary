function res=accurateDotProdK(x,y,K)
%%ACCURATEDOTPRODK Find the dot product of vectors x and y with
%             intermediate results compute at K times double precision.
%
%INPUTS: x,y Real, nX1 or 1XN vectors.
%          K The multiple of the standard double precision desired for
%          intermediate results. If omitted or an empty matrix is passed,
%          the default is K=2. Values larger than n are clipped to n. 
%
%OUTPUTS: res The scalar value fo the dot product.
%
%This implements algorithm 5.10 of [1].
%
%EXAMPLE:
%Here, we have an example, where some large values cancel. In the standard
%dot product, these values cause precision to be lost.
% x=[1e102;-159;11;78;100;9;-37;-148;-4;96];
% y=[100;-43;-162;16;-1e102;-22;-114;202;-235;-50];
% res0=accurateDotProdK(x,y)
% resMatlab=dot(x,y)
%One will get res0=-23433, which is correct. However, resMatlab is -29736,
%which is incorrect (as of Matlab 2018b).
%
%REFERENCES:
%[1] T. Ogita, S. M. Rump, and S. Oishi, "Accurate sum and dot product,"
%    SIAM Journal on Scientific Computing, no. 6, pp. 1955-1988, 2005.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(K))
    K=2;
end

n=length(x);

r=zeros(2*n,1);
[p,r(1)]=exactPairMult(x(1),y(1));
for i=2:n
    [h,r(i)]=exactPairMult(x(i),y(i));
    [p,r(n+i+1)]=exactPairSum(p,h);
end
r(2*n)=p;
res=accurateSumK(r,K);

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
