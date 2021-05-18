function vals=Debye(x,n,kind,ulim)
%%DEBYE A function to compute the Debye functions at a point or vector of
%       points. This computes the integral only. Some people refer to the
%       Debye function as the integral part times n/x^n.
%
%INPUTS:
% x: A vector of values at which the Debye function will be evaluated.
% n: The scalar order of the Debye function. Defaults to 1.
% kind: 1 or 2 specifying first or second kind Debye functions. Defaults to
%       1.
% ulim: The scalar upper limit for evaluating a second kind Debye function.
%       Defaults to 1e6.
%
%OUTPUTS:
% vals: A vector of the same size as x containing the evaluated values.
%
%The Debye function is defined in [1]. Note that the Debye function
%technically is undefined at x=0, but this function defines Debye(0) as 0
%to enforce continuity.
%
%EXAMPLE 1: Plot a few examples of the function.
% x = linspace(-1,10);
% figure()
% hold on
% for n = 1:5
% y = Debye(x,n);
% plot(x,y)
% end
%
%REFERENCES:
%[1] Weisstein, Eric W. "Debye Functions." From MathWorld--A Wolfram Web
%    Resource. https://mathworld.wolfram.com/DebyeFunctions.html
%
%October 2020 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(n))
    n=1;
end

if(nargin<3||isempty(kind))
    kind=1;
end

if(nargin<4||isempty(ulim))
    ulim=1e6;
end

f=@(t)t.^n./(exp(t)-1);

if(kind==1)
    vals=arrayfun(@(t) integral1DAdaptive(f,[0,t]),x);
elseif(kind==2)
    vals=arrayfun(@(t) integral1DAdaptive(f,[t,ulim]),x);
end

%Enforce continuity
vals(x==0)=0;

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
