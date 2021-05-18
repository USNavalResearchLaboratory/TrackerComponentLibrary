function s=doublyCompSum(x)
%%DOUBLYCOMPSUM Sum a vector of numbers using doubly-compensated summation.
%           Compared to using the sum command, this function is less
%           susceptible to finite precision errors when summing many
%           numbers of varying magnitudes. Note that this functions does
%           NOT use extended precision arithmetic to improve the accuracy.
%
%INPUTS: x A 1XN or NX1 vectors of N values to sum.
%
%OUTPUTS: s The sum of the N numbers.
%
%The doubly compensated summation algorithm is given in Chapter 4.3 of [1].
%It oders the number and operations in a certain manner so as to reduce
%finite precision errors.
%
%EXAMPLE:
% x=[1;1e-17*ones(10000,1)];
% doubleCompSum(x)
% sum(x)
%One will see that doubleCompSum gives the sum as being larger than 1,
%whereas the sum command gives the sum as equal to 1.
%
%REFERENCES:
%[1] N. J. Higham, Accuracy and Stability of Numerical Algorithms.
%    Philadelphia: SIAM, 1996.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(x);
x=sort(x,'descend');

s=x(1);
c=0;
for k=2:n
    y=c+x(k);
    u=x(k)-(y-c);
    t=y+s;
    v=y-(t-s);
    z=u+v;
    s=t+z;
    c=z-(s-t);
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
