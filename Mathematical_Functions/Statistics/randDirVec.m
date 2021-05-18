function u=randDirVec(numDim,N)
%%RANDDIRVEC Generate a uniformly distributed random direction vector
%            (unit vector) in an arbitrary number of dimensions.
%
%INPUTS: numDim The number of dimensions of the random unit vector. This
%               is >=1.
%             N The number of random direction vectors to generate. If
%               this parameter is omitted, then N=1 is used.
%
%OUTPUTS: u A numDimXN matrix of N numDim-dimensional unit vectors that
%           are uniformly distributed on the unit-sphere (or hypersphere).
%
%One cannot simply generate a uniformly-distributed direction vector by
%setting each dimensions to a Uniform(-1,1) random variable and then
%normalizing the result. Rather, following the algorithm in Chapter 3.4.1,
%E6 of [1], the elements of the vector must be generated as normal 0-1
%random variables and then the vector normalized.
%
%REFERENCES:
%[1] D. Knuth, The Art of Computer Programming: Seminumerical Algorithms,
%    3rd ed. Reading, MA: Addison-Wesley, 1998, vol. 2.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    N=1;
end

u=zeros(numDim,N);

for curN=1:N
    %Knuth, Chapter 3.4.1, E6 (Random point on an n-dimensional sphere with
    %radius one. 
    u(:,curN)=randn(numDim,1);
    u(:,curN)=u(:,curN)/norm(u(:,curN));

    %Deal with the instance that norm(u)=0, which will occur with an extremely
    %small probability due to the limited precision of the random number
    %generation.
    if(~isfinite(u(:,curN)))
       u(:,curN)=zeros(numDim,1);
       u(1,curN)=1;
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
