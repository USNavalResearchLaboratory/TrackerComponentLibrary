function rou = rootsOfUnity(N)
%%ROOTSOFUNITY Generates 2D vectors corresponding to the Nth roots of
%              unity.
%INPUT: N An integer determing which roots of unity to generate.
%
%OUTPUT: rou A 2-by-numRoots matrix containing the collection of 2D vectors
%            representing the Nth roots of unity.
%
%EXAMPLE:
%Plots the 3rd and 8th roots of unity.
% rou3 = rootsOfUnity(3);
% rou8 = rootsOfUnity(8);
% figure(1); clf
% scatter(rou3(1,:),rou3(2,:),'o','filled')
% hold on 
% scatter(rou8(1,:),rou8(2,:),'d')
% grid on
%
%March 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

rou = zeros(2,N);
deltaPhi = 2*pi/N;
phi = deltaPhi;
for n = 1:N
    rou(:,n) = [cos(phi);sin(phi)];
    phi = phi+deltaPhi;
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
