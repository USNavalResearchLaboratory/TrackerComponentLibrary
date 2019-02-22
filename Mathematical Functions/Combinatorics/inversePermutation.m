function xPerm=inversePermutation(xPerm,algorithm)
%%INVERSEPERMUTATION Given a permutation of the number 1:N, return the
%           inverse of the permutation. If something were rearranged with
%           the original permutation, the inverse permutation puts it back
%           into the original order.
%
%INPUTS: xPerm A 1XN or NX1 vector containing a permutation of the numbers
%              from 1 to N.
%    algorithm An optional input specifying the algorithm to use. Possible
%              values are:
%              0 (The default if omitted or an empty matrix is passed) Use
%                algorithm I from Chapter 1.3.3 of [1]. This is an in-place
%                inversion.
%              1 Use algorithm J from Chapter 1.3.3 of [1]. This is an
%                in-place inversion.
%              3 Explicitely create the inversion as xInvPerm(xPerm(k))=k
%                for k=1:N. This is not an in-place inversion.
%
%OUTPUTS: xInvPerm The NX1 or 1XN (matching the input) inverse of the
%                  permutation in xPerm.
%
%EXAMPLE:
% xPerm=[2;5;1;4;3];
% y=1:5;
% z=y(xPerm);
% xInvPerm=inversePermutation(xPerm);
% all(z(xInvPerm)==y)
%The result of the final line will  be 1 (true).
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 1: Fundamental
%    Algorithms, 3rd Edition, Boston: Addison-Wesley, 2011.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=0;
end

switch(algorithm)
    case 0%Algorithm I in [1].
        n=length(xPerm);

        %Step I1
        m=n;
        j=-1;

        while(m>0)
            %Step I2
            i=xPerm(m);
            if(i>=0)
                while(1)
                    %Step I3
                    xPerm(m)=j;
                    j=-m;
                    m=i;
                    i=xPerm(m);

                    %Step I4
                    if(i>0)
                        continue;
                    end

                    i=j;
                    break;
                end
            end

            %Step I5
            xPerm(m)=-i;
            m=m-1;
        end
    case 1%Algorithm J in [1].
        n=length(xPerm);
        
        %Step J1;
        xPerm=-xPerm;
        m=n;
        
        while(m>0)
            %Step J2
            j=m;

            %Step J3
            i=xPerm(j);
            while(i>0)
                j=i;
                i=xPerm(j);
            end

            %Step J4
            xPerm(j)=xPerm(-i);
            xPerm(-i)=m;

            %Step J5
            m=m-1;
        end
    case 2%A non-in-place algorithm.
        n=length(xPerm);
        xInvPerm=zeros(size(xPerm));
        for k=1:n
            xInvPerm(xPerm(k))=k; 
        end
        xPerm=xInvPerm;
    otherwise
        error('Unknown algorithm specified.')
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
