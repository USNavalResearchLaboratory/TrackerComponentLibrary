function V=crossProdMat(vec,side)
%%CROSSPRODMAT Get a matrix so that the cross product between two vectors
%              can be expressed in terms of a matrix multiplication.
%              if one wishes to evaluate vec X b as V*b, where V is a
%              matrix dependent on vec, then the side option can be
%              omitted or set to left (the default value). If one wishes
%              to evaluate b X vec as V*b, then the side option should be
%              set to 'right'. This produces cross product matrices for
%              right-handed 3D cross products and for 7D cross products
%              using the 7D cross product defined in terms of an
%              orthonormal basis using asymmetry as in [1].
%
%INPUTS: vec  The 3X1 or 7X1 vector in a cross product relation that one
%             wishes to transform into a matrix.
%       side  An optional parameter specifying the side that vec is on in
%             the cross product relation. Possible values are
%             'left' (The default value if omitted) The cross product
%                    relation is of the form vec X b, where b is another
%                    3X1 vector.
%             'right' The cross product relation is of the form b X vec.
%
%OUTPUTS: V  If side is omitted or is 'left', then V*b is the same as v X b
%            (cross(v,b)). If side is 'right, then V*b is the same as bXv.
%
%Cross products for vectors other than 3D and 7D do not exist, as discussed
%in [2]. See the comments to the cross7 function for more information on 7D
%cross products.
%
%REFERENCES:
%[1] P. Lounesto, "Octonians and triality," Advances in Clifford Algebras,
%    vol. 11, no. 2, pp. 191-213, Dec. 2001.
%[2] W. S. Massey, "Cross products of vectors in higher dimensional
%    Euclidean spaces," The American Mathematical Monthly, vol. 90, no. 10,
%    pp. 697-701, Dec. 1983.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    side='left';
end

if(size(vec,1)==3)%3D cross product
    V=[+0,      -vec(3),    +vec(2);
       +vec(3), +0,         -vec(1);
       -vec(2), +vec(1),    +0];
else%7D cross product
    V=[+0,      +vec(4),    +vec(7),    -vec(2),    +vec(6),    -vec(5),    -vec(3);
       -vec(4), +0,         +vec(5),    +vec(1),    -vec(3),    +vec(7),    -vec(6);
       -vec(7), -vec(5),    +0,         +vec(6),    +vec(2),    -vec(4),    +vec(1);
       +vec(2), -vec(1),    -vec(6),    0,          +vec(7),    +vec(3),    -vec(5);
       -vec(6), +vec(3),    -vec(2),    -vec(7),    0,          +vec(1),    +vec(4);
       +vec(5), -vec(7),    +vec(4),    -vec(3),    -vec(1),    0,          +vec(2);
       +vec(3), +vec(6),    -vec(1),    +vec(5),    -vec(4),    -vec(2),    0];
end

switch(side)
    case 'left'
    case 'right'
        V=V';
    otherwise
        error('Incorrect side value specified.');
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
