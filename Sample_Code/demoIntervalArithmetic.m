function demoIntervalArithmetic()
%%DEMOINTERVALARITHMETIC Demonstate how to use the Interval class for
%                        interval arithmetic. An interval defines a range,
%                        or set, of possible values rather than a single
%                        point. This allows calculations to be performed
%                        while keeping track of error bounds.
%
%
%November 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

disp('An example of how to use the Interval class for interval arithmetic') 

%%Creating variables
disp('Intervals define an upper and lower bound of possible values, rather than a single point.')
%Define lower and upper bound
a=Interval(-1,5); %Also notated as [-1,5]
%Use a single input for zero-span Interval, e.g. [1,1]
b=Interval(1);

%% Simple arithmetic
disp('Intervals will calculate the maximum and minimum values possible.')
%[-1,5]*[-1,5] = [-5,25]
disp('[-1,5]*[-1,5]')
a.*a

disp('Intervals calculate upper and lower bounds at floating-point level accuracy.')
disp('Width of [1,1]')
wid(b)
disp('Width of ([1,1] + eps/2)')
wid(b+eps/2)
disp('Interval can also handle infinite bounds when necessary')
disp('1/[0,5]')
1./Interval(0,5)

%Special case when dividing by Interval which crosses 0
%1/[-1,5] should give the bounds [-Inf,-1] and [.2,Inf]
%If only one output is provided, then the entire range is returned
disp('1/[-1,5] with one output')
1./a
%If two outputs are provided, then both Intervals are returned
disp('1/[-1 5] with two outputs')
[c,d]=1./a

%% Other mathematical functions
disp('Power and trig functions')
%Power functions, not intended for values of a<0
a=Interval(3,7); b=Interval(-2,2);
disp('[3,7]^2, [3,7]^[-2,2], [2^[-2,2]')
a.^2
a.^b
2.^b

%Trigonometric functions
disp('sin([-pi,pi/4]), tan([-pi/4,pi/2])')
sin(Interval(-pi,pi/4))
tan(Interval(-pi/4,pi/2))

%% Interval arrays
disp('Intervals can also be arrayed and indexed like any other numeric array.')
%Intervals can also be arrays of values
A=Interval(magic(4),2*magic(4));
%Index
disp('A(2,:)')
A(2,:)
%Assign
disp('A(3,4)=-A(3,4)')
A(3,4)=-A(3,4)
%Concatenate
disp('[A,A]')
[A,A]
%Transpose
disp('A.''')
A.'
%Size
disp('size(A)')
size(A)

disp('The only matrix math currently supported is multiplication.')
disp('A*A')
A*A

%% Set function
disp('Intervals can also be treated like sets to determine intersections and hulls.')
disp('Intersection and hull of [-2,4] and [3,7]')
intersection(Interval(-2,4),a)
hull(Interval(-2,4),a)

disp('If two Intervals do not intersect, an empty Interval is returned, represented by [NaN,NaN].')
disp('Intersection of [-2,2] and [3,7]')
intersection(a,b)
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
