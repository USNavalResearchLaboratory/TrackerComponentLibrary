function xShifted=fracCircShift(x,n0)
%%FRACCIRCSHIFT Perform a circular shift of the elements in a matrix along
% a chosen direction. Unlike the circshift function that is built into
% Matlab, this function allows for non-integer shifts (fractional circular
% shifts). The interpolation needed for this is performed using the
% relationship between time domain shifts and multiplication in the
% frequency domain. Integer shifts produce the same results as circshift,
% subject to finite precision limitations.
%
%INPUTS: x The matrix to which a circular shift will be applied.
%       n0 The real number of steps that the vectors are to be circularly
%          shifted. This can be a non-integer and can be negative. Positive
%          shifts delay elements and negative shifts advance elements.
%
%OUTPUTS: xShifted The elements in all of the vectors in x circularly
%                  shifted by n0 in dimension dim, with appropriate
%                  interpolation for non-integer shifts.
%
%To make the shift direction clear, note that
%fracCircShift([1;2;3;4],1)=[4;1;2;3].
%fracCircShift([1;2;3;4],-1)=[2;3;4;1].
%The relationship between circular shifts and the discrete Fourier
%transform that is used to implement this function, is given in Equation
%5.99 in Chapter 5.7 of [1].
%
%EXAMPLE:
%Here, we use a non-integer shift:
% xShifted=fracCircShift([1;2;3;4],1.1)
%and see that xShifted=[4.1197;1.1932;1.8314;2.8557] when using five
%digits.
%
%REFERENCES:
%[1] S. K. Mitra, Digital Signal Processing: A Computer-Based Approach, 3rd
%    ed. Boston: McGraw Hill, 2006.
%
%September 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(x))
    xShifted=x;
    return;
end

N=size(x,1);
if(N<=1||n0==0)
    xShifted=x;
    return;
end

if(mod(N,2)==0)
    %If even.
    %Frequency bins k.
    k=[0:(N/2-1),(-N/2):(-1)].';
    W=exp(-1j*2*pi*n0*k/N);

    %Bin N/2 + 1 corresponds to the Nyquist frequency at k=-N/2. This is
    %set here to maintain Hermitian symmetry.
    W(N/2+1)=cos(pi*n0);
else
    %Frequency bins in k
    k=[0:((N-1)/2),(-(N-1)/2):-1].';
    W=exp(-1j*2*pi*n0*k/N);
end

%Perform circular shift.
xShifted=ifft(bsxfun(@times,fft(x),W));

if(isreal(x))
    %Deal with finite precision errors. Shifted real signals should remain
    %real.
    xShifted=real(xShifted);
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
