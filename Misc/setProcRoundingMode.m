function retVal=setProcRoundingMode(roundMode)
%%SETPROCROUNDINGMODE Set the rounding mode used in the processor. This
%                     specifies how results should be truncated to a
%                     finite number of bits when performing floating point
%                     arithmetic. Control of the rounding mode is useful
%                     when implementing routines utilizing interval
%                     algebra. Note that this does not appear to change
%                     the rounding mode used in mex files called from
%                     Matlab. When running things in parallel, the rounding
%                     mode will have to be set for each parallel process.
%                     For example, set it at the start of a parfor loop,
%                     not outside the loop.
%
%INPUTS: roundMode An integer specifying the rounding mode to use.
%                  Possible values are
%                  0 Rounding is done towards negative infinity.
%                  1 Rounding is done towards zero.
%                  2 Rounding is done to the nearest value.
%                  3 Rounding is done towards positive infinity.
%
%OUTPUTS: retVal This is zero if the rounding direction was successfully
%                set and a nonzero value otherwise.
%
%The rounding modes are standardized in [1]. The native Matlab version of
%this function calls the undocumented system_dependent('setround',dir)
%function in Matlab. The version that can be compiled in C just calls the
%standard C function fesetround.
%
%EXAMPLE:
% setProcRoundingMode(0);
% 1+eps(0)>1
% %The result should be false (0).
% setProcRoundingMode(3);
% 1+eps(0)>1
% %The result should be true (1).
%
%REFERENCES:
%[1] IEEE Standard 754-1985 for Binary Floating-point Arithmetic, IEEE,
%    (1985).
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

switch(roundMode)
    case 0%Rounding is done towards negative infinity.
        system_dependent('setround',-Inf);
    case 1%Rounding is done towards zero.
        system_dependent('setround',0);
    case 2%Rounding is done to the nearest value.
        system_dependent('setround',0.5);
    case 3%Rounding is done towards positive infinity.
        system_dependent('setround',Inf);
    otherwise
        error('Unknown rounding mode specified.')
end

if(nargout>0)
    actualMode=getProcRoundingMode();
    retVal=(actualMode~=roundMode);
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
