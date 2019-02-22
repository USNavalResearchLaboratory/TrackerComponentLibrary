function convMats=convertCovType(covMats,origType,desType)
%%CONVERTCOVTYPE Covariance matrices used in target tracking tend to take
%          different forms. Some trackers use standard covariance matrices
%          while others use square-root lower triangular versions, inverse
%          covariance matrices, or lower-triangular square-root inverse
%          matrices. This function takes a hypermatrix of many covariance
%          matrices of one type and converts all of them to a different
%          type.
%
%INPUTS: covMats A numDImXnumDimXN set of N covariance matrices of the
%                original type that are to be changed.
% origType,desType These are strings that specify the original and
%                destination types. Possible values are:
%                'LT'    A lower-triangular square root covariance matrix.
%                'std'   A standard covariance matrix (not transformed).
%                'inv'   An inverse covariance matrix.
%                'LTInv' A lower-triangular square root inverse covariance
%                        matrix.
%
%OUTPUTS: convMats The converted covariance matrices.
%
%THis function just determines the proper function for the conversion and
%applies it to all of the matrices using applyFunToEachMatrix.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

switch(origType)
    case 'std'
        switch(desType)
            case 'std'
                convMats=covMats;
                return;
            case 'LT'
                f=@(X)cholSemiDef(X,'lower');
            case 'inv'
                f=@(X)inv(X);
            case 'LTInv'
                f=@(X)chol(inv(X),'lower');
            otherwise
                error('Unknown destination covariance matrix type specified')
        end
    case 'LT'
        switch(desType)
            case 'std'
                f=@(X)X*X';
            case 'LT'
                convMats=covMats;
                return;
            case 'inv'
                f=@(X)inv(X*X');
            case 'LTInv'
                f=@(X)chol(inv(X*X'),'lower');
            otherwise
                error('Unknown destination covariance matrix type specified')
        end
    case 'inv'
        switch(desType)
            case 'std'
                f=@(X)inv(X);
            case 'LT'
                f=@(X)chol(inv(X),'lower');
            case 'inv'
                convMats=covMats;
                return;
            case 'LTInv'
                f=@(X)cholSemiDef(X,'lower');
            otherwise
                error('Unknown destination covariance matrix type specified')
        end
    case 'LTInv'
        switch(desType)
            case 'std'
                f=@(X)inv(X*X');
            case 'LT'
                f=@(X)chol(inv(X*X'),'lower');
            case 'inv'
                f=@(X)X*X';
            case 'LTInv'
                convMats=covMats;
                return;
            otherwise
                error('Unknown destination covariance matrix type specified')
        end
    otherwise
        error('Unknown original covariance matrix type specified')
end

matDims=[size(covMats,1),size(covMats,2)];
convMats=applyFunToEachMatrix(f,covMats,[],matDims);

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
