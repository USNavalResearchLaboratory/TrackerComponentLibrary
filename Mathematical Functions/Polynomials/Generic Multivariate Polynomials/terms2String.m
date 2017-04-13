function curPoly=terms2String(termMat,xVar,maxDigits,complexFormat)
%%TERMS2STRING Convert a matrix of monomial terms for a multivariate
%              polynomial into a string that is suitable for display or for
%              using as an entry in a cell array in the input to the
%              function solvePolySysWithExtProg.
%
%INPUTS: termMat An (n+1)XnumTerms matrix such that
%           termMat(:,i)=[c,a1,a2,...,an] is a monomial term where c is the
%           value of of the monomial coefficient and the monomial is
%           x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1).
%           n is the highest numbered used to index a variable in
%           theString, or is numVars, whichever is larger. Note that there
%           should not be repeated monomials.
%      xVar An nX1 or 1Xn cell array where xVar{j} corresponds to the
%           variable for coefficient aj in termMat. These are the names of
%           the variables given as character strings. The names cannot
%           begin with a number.
% maxDigits The maximum number of digits to display for each number. The
%           default if this parameter is omitted or an empty matrix is
%           passed is 16, which is sufficient for floating point numbers.
% complexFormat A string indicating the value to use to specify an
%           imaginary number. By default, this is i, which is suitable for
%           Bertini and PHCpack. However, one might want to use ii if
%           complex numbers are to be used in Macaulay2.
%
%OUTPUTS: curPoly A string representing the polynomial expressed in termMat
%                 using xVar as the variables. Note that this function does
%                 not combine like terms.
%
%EXAMPLE:
% termMat=[1, -2, -6, -1,-4, 12;
%          1,  0,  1,  0, 0,  0;
%          0,  1,  0,  1, 0,  0;
%          1,  1,  0,  0, 0,  0;
%          0,  0,  0,  0, 1,  0];
% xVar={'xa','xb','xc','xd'};
% curPoly=termMat2String(termMat,xVar)
% %The result is
% curPoly=xa*xc-2*xb*xc-6*xa-xb-4*xd+12
%
%March 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(maxDigits))
    maxDigits=16;
end

if(nargin<4||isempty(complexFormat))
    complexFormat='i';
end

numTerms=size(termMat,2);
numVars=size(termMat,1)-1;

curPoly=[];
for curTerm=1:numTerms
    curCoeff=[];
    for curVar=1:numVars
        exponentVal=termMat(curVar+1,curTerm);
        if(exponentVal~=0)
            if(isempty(curCoeff)&&exponentVal==1)
                curCoeff=xVar{curVar};
            elseif(isempty(curCoeff))
                curCoeff=[xVar{curVar},'^',num2str(exponentVal)];
            elseif(exponentVal==1)
                curCoeff=[curCoeff,'*',xVar{curVar}];
            else
                curCoeff=[curCoeff,'*',xVar{curVar},'^',num2str(exponentVal)];
            end
        end
    end
        
    if(isempty(curPoly))
        if(isempty(curCoeff))%If the term is a constant.
            if(real(termMat(1,curTerm))==0)%If it is purely imaginary
                curPoly=[num2str(imag(termMat(1,curTerm)),maxDigits),'*1i'];
            else
                curPoly=num2str(termMat(1,curTerm),maxDigits);
            end
        else
            if(isreal(termMat(1,curTerm)))%Purely real
                if(termMat(1,curTerm)==1)
                    curPoly=curCoeff;
                elseif(termMat(1,curTerm)==-1)
                    curPoly=['-',curCoeff];
                else
                    curPoly=[num2str(termMat(1,curTerm),maxDigits),'*',curCoeff];
                end
            elseif(real(termMat(1,curTerm))==0)%Purely imaginary
                val=imag(termMat(1,curTerm));
                if(val==1)
                    curPoly=['1i*',curCoeff];
                elseif(termMat(1,curTerm)==-1)
                    curPoly=['-1i*',curCoeff];
                else
                    curPoly=[num2str(termMat(1,curTerm),maxDigits),'*',curCoeff];
                end
            else
                curPoly=['(',num2str(termMat(1,curTerm),maxDigits),')*',curCoeff];
            end
        end
    else
        if(isempty(curCoeff))%If the term is a constant
            if(isreal(termMat(1,curTerm)))
                if(termMat(1,curTerm)>0)
                    curPoly=[curPoly,'+',num2str(termMat(1,curTerm),maxDigits)];
                else
                    curPoly=[curPoly,'-',num2str(-termMat(1,curTerm),maxDigits)];
                end
            elseif(real(termMat(1,curTerm))==0)%Purely imaginary
                val=imag(termMat(1,curTerm));
                
                if(val>0)
                    curPoly=[curPoly,'+',num2str(val,maxDigits),'*1i'];
                else
                    curPoly=[curPoly,'-',num2str(-val,maxDigits),'*1i'];
                end
            else
                curPoly=[curPoly,'+(',num2str(termMat(1,curTerm),maxDigits),')'];
            end
        else
            if(isreal(termMat(1,curTerm)))
                if(termMat(1,curTerm)==1)
                    curPoly=[curPoly,'+',curCoeff];
                elseif(termMat(1,curTerm)==-1)
                    curPoly=[curPoly,'-',curCoeff];
                else
                    if(termMat(1,curTerm)>0)
                        curPoly=[curPoly,'+',num2str(termMat(1,curTerm),maxDigits),'*',curCoeff];
                    else
                        curPoly=[curPoly,'-',num2str(-termMat(1,curTerm),maxDigits),'*',curCoeff];
                    end
                end
            elseif(real(termMat(1,curTerm))==0)%Purely imaginary
                val=imag(termMat(1,curTerm));
                
                if(val==1)
                    curPoly=[curPoly,'+1i*',curCoeff];
                elseif(val==-1)
                    curPoly=[curPoly,'-1i*',curCoeff];
                else
                    if(val>0)
                        curPoly=[curPoly,'+',num2str(val,16),'*1i*',curCoeff];
                    else
                        curPoly=[curPoly,'-',num2str(-val,16),'*1i*',curCoeff];
                    end
                end
            else
                curPoly=[curPoly,'+(',num2str(termMat(1,curTerm),16),')*',curCoeff];
            end
        end
    end
end

%Replace any complex numbers with the proper string.
curPoly=strrep(curPoly,'1i',complexFormat);
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
