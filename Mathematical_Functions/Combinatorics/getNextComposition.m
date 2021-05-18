function [theComp,recurData]=getNextComposition(param)
%%GETNEXTCOMPOSITION Get the next composition of the integer n. A
%           composition is a method of putting n unlabeled items into 1...n
%           labeled slots. Compositions are tuples of integers >=1 that sum
%           to n.
%
%INPUTS: param If the first composition is desired, then param is n.
%              Otherwise, param is the output recurData from the previous
%              call to this function.
%
%OUTPUTS: theComp The tX1 composition or an empty matrix is the final
%                 composition has been reached.
%       recurData A data structure that can be passed back to this function
%                 to get subsequent compositions.
%
%This function implements Problem 12a of Section 7.2.1.1 of [1].
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(~isstruct(param))
    n=param;
    
    s=zeros(n,1);%Maximum length.
    %Step C1, initialize.
    t=1;
    s(1)=n;
    
    theComp=s(1);   
    recurData.s=s;
    recurData.t=t;
    return
else
    recurData=param;
    s=recurData.s;
    t=recurData.t;
end

if(mod(t,2)==0)%Even step
    %Step C4
    
    if(s(t-1)>1)
        s(t-1)=s(t-1)-1;
        s(t+1)=s(t);
        s(t)=1;
        t=t+1;
    else
        t=t-1;
        if(t==1)
            theComp=[];
            return
        end
        
        s(t)=s(t+1);
        s(t-1)=s(t-1)+1;
    end
else%Odd Step
    %Step C3
    
    if(s(t)>1)
        s(t)=s(t)-1;
        s(t+1)=1;
        t=t+1;
    else
        t=t-1;
        s(t)=s(t)+1;
    end
end

%Step C2
theComp=s(1:t);   
recurData.s=s;
recurData.t=t;

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
