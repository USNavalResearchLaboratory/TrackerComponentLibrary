classdef KeyVal
%KEYVAL     A class that encapsulates a key and a value, which could be
%           an instance of a handle class. It is assumed that standard
%           comparison operations can be performed on the keys. This class
%           is useful for inserting elements into a BinaryHeap or other
%           data structure that requires a key and a value. 
%
%DEPENDENCIES: NONE
%
%Instances of the KeyVal class can be compared to each other by key or to
%other objects that have the same type of key. For example,
%keyVal1<keyVal2 will compare the two objects by their keys.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    properties
       key
       value
    end
    methods
        function newKeyVal=KeyVal(key,val)
        %KEYVAL The constructor method.
            if(nargin<2)
                val=[];
            end
            if(nargin<1)
                key=[]; 
            end
            
            newKeyVal.key=key;
            newKeyVal.value=val;
        end
        
        function val=lt(obj1,obj2)
        %%LT Less than comparison  
          
            if(isa(obj1,'KeyVal')==true&&isa(obj2,'KeyVal')==true)
                val=obj1.key<obj2.key;
            elseif(isa(obj1,'KeyVal')==true)
                val=obj1.key<obj2;
            else
                val=obj1<obj2.key;
            end
        end
      
        function val=gt(obj1,obj2)
        %%GT Greater than comparison
        
            if(isa(obj1,'KeyVal')==true&&isa(obj2,'KeyVal')==true)
                val=obj1.key>obj2.key;
            elseif(isa(obj1,'KeyVal')==true)
                val=obj1.key>obj2;
            else
                val=obj1>obj2.key;
            end
        end
      
        function val=le(obj1,obj2)
        %%LE Less than or equal to comparison
        
            if(isa(obj1,'KeyVal')==true&&isa(obj2,'KeyVal')==true)
                val=obj1.key<=obj2.key;
            elseif(isa(obj1,'KeyVal')==true)
                val=obj1.key<=obj2;
            else
                val=obj1<=obj2.key;
            end
        end
      
        function val=ge(obj1,obj2)
        %%GE Greater than or equal to comparison
        
            if(isa(obj1,'KeyVal')==true&&isa(obj2,'KeyVal')==true)
                val=obj1.key>=obj2.key;
            elseif(isa(obj1,'KeyVal')==true)
                val=obj1.key>=obj2;
            else
                val=obj1>=obj2.key;
            end
        end
      
        function val=ne(obj1,obj2)
        %%NE Not equal to comparison
        
            if(isa(obj1,'KeyVal')==true&&isa(obj2,'KeyVal')==true)
                val=obj1.key~=obj2.key;
            elseif(isa(obj1,'KeyVal')==true)
                val=obj1.key~=obj2;
            else
                val=obj1~=obj2.key;
            end
        end
      
        function val=eq(obj1,obj2)
        %%EQ Equal to comparison
        
            if(isa(obj1,'KeyVal')==true&&isa(obj2,'KeyVal')==true)
                val=obj1.key==obj2.key;
            elseif(isa(obj1,'KeyVal')==true)
                val=obj1.key==obj2;
            else
                val=obj1==obj2.key;
            end
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
