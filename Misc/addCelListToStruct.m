function opts=addCelListToStruct(opts,args)
%%ADDCELLLISTTOSTRUCT Given a cell array full of string-value pairs, such
%          as one might obtain using a varargin input to a function, take
%          the structure opts and set opts.stringName=value for all string
%          names and value pairs in args.
%
%INPUTS: opts A structure to which the string name-value pairs in args are
%             to be added. opts.stringName must already exist for all
%             possible string names. Typically, these will be set with
%             default values.
%        args A cell array where ars{2*i+1} is a string and args{2*i+2} is
%             the value associated with the string.
%
%OUTPUTS: opts The structure with the elements corresponding to string
%              names in args changed to the values that are also provided.
%
%This function is most useful if one stores default values for a function
%in a structure opts and then lets the user change individual entries with
%a varargin input.
%
%EXAMPLE:
%Consider the following structure, which holds default values where the
%second one is to be changed to 2:
% opts=[];
% opts.one=1;
% opts.two=20;
% args=cell(2,1);
% args{1}='two';
% args{2}=2;
% opts=addCelListToStruct(opts,args)
%The new opts variable has opts.one=1 and opts.two=2.
%
%February 2016 Hatim F. Alqadah, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if nargin < 2
        error('Invalid number of arguments');
    end
    
    if(isempty(args))
       return; 
    end
    
    assert(isstruct(opts),'opts should be a struct');
    assert(iscell(args),'args should be a cell array');
    assert(isvector(args),'args should be a vector');
    n = length(args);
    if(mod(n,2)~=0)
        error('args should be given as field-value pairs');
    end
    
    %Get the option names as specifed in opts
    optnames=fieldnames(opts);
    for arg=1:2:n
        name = args{arg};
        if(~ischar(name))
             warning('MATLAB:invalid_fieldname','encountered a non-string option name');
        end
        %Look for a match in the field names
        idx = find(strcmpi(optnames,name),1);
        if(~isempty(idx))
            opts.(optnames{idx})=args{arg+1};
        else
            warning('MATLAB:invalid_option',['unrecognized option ',name,' skipping...']);
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
