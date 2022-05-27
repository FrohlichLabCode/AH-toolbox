function getStructVars(struct)
% This is normally not suggested, but could be convenient when there are a
% lot of variables to load from a struct.
% AH 2020/11

fieldName = fieldnames(struct);
for ifield=1:length(fieldName)
    eval([fieldName{ifield} '=' ...,
    'struct.(fieldName{ifield});' ]);
end
end