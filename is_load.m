function varargout = is_load(fileName,varargin)
% This file loads data and can be used in parfor loops
% I.S. 2013

for numoutputs = 1:length(varargin);
    varargout{numoutputs} = load(fileName,varargin{numoutputs});
    eval(['varargout{numoutputs} = varargout{numoutputs}.' varargin{numoutputs} ';'])
end

