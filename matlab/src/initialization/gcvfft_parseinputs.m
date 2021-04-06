% Adapted from code originally written by Dr. Michael Liebling, UCSB
% Nikhil Chacko
% nchacko@ece.ucsb.edu

% $Author: nchacko $
% $LastChangedRevision: 5893 $
% $LastChangedDate: 2013-08-13 00:45:02 -0700 (Tue, 13 Aug 2013) $
% This file is distributed under the following license:
% ===========================================================================================
% Copyright 2013. The Regents of the University of California (Regents). All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification, are
% permitted for non-commercial educational, non-commercial research, and not-for-profit
% purposes only, provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this list of
% conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice, this list
% of conditions and the following disclaimer in the documentation and/or other materials
% provided with the distribution.
%
% 3. Neither the name of The Regents or the University of California nor the names of its
% contributors may be used to endorse or promote products derived from this software without
% specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
% EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
% SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
% BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
% WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% ===========================================================================================
   
function exp_params = gcvfft_parseinputs(varargin)
switch nargin
    case 1
        exp_params = gcvfft_initializeparams(); % initialize structure with default values
        if ~isempty(varargin{1}) % if user provides a structure with existing fields (even partial),
            % we replace the fields by those provided by user
            
            if ~iscell(varargin{1})
                error('gcvfft_parseinputs:inputisnotcell','The input to gcvfft_parseinputs must be a cell (e.g. varargin, or {p})')
            end
            userprovided_params = cell2mat(varargin{1}); % read the parameters provided by the user
            if ~isstruct(userprovided_params)
                error('gcvfft_parseinputs:inputisnotstructure','The provided input must be a structure(e.g. p, with p.xmin, etc.)')
            end
            exp_params = gcvfft_overwritestructfields(exp_params,userprovided_params);
        end
    otherwise
        error('gcvfft_parseinputs:wrongnumberofargs',['gcvfft_parseinputs expects one argument but ' num2str(nargin) ' provided'] )
end
end
