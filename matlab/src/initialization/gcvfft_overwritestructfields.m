% OVERWRITESTRUCTFIELDS
% overwrites the fields in a structure with values provided in corresponding fields in another structure
% Usage:
%    sout = overwritestructfields(s,snew)
%    overwrites the values of the subset of fields in s that are in snew.
%    All fields of snew must exist in s, but s can have more fields than snew.
% Example:
%    s.a = 1;
%    s.b = 'a';
%    s.c = [1 2 3];
%    snew.b = 'q';
%    sout = overwritestructfields(s,snew)
%    sout = 
%
%       a: 1
%       b: 'q'
%       c: [1 2 3]
% 
%    Michael Liebling 
%    UCSB

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
   
function s = gcvfft_overwritestructfields(s,snew)
	if ~isstruct(s)
	    error('overwritestructfields:inputoldisnotstructure','The provided input must be a structure')	
	end
	if ~isstruct(snew)
	    error('overwritestructfields:inputnewisnotstructure','The provided input must be a structure')	
	end
	fieldsToUpdate = fieldnames(snew); % retrieve field names in the new field 
	for k = 1:length(fieldsToUpdate) % overwrite the default parameters by user provided ones
		if isfield(s,fieldsToUpdate{k})
	    	s.(fieldsToUpdate{k}) = snew.(fieldsToUpdate{k});
		else
			warning('overwritestructfields:providedfielddoesnotexist',...
                 ['The provided field ''' ...
                   fieldsToUpdate{k} ... 
                   ''' is not present in structure to be overwritten. Only existing fields will be overwritten.'] )	
		end
	end
end
