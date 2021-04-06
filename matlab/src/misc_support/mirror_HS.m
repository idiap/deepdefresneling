% A function that creates a HS symmetric version of input.
% Primarily used during debugging and testing phase.
% 
% Nikhil Chacko,
% August 10, 2013

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
   
function cm = mirror_HS(c)

if isvector(c)
    if isrow(c)
        cm = [c, fliplr(c)];
    else
    cm = [c; flipud(c)];
    end
    
else
    cm = [c,            fliplr(c);
          flipud(c),    rot90(c,2)];
end
end
