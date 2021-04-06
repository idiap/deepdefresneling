function M2=stretchcontrast(M1,newMin,newMax,varargin)
%STRETCHCONTRAST stretches the contrast to the specified values.
%	Usage:
%		M2=stretchcontrast(M1,mi,ma);
%		Computes the min and max values of the image M1 and rescales it
%		linearly in order that it covers the range [mi,max];
%		
%		M2=stretchcontrast(M1,newmin,newmax,refmin,refmax);
%		Rescales the image linearly such that the refmin value
%		becomes newmin and refmax becomes newmax.
%	
%	Send questions and comments to:
%	Michael Liebling 
%	Biomedical Imaging Group
%	IOA-STI-EPFL
%	CH-1015 Lausanne, Switzerland
%	michael.liebling AT epfl.ch
%	Michael Liebling, 7 november 2002

%   $LastChangedRevision: 3104 $
%   $LastChangedDate: 2011-05-31 12:09:16 -0700 (Tue, 31 May 2011) $
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

  
nin=nargin-3;

switch nin
	case 0
		oldMin=min(M1(:));
		oldRange=range(M1(:));
		if oldRange==0
			%	error('Range of input image is zero.')
			M2=M1;% if range is zero: do nothing. MLg 1 dec 2003
			return; 
		end
	case 2
		oldMin=varargin{1};
		oldMax=varargin{2};
		if (oldMax>=oldMin)
			oldRange=oldMax-oldMin;
		else
			error('The provided old max value is less than the min value.');
		end
	otherwise
		error('Wrong number of input arguments');
end
if (newMax>=newMin)
	newRange=newMax-newMin;
else
	error('The requested final max value is less than the min value.');
end
M2=((M1-oldMin)/oldRange)*newRange+newMin;

end % end stretchcontrast.m
