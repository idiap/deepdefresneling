% Adapted from code originally written by Dr. Michael Liebling, UCSB
% Nikhil Chacko
% nchacko@ece.ucsb.edu

% $Author: nchacko $
% $LastChangedRevision: 5924 $
% $LastChangedDate: 2013-08-22 23:33:49 -0700 (Thu, 22 Aug 2013) $
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
   
function exp_params = gcvfft_initializeparams()
% Algorithm parameters

exp_params.T1               = [1,1];
exp_params.h                = 'FrT';
exp_params.h_arg            = [1,0];
exp_params.phi1             = 'cardinal-bspline';
exp_params.deg1             = 3;
exp_params.phi2             = 'delta';
exp_params.deg2             = 1;
exp_params.p                = 1;
exp_params.q                = 1;
exp_params.gamma            = 1;
exp_params.K                = 500;

end
