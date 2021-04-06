% FFT_MIRROR is a function designed to compute the 2N-point FFT of
% a HS symmetric vector using a single N-point FFT.
% 
% Usage:
% C_mirrorFFT   = fft_mirror(c);
% 
% where
% c is the non-redundant half of the HS symmetric input coefficients
% C_mirrorFFT is the 2N-point FFT of the HS symmetric sequence
% See also: ifft_mirror, fft2_mirror, ifft2_mirror
% 
% Nikhil Chacko
% August 10, 2013
% nchacko@ece.ucsb.edu

% $Author: nchacko $
% $LastChangedRevision: 5914 $
% $LastChangedDate: 2013-08-19 19:37:40 -0700 (Mon, 19 Aug 2013) $
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

function C_mirrorFFT = fft_mirror(c)

%% Input check
if ~isvector(c) || ~isnumeric(c)
    error('Input should be a 1D numerical vector');
end

N = length(c);

if N<2
    C_mirrorFFT = c;
    return;
end

reshapeNeeded=0;
if ~isrow(c)
    c = transpose(c);
    reshapeNeeded=1;
end


%% Re-arrange data
if (2*round(N/2) == N)          % N is even
    c_mirrorEven(1:N/2)          = downsample(c,2);
    c_mirrorEven(N/2+1:N)        = downsample(fliplr(c),2);
    
else                            % N is odd
    c_mirrorEven(1:(N+1)/2)      = downsample(c,2);
    c_mirrorEven((N+1)/2+1:N)    = downsample(fliplr(c(1:N-1)),2);
end


%% 2N-point fft using a single N-point fft
k0 = 0:2*N-1;
C_mirrorEvenFFT          = fft(c_mirrorEven, N);
C_mirrorFFT              =  repmat(C_mirrorEvenFFT,1,2) + exp(1j*pi.*k0/N).* ...
                            repmat([C_mirrorEvenFFT(1), fliplr(C_mirrorEvenFFT(2:N))],1,2);


%% Final reshape, if necessary
if reshapeNeeded
    C_mirrorFFT = transpose(C_mirrorFFT);
end
end
