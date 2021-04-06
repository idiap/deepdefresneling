% IFFT_MIRROR is the inverse counterpart to fft_mirror 
% Usage:
% d_mirrorIFFT   = ifft_mirror(D);
% 
% where 
% D:            the 2N-point FFT of a HS symmetric sequence
% d_mirrorIFFT: the N-point non-redundant half of the 2N-point IFFT of D
% 
% See also: fft_mirror, fft2_mirror, ifft2_mirror
% 
% Nikhil Chacko
% August 10, 2013
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

function d_mirrorIFFT = ifft_mirror(D)

%% Input check
if ~isvector(D) || ~isnumeric(D)
    error('Input should be a 1D numerical vector');
end

twoN = length(D);

if (2*round(twoN/2) ~= twoN) 
    error('Length of input vector should be even')
end

N = twoN/2;

if N<2
    d_mirrorIFFT = D;
    return;
end

if abs(D(N+1)) > eps*10^(6);
    warning('The FFT coefficient in the middle index should be zero for half-sample mirror-symmetric signals')
end

reshapeNeeded=0;
if ~isrow(D)
    D = transpose(D);
    reshapeNeeded=1;
end


%% 2N-point ifft using a single N-point ifft
D = 0.5*(D(1:N) + D(N+1:end));
d_mirrorEven = ifft(D, N);


%% Re-arrange data
if (2*round(N/2) == N)      % N is even
    d_mirrorIFFT(1:2:N)     = d_mirrorEven(1:N/2);
    d_mirrorIFFT(2:2:N)     = d_mirrorEven(N:-1:N/2+1);
    
else                        % N is odd
    d_mirrorIFFT(1:2:N)     = d_mirrorEven(1:(N+1)/2);
    d_mirrorIFFT(2:2:N)     = d_mirrorEven(N:-1:(N+1)/2+1);
end


%% Final reshape, if necessary            
if reshapeNeeded
    d_mirrorIFFT = transpose(d_mirrorIFFT);
end
end
