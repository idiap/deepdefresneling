% IFFT2_MIRROR is the inverse counterpart to fft2_mirror 
% Usage:
% d_mirrorIFFT2   = ifft2_mirror(D);
% 
% where 
% D:                the 2N-point FFT of a HS symmetric 2D matrix
% d_mirrorIFFT2:    the N-point non-redundant quadrant of the 2N-point IFFT2 of D
% 
% See also: fft2_mirror, fft_mirror, ifft_mirror
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

function d_mirrorIFFT2 = ifft2_mirror(D)

%% Input check
if isvector(D) || ~isnumeric(D) || ~ismatrix(D)
    error('Input should be a 2D numerical matrix');
end

[twoNr,twoNc] = size(D);

if ((2*round(twoNr/2) ~= twoNr) || (2*round(twoNc/2) ~= twoNc))
    error('The number of rows and columns in input should be even')
end

Nr = twoNr/2;
Nc = twoNc/2;

if (any(abs(D(Nr+1,:)) > eps*10^(6)) || (any(abs(D(:,Nc+1)) > eps*10^(6))))
    warning('The FFT2 coefficients in the middle row and column should be zero for half-sample mirror-symmetric signals')
end


%% 2N-point ifft2 using a single N-point ifft2
E               = 0.5*(D + circshift(D,[Nr,Nc]));
D               = 0.5*(E(1:Nr,1:Nc) + E(Nr+1:end, 1:Nc));
d_mirrorEven    = ifft2(D, Nr, Nc);


%% Re-arrange data
if (2*round(Nr/2) == Nr)        % Nr is even
    d_mirrorIFFT2_row(1:2:Nr,:)     = d_mirrorEven(1:Nr/2,:);
    d_mirrorIFFT2_row(2:2:Nr,:)     = d_mirrorEven(Nr:-1:Nr/2+1,:);
    
else                            % Nr is odd
    d_mirrorIFFT2_row(1:2:Nr,:)     = d_mirrorEven(1:(Nr+1)/2,:);
    d_mirrorIFFT2_row(2:2:Nr,:)     = d_mirrorEven(Nr:-1:(Nr+1)/2+1,:);
end

if (2*round(Nc/2) == Nc)        % Nc is even
    d_mirrorIFFT2(:,1:2:Nc)     = d_mirrorIFFT2_row(:,1:Nc/2);
    d_mirrorIFFT2(:,2:2:Nc)     = d_mirrorIFFT2_row(:,Nc:-1:Nc/2+1);
    
else                            % Nc is odd
    d_mirrorIFFT2(:,1:2:Nc)     = d_mirrorIFFT2_row(:,1:(Nc+1)/2);
    d_mirrorIFFT2(:,2:2:Nc)     = d_mirrorIFFT2_row(:,Nc:-1:(Nc+1)/2+1);
end

end
