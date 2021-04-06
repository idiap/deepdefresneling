% FFT2_MIRROR is a function designed to compute the 2N-point FFT2 of
% a HS symmetric 2D matrix using a single N-point FFT2.
% 
% Usage:
% C_mirrorFFT2   = fft2_mirror(c);
% 
% where
% c:            the non-redundant quadrant of a HS symmetric 2D matrix
% C_mirrorFFT:  the 2N-point FFT2 of the HS symmetric 2D matrix
% See also: ifft2_mirror, fft_mirror, ifft_mirror
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
   

function C_mirrorFFT2 = fft2_mirror(c)

%% Input check
if isvector(c) || ~isnumeric(c) || ~ismatrix(c)
    error('Input should be a 2D numerical matrix');
end

[Nr,Nc] = size(c);


%% Re-arrange data
if (2*round(Nr/2) == Nr)                % Nr is even
    c_mirrorRowEven(1:Nr/2,:)          = downsample(c,2);
    c_mirrorRowEven(Nr/2+1:Nr,:)       = downsample(flipud(c),2);
    
else                                    % Nr is odd
    c_mirrorRowEven(1:(Nr+1)/2,:)       = downsample(c,2);
    c_mirrorRowEven((Nr+1)/2+1:Nr,:)    = downsample(flipud(c(1:Nr-1,:)),2);
end

if (2*round(Nc/2) == Nc)                % Nc is even
    c_mirrorEven(:,1:Nc/2)              = transpose(downsample(transpose(c_mirrorRowEven),2));
    c_mirrorEven(:,Nc/2+1:Nc)           = transpose(downsample(transpose(fliplr(c_mirrorRowEven)),2));
    
else                                    % Nr is odd
    c_mirrorEven(:,1:(Nc+1)/2)          = transpose(downsample(transpose(c_mirrorRowEven),2));
    c_mirrorEven(:,(Nc+1)/2+1:Nc)       = transpose(downsample(transpose(fliplr(c_mirrorRowEven(:,1:Nc-1))),2));
end

[kc0,kr0] = meshgrid(0:2*Nc-1,0:2*Nr-1);


%% 2N-point fft2 using a single N-point fft2
C_mirrorEvenFFT2          = fft2(c_mirrorEven, Nr, Nc);
C_mirrorFFT2              =  repmat(C_mirrorEvenFFT2,2,2) + ...
                            exp(1j*pi.*kc0/Nc).*repmat([C_mirrorEvenFFT2(:,1), fliplr(C_mirrorEvenFFT2(:,2:Nc))],2,2) + ...
                            exp(1j*pi.*kr0/Nr).*repmat([C_mirrorEvenFFT2(1,:); flipud(C_mirrorEvenFFT2(2:Nr,:))],2,2) + ...
                            exp(1j*pi.*kr0/Nr).*exp(1j*pi.*kc0/Nc).*...
                            repmat([C_mirrorEvenFFT2(1,1), fliplr(C_mirrorEvenFFT2(1,2:Nc));
                                    flipud(C_mirrorEvenFFT2(2:Nr,1)), rot90(C_mirrorEvenFFT2(2:Nr,2:Nc),2)],2,2);

end
