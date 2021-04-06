% Function that returns the frequency response of convolution kernels
% Nikhil Chacko
% nchacko@ece.ucsb.edu

% $Author: $
% $LastChangedRevision: $
% $LastChangedDate: $
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
   
function h_hat = filter_hat(h,nu,h_arg)

if ~ischar(h)
    error('MATLAB:filter_hat:NonStringConvKernelType', ...
        'Usage: Convolution kernel type should be string');
else
    h=lower(h);
end

switch h
    case {'frt'}
        
        if ischar(h_arg)
            error('MATLAB:filter_hat:InvalidHArg', ...
                'Usage: h_arg should be chosen corresponding to h');
        else
            tau=h_arg(1);
            h_hat=exp(1j*pi/4)*(exp(-1j*pi*(tau*nu).^2));
        end
        
    case 'ifrt'
        
        if ischar(h_arg)
            error('MATLAB:filter_hat:InvalidHArg', ...
                'Usage: h_arg should be chosen corresponding to h');
        else
            tau=h_arg(1);
            h_hat=exp(-1j*pi/4)*(exp(1j*pi*(tau*nu).^2));
        end
        
    case 'rs'
        lambda=h_arg(1);
        z=h_arg(2);
        
        if (ischar(h_arg) || ~isreal(lambda) || ~isreal(z) || abs(lambda)<eps)
            error('MATLAB:filter_hat:InvalidHArg', ...
                'Usage: h_arg should be chosen corresponding to h');
        else
            h_hat=exp((1j*2*pi*z/lambda)*sqrt(1-((lambda*nu).^2)));
        end
        
    case {'identity','delta'}
        h_hat=1;
        
    case 'shift'
        
        if ~isscalar(h_arg) || ischar(h_arg) || ~isreal(h_arg)
            error('MATLAB:filter_hat:InvalidHArg', ...
                'Usage: h_arg should be chosen corresponding to h');
        else
            x0=h_arg(1);
            h_hat=exp(-1j*2*pi.*nu*x0);
        end
        
        
    otherwise
        error('MATLAB:filter_hat:UndefinedConvKernel', ...
            'Usage: Undefined convolution kernel');
end
