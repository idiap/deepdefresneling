% Function that returns the frequency response of basis functions.
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
   
function y=phi_hat(phi,nu,varargin)

numvarargs = length(varargin);

%Set defaults for optional inputs: deg, gamma
optargs={0, 1};

%Over-ride defaults if specified explicitly in function call
optargs(1:numvarargs)=varargin;

[deg, gamma]=optargs{:};

if ~ischar(phi)
    error('MATLAB:phi_hat:NonStringBasisType', ...
        'Usage: Basis type should be string');
else
    phi=lower(phi);
end


switch phi
    
    case {'causal-bspline'}
        y=(exp(-1j*pi*nu).*sinc(nu)).^(deg+1);
               
    case {'b-spline', 'bspline'}        
        y=(sinc(nu)).^(deg+1);
    
    case {'bspline-derivative'}
        y=2*1j*sin(pi*nu).*((sinc(nu)).^deg);
        
    case 'sinc'
        y=zeros(size(nu));
        y(abs(nu)<0.5)=1;
        y(abs(nu)==0.5)=0.5;
               
    case {'delta'}
        y=1;
        
    case 'ccd'
        y=(sinc(gamma*nu));
        
    
    otherwise
        error('MATLAB:phi_hat:UndefinedBasis', ...
            'Usage: Undefined basis function');
end
