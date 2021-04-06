% GCVFFT_MIRROR is a MATLAB function designed for HS symmetric input signals
% to enable the discretization of continuous convolution operations using FFTs, 
% without aliasing, even when the signals are not band-limited. The algorithm 
% exploits an underlying discrete relation between two sets of expansion coefficients 
% (or discrete samples) that uniquely characterize the continuous 
% input and convolved output functions, respectively.
% 
% Usage:
% d       = gcvfft_mirror(c,s);
% [d,s]   = gcvfft_mirror(c,s);
% 
% where
% c is the non-redundant half of the HS symmetric input coefficients characterizing the continuous input
% s is a structural array that contains the algorithm parameters
% d is the output coefficients characterizing the convolution result
% 
% See also: gcvfft_mirror_inv
% 
% Nikhil Chacko
% August 10, 2013
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

function [d, s] = gcvfft_mirror(c, varargin)

s = gcvfft_parseinputs(varargin);

T1              = s.T1;
h               = s.h;
h_arg           = s.h_arg;
phi1            = s.phi1;
deg1            = s.deg1;
phi2            = s.phi2;
deg2            = s.deg2;
gamma           = s.gamma;
K               = s.K;


[Nr,Nc]=size(c);

if isvector(s.T1) && max(size(T1))==2
        Tx=T1(1);
        Ty=T1(2);
    
    elseif isscalar(s.T1)
        Tx=T1;
        Ty=Tx;
        
    else
        error('MATLAB:gcvfft_mirror:SizeOfT1', 'Usage: T1=Tx or [Tx Ty] (Sampling Step)');
end

if ischar(c)
    error('MATLAB:gcvfft_mirror:InvalidInput',...
        'Usage:Input coefficients cannot be string');
end


if ~isreal(T1) || Tx<0 || Ty<0
    error('MATLAB:gcvfft_mirror:InvalidT1',...
        'Usage:Value of T1(sampling step) should be real and positive');
end



if ~isreal(gamma) || ~isscalar(gamma)
    error('MATLAB:gcvfft_mirror:InvalidGamma',...
        'Usage:Value of fill-factor(gamma) should be a real number and 0 < gamma <= 1');
elseif gamma<=0 || gamma>1 
    error('MATLAB:gcvfft_mirror:InvalidGamma',...
        'Usage:Value of fill-factor(gamma) should be a real number and 0 < gamma <= 1');
end

    


if ~isreal(K) || round(K)~=K || K<0 || isinf(K)
    error('MATLAB:gcvfft_mirror:InvalidSummationLimit',...
        'Usage:Value of infinite-sum limits should be real and positive');
end


if round(deg1)~=deg1 || deg1<0 || deg1>7
    error('MATLAB:gcvfft_mirror:InvalidDeg1', ...
        'Usage: Basis degree should be an integer between 0 and 7');
end

if round(deg2)~=deg2 || deg2<0 || deg2>7
    error('MATLAB:gcvfft_mirror:InvalidDeg2', ...
        'Usage: Basis degree should be an integer between 0 and 7');
end

if strcmpi(h, 'shift')
    error('MATLAB:gcvfft_mirror:InvalidFilter',...
        'Shift cannot be used with gcvfft_mirror/igcvfft_mirror (convolution kernel should be symmetric). Use gcvfft/igcvfft instead.')
end
    
if strcmpi(phi1,'cardinal-bspline')
    phi1='bspline';
    c=changebasis_HS(c,'Cardinal','B-spline',deg1);
end

changePhi2toCardinal=0;
if strcmpi(phi2,'cardinal-bspline')
    phi2='bspline';
    changePhi2toCardinal=1;
end


if Nr==1 || Nc==1 %One-dimensional case
    N=length(c);
    N = 2*N;
   
    reShapeNeeded = 0;
    if ~isrow(c)
        c = transpose(c);
        reShapeNeeded = 1;
    end
    
    nu=(rem((0:N-1)+N/2,N)-N/2)/N;
    U = zeros(1,N);
    
    
    if strcmpi(phi2,'causal-bspline')
        error('MATLAB:gcvfft_mirror:InvalidPhi',...
            'Basis function should be even in nature for exploiting mirror-symmetry advantages.')
        
    elseif strcmpi(phi2,'bspline')
        denr=zeros(1,length(nu));
        for n=-K:K
            denr=denr+(sinc(nu-n)).^(2*deg2+2);
        end
        
    elseif strcmpi(phi1,'sinc') && strcmpi(phi2,'sinc')
        phi2='delta';
        
    elseif strcmpi(phi2,'CCD') && gamma==1
        phi2='bspline';
        deg2=0;
    end
    

    
    if strcmpi(phi2,'delta')
        for m=-K:K
            phi1_hat=phi_hat(phi1,(nu-m)*1,deg1,gamma);
            h_hat=filter_hat(h,(nu-m)*1/Tx,h_arg);
            
            U=U+(phi1_hat.*h_hat);
        end

    else
        for m=-K:K           
            phi1_hat=phi_hat(phi1,(nu-m)*1,deg1,gamma);
            phi2_hat=conj(phi_hat(phi2,(nu-m)*1,deg2,gamma));
            h_hat=filter_hat(h,(nu-m)*1/Tx,h_arg);
            
            U=U+(phi1_hat.*phi2_hat.*h_hat);
        end
               
        if ~strcmpi(phi2,'CCD') && ~strcmpi(phi2,'sinc') && deg2~=0
            U = U./denr;
        end
        
    end
    
    if strcmpi(phi1, 'CCD')
            U=gamma*U;
    end
        
    
    
    d=ifft_mirror(fft_mirror(c).*U);
    
    if changePhi2toCardinal
        d=changebasis_HS(d,'B-spline','Cardinal',deg2);
    end
    
    if reShapeNeeded
        d = transpose(d);
    end
    
elseif Nr==Nc && Tx==Ty %2-D case (square images with equal sampling steps along x and y direction)
    
    N=2*Nr;
    nu=(rem((0:N-1)+N/2,N)-N/2)/N;
    U=zeros(1,N);
    
    if strcmpi(phi2,'causal-bspline')
        error('MATLAB:gcvfft_mirror:InvalidPhi',...
            'Basis function should be even in nature for exploiting mirror-symmetry advantages.')
        
    elseif strcmpi(phi2,'bspline')
        denr=zeros(1,length(nu));
        for n=-K:K
            denr=denr+(sinc(nu-n)).^(2*deg2+2);
        end
        
    elseif strcmpi(phi1,'sinc') && strcmpi(phi2,'sinc')
        phi2='delta';
        
    elseif strcmpi(phi2,'CCD') && gamma==1
        phi2='bspline';
        deg2=0;
    end
    
    
    if strcmpi(phi2,'delta')
        for m=-K:K
            phi1_hat=phi_hat(phi1,(nu-m)*1,deg1,gamma);
            h_hat=filter_hat(h,(nu-m)*1/Tx,h_arg);
            
            U=U+(phi1_hat.*h_hat);
        end

    else
        for m=-K:K           
            phi1_hat=phi_hat(phi1,(nu-m)*1,deg1,gamma);
            phi2_hat=conj(phi_hat(phi2,(nu-m)*1,deg2,gamma));
            h_hat=filter_hat(h,(nu-m)*1/Tx,h_arg);
            
            U=U+(phi1_hat.*phi2_hat.*h_hat);
        end
        
        if ~strcmpi(phi2,'CCD') && ~strcmpi(phi2,'sinc') && deg2~=0
            U = U./denr;
        end
        
    end
    
    if strcmpi(phi1, 'CCD')
        U=gamma*U;
    end
    
    
    d=ifft2_mirror(fft2_mirror(c).*(U.'*U));
    
    if changePhi2toCardinal
        d=changebasis_HS(d,'B-spline','Cardinal',deg2);
    end
    
    
else %2-D case (rectangular images or square images with unequal sampling steps along x and y directions)
    
    Nr = 2*Nr;
    Nc = 2*Nc;
    
    nuc=(rem((0:Nc-1)+Nc/2,Nc)-Nc/2)/Nc;
    nur=((rem((0:Nr-1)+Nr/2,Nr)-Nr/2)/Nr);
    
    Ur=zeros(1,Nr);
    Uc=zeros(1,Nc);
       
    
    if strcmpi(phi2,'causal-bspline')
        error('MATLAB:gcvfft_mirror:InvalidPhi',...
            'Basis function should be even in nature for exploiting mirror-symmetry advantages.')
        
    elseif strcmpi(phi2,'bspline')
        denrr=zeros(1,length(nur));
        denrc=zeros(1,length(nuc));
        for n=-K:K
            denrr=denrr+(sinc(nur-n)).^(2*deg2+2);
            denrc=denrc+(sinc(nuc-n)).^(2*deg2+2);
        end
        
    elseif strcmpi(phi1,'sinc') && strcmpi(phi2,'sinc')
        phi2='delta';
        
    elseif strcmpi(phi2,'CCD') && gamma==1
        phi2='bspline';
        deg2=0;
        
    end
    
    
    if strcmpi(phi2,'delta')
        for m=-K:K
            phi1r=phi_hat(phi1,(nur-m)*1,deg1,gamma);
            phi1c=phi_hat(phi1,(nuc-m)*1,deg1,gamma);
            
            if strcmpi(h,'shift')
                if (length(h_arg)>1)
                    hr=filter_hat(h,(nur-m)*1/Ty,h_arg(2));
                else
                    hr=ones(1,Nrq);
                end
                hc=filter_hat(h,(nuc-m)*1/Tx,h_arg(1));
            else
                hr=filter_hat(h,(nur-m)*1/Ty,h_arg);
                hc=filter_hat(h,(nuc-m)*1/Tx,h_arg);
            end
            
            Ur=Ur+(phi1r.*hr);
            Uc=Uc+(phi1c.*hc);
        end
        
    else
        for m=-K:K
            phi1r=phi_hat(phi1,(nur-m)*1,deg1,gamma);
            phi1c=phi_hat(phi1,(nuc-m)*1,deg1,gamma);
            
            phi2r=conj(phi_hat(phi2,(nur-m)*1,deg2,gamma));
            phi2c=conj(phi_hat(phi2,(nuc-m)*1,deg2,gamma));
            
            if strcmpi(h,'shift')
                if (length(h_arg)>1)
                    hr=filter_hat(h,(nur-m)*1/Ty,h_arg(2));
                else
                    hr=ones(1,Nrq);
                end
                hc=filter_hat(h,(nuc-m)*1/Tx,h_arg(1));
            else
                hr=filter_hat(h,(nur-m)*1/Ty,h_arg);
                hc=filter_hat(h,(nuc-m)*1/Tx,h_arg);
            end
            
            Ur=Ur+(phi1r.*phi2r.*hr);
            Uc=Uc+(phi1c.*phi2c.*hc);
            
        end
        
        if ~strcmpi(phi2,'CCD') && ~strcmpi(phi2,'sinc') && deg2~=0
            Ur = Ur./denrr;
            Uc = Uc./denrc;
        end
        
    end
    
    if strcmpi(phi1, 'CCD')
        Ur=gamma*Ur;
        Uc=gamma*Uc;
    end

    
    d=ifft2_mirror(fft2_mirror(c).*(Ur.'*Uc));

    
    if changePhi2toCardinal
        d=changebasis_HS(d,'B-spline','Cardinal',deg2);
    end
    
end



