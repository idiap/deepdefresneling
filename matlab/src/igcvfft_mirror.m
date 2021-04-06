% IGCVFFT_MIRROR is the inverse counterpart to gcvfft_mirror and is designed to 
% enable the reconstruction of the input HS symmetric expansion coefficients 
% (or discrete samples) from the output by a filtering operation 
% equivalent to least-squares computation
% 
% Usage:
% c       = igcvfft_mirror(d,s);
% [c,s]   = igcvfft_mirror(d,s);
% 
% where
% d is the HS symmetric coefficients characterizing the given continuous signal
% s is a structural array that contains the algorithm parameters
% c is the reconstructed coefficients
% 
% See also: gcvfft_mirror
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
   
function [c, s] = igcvfft_mirror(d, varargin)

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


[Nr,Nc]=size(d);

if isvector(s.T1) && max(size(T1))==2
        Tx=T1(1);
        Ty=T1(2);
    
    elseif isscalar(s.T1)
        Tx=T1;
        Ty=Tx;
        
    else
        error('MATLAB:igcvfft_mirror:SizeOfT1', 'Usage: T1=Tx or [Tx Ty] (Sampling Step)');
end

if ischar(d)
    error('MATLAB:igcvfft_mirror:InvalidInput',...
        'Usage:Input coefficients cannot be string');
end


if ~isreal(T1) || Tx<0 || Ty<0
    error('MATLAB:igcvfft_mirror:InvalidT1',...
        'Usage:Value of T1(sampling step) should be real and positive');
end



if ~isreal(gamma) || ~isscalar(gamma)
    error('MATLAB:igcvfft_mirror:InvalidGamma',...
        'Usage:Value of fill-factor(gamma) should be a real number and 0 < gamma <= 1');
elseif gamma<=0 || gamma>1 
    error('MATLAB:igcvfft_mirror:InvalidGamma',...
        'Usage:Value of fill-factor(gamma) should be a real number and 0 < gamma <= 1');
end

    


if ~isreal(K) || round(K)~=K || K<0 || isinf(K)
    error('MATLAB:igcvfft_mirror:InvalidSummationLimit',...
        'Usage:Value of infinite-sum limits should be real and positive');
end


if round(deg1)~=deg1 || deg1<0 || deg1>7
    error('MATLAB:igcvfft_mirror:InvalidDeg1', ...
        'Usage: Basis degree should be an integer between 0 and 7');
end

if round(deg2)~=deg2 || deg2<0 || deg2>7
    error('MATLAB:igcvfft_mirror:InvalidDeg2', ...
        'Usage: Basis degree should be an integer between 0 and 7');
end
    
    
changePhi1toCardinal=0;
if strcmpi(phi1,'cardinal-bspline')
    phi1='bspline';
    changePhi1toCardinal=1;
end


if strcmpi(phi2,'cardinal-bspline')
    phi2='bspline';
    d=changebasis_HS(d,'Cardinal','B-spline',deg2);
end


if Nr==1 || Nc==1 %One-dimensional case
    N=length(d);
    N = 2*N;
    
    reShapeNeeded = 0;
    if ~isrow(d)
        d = transpose(d);
        reShapeNeeded = 1;
    end
    
    nu=(rem((0:N-1)+N/2,N)-N/2)/N;
    U = zeros(1,N);
    
    if strcmpi(phi2,'causal-bspline')
        error('MATLAB:igcvfft_mirror:InvalidPhi',...
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
    

    
    for m=-K:K
        if strcmpi(phi1, 'CCD')
            phi1_hat=gamma*phi_hat(phi1,(nu-m),deg1,gamma);
        else
            phi1_hat=phi_hat(phi1,(nu-m),deg1,gamma);
        end
        
        if strcmpi(phi2,'delta')
            phi2_hat=1;
        elseif strcmpi(phi2,'CCD') || strcmpi(phi2,'sinc') || deg2==0
            phi2_hat=conj(phi_hat(phi2,(nu-m),deg2,gamma));
        else
            phi2_hat=conj(phi_hat(phi2,(nu-m),deg2,gamma))./denr;
        end
        h_hat=filter_hat(h,(nu-m)/Tx,h_arg);
        
        
        U=U+(phi1_hat.*phi2_hat.*h_hat);
        
    end
    
    V=zeros(1,N);
    for k=0:N-1
        Umatrix=U(k+1);
        Vvector=pinv(Umatrix);
        V(k+1)=Vvector;
    end
    
    c=ifft_mirror(fft_mirror(d).*V);
    
    if changePhi1toCardinal
        c=changebasis_HS(c,'B-spline','Cardinal',deg1);
    end
    
    if reShapeNeeded
        c = transpose(c);
    end

    
elseif Nr==Nc && Tx==Ty %2-D case (square images with equal sampling steps along x and y direction)
    
    N=2*Nr;
    nu=(rem((0:N-1)+N/2,N)-N/2)/N;
    U=zeros(1,N);
    
    if strcmpi(phi2,'causal-bspline')
        error('MATLAB:igcvfft_mirror:InvalidPhi',...
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
    
    
    for m=-K:K
        if strcmpi(phi1, 'CCD')
            phi1_hat=gamma*phi_hat(phi1,(nu-m),deg1,gamma);
        else
            phi1_hat=phi_hat(phi1,(nu-m),deg1,gamma);
        end
        
        if strcmpi(phi2,'delta')
            phi2_hat=1;
        elseif strcmpi(phi2,'CCD') || strcmpi(phi2,'sinc') || deg2==0
            phi2_hat=conj(phi_hat(phi2,(nu-m),deg2,gamma));
        else
            phi2_hat=conj(phi_hat(phi2,(nu-m),deg2,gamma))./denr;
        end
        h_hat=filter_hat(h,(nu-m)/Tx,h_arg);
        
        
        U=U+(phi1_hat.*phi2_hat.*h_hat);
        
    end
    
    V=zeros(1,N);
    for k=0:N-1
        Umatrix=U(k+1);
        Vvector=pinv(Umatrix);
        V(k+1)=Vvector;
    end
    
    c=ifft2_mirror(fft2_mirror(d).*(V.'*V));
    
    if changePhi1toCardinal
        c=changebasis_HS(c,'B-spline','Cardinal',deg1);
    end
    
    
else %2-D case (rectangular images or square images with unequal sampling steps along x and y directions)
    
    Nr = 2*Nr;
    Nc = 2*Nc;
    
    nuc=(rem((0:Nc-1)+Nc/2,Nc)-Nc/2)/Nc;
    nur=((rem((0:Nr-1)+Nr/2,Nr)-Nr/2)/Nr);
    
    Ur=zeros(1,Nr);
    Uc=zeros(1,Nc);
       
    
    if strcmpi(phi2,'causal-bspline')
        error('MATLAB:igcvfft_mirror:InvalidPhi',...
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
    
    Vr=zeros(1,Nr);
    for k=0:Nr-1
        Umatrix=Ur(k+1);
        Vvector=pinv(Umatrix);
        Vr(k+1)=Vvector;
    end
    
    Vc=zeros(1,Nc);
    for k=0:Nc-1
        Umatrix=Uc(k+1);
        Vvector=pinv(Umatrix);
        Vc(k+1)=Vvector;
    end
    
    
    c=ifft2_mirror(fft2_mirror(d).*(Vr.'*Vc));

    
    if changePhi1toCardinal
        c=changebasis_HS(c,'B-spline','Cardinal',deg1);
    end
    
end



