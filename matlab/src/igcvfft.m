% IGCVFFT is the inverse counterpart to gcvfft and is designed to 
% enable the reconstruction of the input expansion coefficients 
% (or discrete samples) from the output by a filtering operation 
% equivalent to least-squares computation
% 
% Usage:
% c       = igcvfft(d,s);
% [c,s]   = igcvfft(d,s);
% 
% where
% d is the coefficients characterizing the given continuous signal
% s is a structural array that contains the algorithm parameters
% c is the reconstructed coefficients
% 
% See also: gcvfft
% 
% Nikhil Chacko
% April 27, 2013
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

function [c, s] = igcvfft(d, varargin)


s = gcvfft_parseinputs(varargin);

T1              = s.T1;
h               = s.h;
h_arg           = s.h_arg;
phi1            = s.phi1;
deg1            = s.deg1;
phi2            = s.phi2;
deg2            = s.deg2;
p               = s.p;
q               = s.q;
gamma           = s.gamma;
K               = s.K;

[Nr,Nc]=size(d);

if isvector(T1) && max(size(T1))==2
        Tx=T1(1);
        Ty=T1(2);
    
    elseif isscalar(T1)
        Tx=T1;
        Ty=Tx;
        
    else
        error('MATLAB:igcvfft:SizeOfT1', 'Usage: T1=Tx or [Tx Ty] (Sampling Step)');
end


if ischar(d)
    error('MATLAB:igcvfft:InvalidInput',...
        'Usage:Input coefficients cannot be string');
end


if ~isreal(T1) || Tx<0 || Ty<0
    error('MATLAB:igcvfft:InvalidT1',...
        'Usage:Value of T1(sampling step) should be real and positive');
end


if ~isreal(p) || ~isreal(q) || ~isscalar(p) || ~isscalar(q)
    error('MATLAB:igcvfft:IntMultiRateSampling',...
        'Usage:Values of p and q should be positive integers (T2/T1=p/q)');
elseif round(p)~=p || round(q)~=q || p<=0 || q<=0
    error('MATLAB:igcvfft:IntMultiRateSampling',...
        'Usage:Values of p and q should be positive integers (T2/T1=p/q)');
end


if ~isreal(gamma) || ~isscalar(gamma)
    error('MATLAB:igcvfft:InvalidGamma',...
        'Usage:Value of fill-factor(gamma) should be a real number and 0 < gamma <= 1');
elseif gamma<=0 || gamma>1 
    error('MATLAB:igcvfft:InvalidGamma',...
        'Usage:Value of fill-factor(gamma) should be a real number and 0 < gamma <= 1');
end


if isvector(d)
    if round(max(Nr,Nc)/q)~=max(Nr,Nc)/q
    error('MATLAB:igcvfft:SizeOfInputCoeff',...
        'Usage:Size of given signal should be divisible by q so that the reconstructed size is an integer (T2/T1=q/p)');
    end    
else    
    if round(Nr/q)~=Nr/q || round(Nc/q)~=Nc/q
    error('MATLAB:igcvfft:SizeOfInputCoeff',...
        'Usage:Size of given signal should be divisible by q so that the reconstructed size is an integer (T2/T1=p/q)');
    end
end

if ~isreal(K) || round(K)~=K || K<0 || isinf(K)
    error('MATLAB:igcvfft:InvalidSummationLimit',...
        'Usage:Value of infinite-sum limits should be real and positive');
end


if round(deg1)~=deg1 || deg1<0 || deg1>7
    error('MATLAB:igcvfft:InvalidDeg1', ...
        'Usage: Basis degree should be an integer between 0 and 7');
end

if round(deg2)~=deg2 || deg2<0 || deg2>7
    error('MATLAB:igcvfft:InvalidDeg2', ...
        'Usage: Basis degree should be an integer between 0 and 7');
end

changePhi1toCardinal=0;
if strcmpi(phi1,'cardinal-bspline')
    phi1='bspline';
    changePhi1toCardinal=1;
end


if strcmpi(phi2,'cardinal-bspline')
    phi2='bspline';
    d=changebasis_PeriodicBoundary(d,'Cardinal','B-spline',deg2);
end


if (Nr==1 || Nc==1) && ~strcmpi(h,'shift') %One-dimensional case
    N=length(d)*p/q;
    Nq=N*q;
    nu=(rem((0:Nq-1)+Nq/2,Nq)-Nq/2)/Nq;
    U=zeros(1,Nq);
    
    reShapeNeeded = 0;
    if ~isrow(d)
        d = transpose(d);
        reShapeNeeded = 1;
    end
    
    
    if strcmpi(phi2,'causal-bspline')
        denr=zeros(1,length(nu));
        for n=-K:K
            denr=abs(denr+(exp(-1j*pi.*(nu*p-n)).*sinc(nu*p-n))).^(2*deg2+2);
        end
        
    elseif strcmpi(phi2,'bspline')
        denr=zeros(1,length(nu));
        for n=-K:K
            denr=denr+(sinc(nu*p-n)).^(2*deg2+2);
        end
        
    elseif strcmpi(phi1,'sinc') && strcmpi(phi2,'sinc') && p==1 && q==1
        phi2='delta';
        
    elseif strcmpi(phi2,'CCD') && gamma==1
        phi2='bspline';
        deg2=0;
    end
    
    if strcmpi(phi2,'delta')
        for m=-K:K
            phi1_hat=phi_hat(phi1,(nu-m)*q,deg1,gamma);
            h_hat=filter_hat(h,(nu-m)*q/Tx,h_arg);
            
            U=U+(phi1_hat.*h_hat);
        end
    else
        for m=-K:K           
            phi1_hat=phi_hat(phi1,(nu-m)*q,deg1,gamma);
            phi2_hat=conj(phi_hat(phi2,(nu-m)*p,deg2,gamma));
            h_hat=filter_hat(h,(nu-m)*q/Tx,h_arg);
            
            U=U+(phi1_hat.*phi2_hat.*h_hat);
        end
        
        if ~strcmpi(phi2,'CCD') && ~strcmpi(phi2,'sinc') && deg2~=0
            U = U./denr;
        end
        
    end
    
    if strcmpi(phi1, 'CCD')
        U=q*gamma*U;
    else
        U=q*U;
    end
    
    
    V=zeros(Nq,1);
    for k=0:N*q-1
        Umatrix=zeros(p,q);
        b=[1;zeros(p-1,1)];
        for l=0:p-1
            for m=0:q-1
                Umatrix(l+1,m+1)=U(mod(k+N*m+N*q*l/p,N*q)+1);
            end
        end
        Vvector=pinv(Umatrix)*b;
        V(mod(k:N:k+(q-1)*N,N*q)+1)=Vvector;
    end
    V=transpose(V);

    
    V=V*p*q;
    c=ifft(repmat(fft(d),[1,p]).*V);
    c=downsample(c,q);
    
    
    if changePhi1toCardinal
        c=changebasis_PeriodicBoundary(c,'B-spline','Cardinal',deg1);
    end
    
    if reShapeNeeded
        c = transpose(c);
    end
    
elseif Nr==Nc && Tx==Ty && ~strcmpi(h,'shift') %2-D case (square images with equal sampling steps along x and y direction)
    N=length(d)*p/q;
    Nq=N*q;
    nu=(rem((0:Nq-1)+Nq/2,Nq)-Nq/2)/Nq;
    U=zeros(1,N*q);
    
    if strcmpi(phi2,'causal-bspline')
        denr=zeros(1,length(nu));
        for n=-K:K
            denr=abs(denr+(exp(-1j*pi.*(nu*p-n)).*sinc(nu*p-n))).^(2*deg2+2);
        end
        
    elseif strcmpi(phi2,'bspline')
        denr=zeros(1,length(nu));
        for n=-K:K
            denr=denr+(sinc(nu*p-n)).^(2*deg2+2);
        end
        
    elseif strcmpi(phi1,'sinc') && strcmpi(phi2,'sinc') && p==1 && q==1
        phi2='delta';
        
    elseif strcmpi(phi2,'CCD') && gamma==1
        phi2='bspline';
        deg2=0;
    end
    
    if strcmpi(phi2,'delta')
        for m=-K:K
            phi1_hat=phi_hat(phi1,(nu-m)*q,deg1,gamma);
            h_hat=filter_hat(h,(nu-m)*q/Tx,h_arg);
            
            U=U+(phi1_hat.*h_hat);
        end
        
    else
        for m=-K:K           
            phi1_hat=phi_hat(phi1,(nu-m)*q,deg1,gamma);
            phi2_hat=conj(phi_hat(phi2,(nu-m)*p,deg2,gamma));
            h_hat=filter_hat(h,(nu-m)*q/Tx,h_arg);
            
            U=U+(phi1_hat.*phi2_hat.*h_hat);
        end
        
        if ~strcmpi(phi2,'CCD') && ~strcmpi(phi2,'sinc') && deg2~=0
            U = U./denr;
        end
        
    end
    
    if strcmpi(phi1, 'CCD')
        U=q*gamma*U;
    else
        U=q*U;
    end
    
    V=zeros(Nq,1);
    for k=0:N*q-1
        Umatrix=zeros(p,q);
        b=[1;zeros(p-1,1)];
        for l=0:p-1
            for m=0:q-1
                Umatrix(l+1,m+1)=U(mod(k+N*m+N*q*l/p,N*q)+1);
            end
        end
        Vvector=pinv(Umatrix)*b;
        V(mod(k:N:k+(q-1)*N,N*q)+1)=Vvector;
    end
    V=transpose(V);
    V=V*q*p;
    
    c=ifft2(repmat(fft2(d),[p,p]).*(V.'*V));
    c=downsample(c,q);
    c=downsample(c',q);
    c=c';
    
    if changePhi1toCardinal
        c=changebasis_PeriodicBoundary(c,'B-spline','Cardinal',deg1);
    end
    
else %2-D case (rectangular images)
    
    Nr=Nr*p/q;
    Nc=Nc*p/q;
    
    Nrq=Nr*q;
    Ncq=Nc*q;
    
    reShapeNeeded = 0;
    if isvector(d)
        if ~isrow(d)
            d = transpose(d);
            reShapeNeeded = 1;
            Ncq=Nrq;
        end
        Nrq=1;
    end
    
    nuc=(rem((0:Ncq-1)+Ncq/2,Ncq)-Ncq/2)/Ncq;
    nur=((rem((0:Nrq-1)+Nrq/2,Nrq)-Nrq/2)/Nrq);
    
    Ur=zeros(1,Nrq);
    Uc=zeros(1,Ncq);
    
    if strcmp(phi2,'causal-bspline')
        denrr=zeros(1,length(nur));
        denrc=zeros(1,length(nuc));
        for n=-K:K
            denrr=abs(denrr+(exp(-1j*pi.*(nur*p-n)).*sinc(nur*p-n))).^(2*deg2+2);
            denrc=abs(denrc+(exp(-1j*pi.*(nuc*p-n)).*sinc(nuc*p-n))).^(2*deg2+2);
        end
        
    elseif strcmp(phi2,'bspline')
        denrr=zeros(1,length(nur));
        denrc=zeros(1,length(nuc));
        for n=-K:K
            denrr=denrr+(sinc(nur*p-n)).^(2*deg2+2);
            denrc=denrc+(sinc(nuc*p-n)).^(2*deg2+2);
        end
        
    elseif strcmp(phi2,'sinc') && strcmp(phi2,'sinc') && p==1 && q==1
        phi2='delta';     
    end
    
    
    if strcmpi(phi2,'delta')
        for m=-K:K
            phi1r=phi_hat(phi1,(nur-m)*q,deg1,gamma);
            phi1c=phi_hat(phi1,(nuc-m)*q,deg1,gamma);
            
            if strcmpi(h,'shift')
                if (length(h_arg)>1)
                    hr=filter_hat(h,(nur-m)*q/Ty,h_arg(2));
                else
                    hr=ones(1,Nrq);
                end
                hc=filter_hat(h,(nuc-m)*q/Tx,h_arg(1));
            else
                hr=filter_hat(h,(nur-m)*q/Ty,h_arg);
                hc=filter_hat(h,(nuc-m)*q/Tx,h_arg);
            end
            
            Ur=Ur+(phi1r.*hr);
            Uc=Uc+(phi1c.*hc);
        end
        
    else
        for m=-K:K
            phi1r=phi_hat(phi1,(nur-m)*q,deg1,gamma);
            phi1c=phi_hat(phi1,(nuc-m)*q,deg1,gamma);
            
            phi2r=conj(phi_hat(phi2,(nur-m)*p,deg2,gamma));
            phi2c=conj(phi_hat(phi2,(nuc-m)*p,deg2,gamma));
            
            if strcmpi(h,'shift')
                if (length(h_arg)>1)
                    hr=filter_hat(h,(nur-m)*q/Ty,h_arg(2));
                else
                    hr=ones(1,Nrq);
                end
                hc=filter_hat(h,(nuc-m)*q/Tx,h_arg(1));
            else
                hr=filter_hat(h,(nur-m)*q/Ty,h_arg);
                hc=filter_hat(h,(nuc-m)*q/Tx,h_arg);
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
        Ur=q*gamma*Ur;
        Uc=q*gamma*Uc;
    else
        Ur=q*Ur;
        Uc=q*Uc;
    end
    
    Vr=zeros(Nrq,1);
    Vc=zeros(Ncq,1);
    
    for k=0:Nrq-1
        Umatrixl=zeros(p,q);
        b=[1;zeros(p-1,1)];
        for l=0:p-1
            for m=0:q-1
                Umatrixl(l+1,m+1)=Ur(mod(k+Nr*m+Nr*q*l/p,Nr*q)+1);
            end
        end
        Vvectorl=pinv(Umatrixl)*b;
        Vr(mod(k:Nr:k+(q-1)*Nr,Nr*q)+1)=Vvectorl;
    end
    Vr=transpose(Vr);
    Vr=Vr*p*q;
    
    
    for k=0:Ncq-1
        Umatrixc=zeros(p,q);
        b=[1;zeros(p-1,1)];
        for l=0:p-1
            for m=0:q-1
                Umatrixc(l+1,m+1)=Uc(mod(k+Nc*m+Nc*q*l/p,Nc*q)+1);
            end
        end
        Vectorc=pinv(Umatrixc)*b;
        Vc(mod(k:Nc:k+(q-1)*Nc,Nc*q)+1)=Vectorc;
    end
    Vc=transpose(Vc);
    Vc=Vc*p*q;
    
    
    if ~isvector(d)
        c=ifft2(repmat(fft2(d),[p,p]).*(Vr.'*Vc));
        c=downsample(c,q);
        c=downsample(c',q);
        c=c';
    else
        c=ifft(repmat(fft(d),[1,p]).*Vc);
        c=downsample(c,q);
    end
    
    
    if changePhi1toCardinal
        c=changebasis_PeriodicBoundary(c,'B-spline','Cardinal',deg1);
    end
    
    if reShapeNeeded
        c = transpose(c);
    end
    
end

