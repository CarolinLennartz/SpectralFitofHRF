function x=deconv_wiener(data,H)

% H needs to be in the form  length(time series)x #timeseries
% data needs to be in the form #timeseries x length(time series)
NFFT=128;
TR=0.1;

% find noise level from data
Fs=1/TR;
f=Fs/NFFT*(0:floor(NFFT/2));

[pad,s]=size(H);

Sf2=abs(fft(data,NFFT,2)).^2;
N0=mean(Sf2(:,find(f>2,1,'first'):floor(NFFT/2)),2);
%Nv=var((Sf2(:,find(f>2,1,'first'):floor(NFFT/2))),[],2);

[m,n]=size(data);
if m<s
    data=data(1:m-1,:);
    [m,n]=size(data);
end

% deconvolve data 
Hf=transpose(fft(padarray(H,[n+2*pad-pad,0],'post')));  %adapt length of HRFs to length of data + 2*length of HRF to adapt for transient effects
Sf_data=fft(padarray(data,[0,pad],'symmetric'),[],2);   %adapt length of data to length of data +2*length of HRF
Sf1=1;
% f=1/(n+2*pad)*(0:floor((n+2*pad-1)));
% Sf1=repmat(sqrt(Nv)*0.001,[1,n+2*pad])./repmat((1+2*(0.3)^2-2*0.3*(1-0.3)*cos(2*pi*f)-2*0.3*cos(4*pi*f)),[m,1]);%(1./f);%AR(2)mit phi1=phi2=0.3
% Sf1(1)=3000;
% Sf1(end)=3000;
Nf1=N0;

% Yfirst=conj(Hf).*Sf1.*Sf_data./(conj(Hf).*Hf.*Sf1+repmat(Nf1,1,n+2*pad));
% x=real(ifft(Yfirst,[],2)); %eigentlich real(ifft(Yfirst,[],2));

x=real(ifft(conj(Hf).*Sf1.*Sf_data./(conj(Hf).*Hf.*Sf1+repmat(Nf1,[1,n+2*pad])),[],2));
x=x(:,pad+1:end-pad);
