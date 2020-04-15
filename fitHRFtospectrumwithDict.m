function [HRF_est,pHRF]=fitHRFtospectrumwithDict(data,DICT,TR,freq_range,pHRF)

%data has size number of observations x number of time points
% Dict is a Dictionary of HRFs with size length of timeseries x a x b x c

%define frequency range [d e f] with d the lower limit of the respiratory peak and e
%the upper limit of the respiratory peak, f is the upper border of the
%frequency band chosen for the fit
if nargin<4
    freq_range=[0.25, 0.35, 0.4];
end

TR_dict=0.1;

[l,p,q,r]=size(DICT);

[m,n]=size(data);
n=n*TR/TR_dict;

NFFT=floor(n/2)+1;

HRF_est=zeros(l,m);

% calculate the spectrum dictionary pHRF if not given
if nargin< 5 || isempty(pHRF)
    disp('calculate Spectral dict')
    pHRF=zeros(NFFT,p,q,r);
    for i=1:p
        for j=1:q
            D=padarray(squeeze(DICT(:,i,j,:))./repmat(max(squeeze(DICT(:,i,j,:))),[l,1]),[n-l,0],'post');
            pH=abs(fft(D));
            pHRF(1,i,j,:)=pH(1,:);
            pHRF(2:NFFT,i,j,:)=2*pH(2:NFFT,:);
        end
    end
end

%% Define frequency band for fit

% frequency bins
f=(1/(TR_dict))/n*(0:(NFFT-1));

id1=find(0.005<=f,1);
id2=find(f<=freq_range(1),1,'last');

if freq_range(1)<0.2
    id3=find(f>=freq_range(2),1);
    id4=find(f>=freq_range(3),1);
    id=[id1:id2,id3:id4];
else
    id=id1:id2;
end
%%
disp('fit HRF') % fit the spectum dictionary to the data

    clear pdata
    
    C=zeros(m,p,q,r);
    FAC=zeros(m,p,q,r);
    data_s=zscore(data')';
    
    pdata=abs(fft(rescale(data_s),[],2));
    %pdata=pdata(:,1:NFFT);
    pdata=pdata';
  
    
    for i=1:p
        for j=1:q
            for k=1:r        
                Hf1=pHRF(id,i,j,k);
                H1=(Hf1'*Hf1)\Hf1';
                fac=H1*pdata(id,:);
                C(:,i,j,k)=sum(abs(fac'*Hf1'-pdata(id,:)'),2);
                FAC(:,i,j,k)=fac';
            end
        end
    end

    a=zeros(1,m); b=zeros(1,m);c=zeros(1,m);
    for numreg=1:m
        if squeeze(C(numreg,:,:,:))==0
            disp(['timeseries ', num2str(numreg), ' contains only zeros'])
            HRF_est(:,numreg)=zeros(l,1);
            numzeros=numreg;
        else
            [a(:,numreg),b(:,numreg),c(:,numreg)] = ind2sub(size(squeeze(C(numreg,:,:,:))),find(squeeze(C(numreg,:,:,:)) == min(min(min(squeeze(C(numreg,:,:,:)))))));
            HRF_est(:,numreg)=FAC(numreg,a(1,numreg),b(1,numreg),c(1,numreg))*DICT(:,a(1,numreg),b(1,numreg),c(1,numreg));
        end
    end
    
    if exist('numzeros','var')==1
        HRF_est(:,numzeros)=[];
        %data_s(numzeros,:)=[];
    end
    
HRF_est=HRF_est(3:end,:);
end