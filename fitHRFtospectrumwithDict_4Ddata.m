function [HRFest,pHRF]=findHRFfromspectrumwithDict_4Ddata(data,DICT,TR,freq_range,pHRF)

if nargin<4
    freq_range=[0.25, 0.35, 0.8];
end

TR_dict=0.1; %TR used for dictionary

[numvox1,numvox2,numvox3,n]=size(data);
[l,p,q,r]=size(DICT);

HRFest=zeros(numvox1,numvox2,numvox3,l-2);

n=n*floor(TR/TR_dict); %adapt dictionary to actually used TR
NFFT=floor(n/2)+1; %number of fast fourier bins


if nargin< 5 || isempty(pHRF)
    disp('calculate Spectral dict')
    pHRF=zeros(NFFT,p,q,r);
for i=1:p
    for j=1:q
        D=padarray(squeeze(DICT(:,i,j,:))./repmat(max(squeeze(DICT(:,i,j,:))),[l,1]),[n-l,0],'post');%padarray(squeeze(DICT(:,i,j,:)),[n-l,0],'post');
        pH=abs(fft(D));
        pHRF(1,i,j,:)=pH(1,:);
        pHRF(2:NFFT,i,j,:)=2*pH(2:NFFT,:);
    end
end
end

%% Define frequency band for fit

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
%% fit spectrum of HRFs to spectrum of data

disp('fit HRF')
for S=1:numvox1
    clear pdata
    data_flat=reshape(squeeze(data(S,:,:,:)),numvox2*numvox3,n); %reshape 3D data to 2D
    [m,~]=size(data_flat);
    
    if m==0
        disp(['region ',num2str(S),' contains no signal'])
        continue
    end
    
    if all(reshape(data_flat,m*n,1)==0) %if slice contains no brain signal skip this slice
        HRFest(S,:,:,:)=zeros(numvox2,numvox3,l-2);
        disp([num2str(round((S)/numvox1*100)),'% done'])
        continue
    end
    
    C=zeros(m,p,q,r);
    HRF_fin=zeros(l,m);
    FAC=zeros(m,p,q,r);
    data_s=zscore(data_flat')';

    pdata=abs(fft(rescale(data_s),[],2)); %calculate spectrum of data
    %pdata=pdata(:,1:NFFT);
    pdata=(pdata)';
    
    %fit HRF dictionary to data spectrum with least squares fit
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

    %check wich HRF fits best
    a=zeros(1,m); b=zeros(1,m);c=zeros(1,m);
    for numreg=1:m
        if squeeze(C(numreg,:,:,:))==0
            %disp('timeseries contains only zeros')
            HRF_fin(:,numreg)=zeros(l,1);
            %numzeros=numreg;
        else
            [a(:,numreg),b(:,numreg),c(:,numreg)] = ind2sub(size(squeeze(C(numreg,:,:,:))),find(squeeze(C(numreg,:,:,:)) == min(min(min(squeeze(C(numreg,:,:,:)))))));
            HRF_fin(:,numreg)=FAC(numreg,a(1,numreg),b(1,numreg),c(1,numreg))*DICT(:,a(1,numreg),b(1,numreg),c(1,numreg));
        end
    end
    
%     if exist('numzeros','var')==1
%     HRF_fin(:,numzeros)=[];
%     data_s(numzeros,:)=[];
%     end
HRFest(S,:,:,:)=reshape(HRF_fin(3:end,:)',numvox2,numvox3,l-2); %delete first 2 time bins as stimulus for HRFs started at time point 3 (i.e. 0.2s)
disp([num2str(round((S)/numvox1*100)),'% done'])
end
