function data_corr=correctHRFlatency(data,HRFest,modality)

%data are the voxel timeseries

if nargin<3
    modality='peakshift';
end

switch(modality)
    %% correct data by shifting data using the relative time to peak difference
	case 'peakshift'

        [~,peak]=max(HRFest);
        [reg,l]=size(data);
        
        %define max positive and negative shift
        meanpeak=floor(median(peak));
        maxnegshift=abs(max(peak)-meanpeak); %late HRF gives negative shift of time series to correct
        maxposshift=abs(min(peak)-meanpeak); % early HRF gives positive shift

        data_corr=zeros(reg,l-maxnegshift-maxposshift);
        % correct data for time shift 
        for j=1:reg
            tau=-(peakdel(j)-meanpeak);
            %data_shift=circshift(data{i}(j,:),tau,2);
            data_corr(j,:)=data(j,maxposshift-tau+1:end-(maxnegshift+tau));
        end

    %% correct data using a Wiener filter
    case 'Wiener'
        data_corr=deconv_wiener(data,HRFest);
end
