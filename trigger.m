load('D:\FerretData\0168\Preprocessed\0168_LateralVideo_022_20180724\adc_data.mat')
load('D:\FerretData\0168\Preprocessed\0168_LateralVideo_022_20180724\triggerData.mat')

photo = adc_data(4,:);

section = adc_data(:,1:30000*10);
temp = zeros(1,size(section,2));
temp(section(4,:)>3.2) = 1;
photodiode = zeros(1,size(section,2));
stimDuration = 5; % in sec

for i = 1:(length(photodiode)-30000*stimDuration-1) % 5sec photodiode data is 96700
    if photodiode(i)==1 && sum(photodiode(i:(i+30000*stimDuration))) > 96000
        photodiode(i:(i+30000*stimDuration))=1;
    end
end

plot(temp);
plot(photodiode);