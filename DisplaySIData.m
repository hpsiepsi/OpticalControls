% function DisplayShaperData(ShaperSet = me)

data2 = data - removeAverage(data',40)';
%data2    = data;
lambdas2 = lambdas;
%% plots of Data
%figure(20)
%plot(data2(200:500,:)')

data2 = flipdim(data2,1);
data2 = flipdim(data2,2);
lambdas2  = flipdim(lambdas2,1);

figure(21)
    plot(data2(200:500,:)')
    title( ' Processed intensity of spec pixels as a unction of t_1 ' )

figure(22)
    imagesc(fliplr(data2))
    title('Processed data')
    
    
%% Setup Axies
wnwidth = 25;

% Spectrometer / detect axis
%w_d = 10000000*lambdas2.^-1 - 12478;  % 12478 is an estimate
%pixlim = [find( (w_d-shapeSet.w_trans+wnwidth > 0), 1,'first'), ...
%    find( (w_d-shapeSet.w_trans-wnwidth > 0), 1,'first')];
%pixlim = [1 700]
%data2 = data2(pixlim(1):pixlim(2),:);
%w_d = w_d(pixlim(1):pixlim(2));


%wdlim = [min(w_d) max(w_d)];

% fft/ excitation axis
fftsize = 2^14;
%w_e = 33.3*1/(.001*shapeSet.deltaT*2)*linspace(0,1,fftsize/2)+shapeSet.rotFrame; % indata2 THz
%flim = [1 fftsize/2];

% Zoom in on Excitation frequencies within 75cm^-1 of the transition
%w_trans = 1980;
%flim = [find( (w_e-shapeSet.w_trans+wnwidth > 0), 1,'first'), ...
%    find( (w_e-shapeSet.w_trans-wnwidth > 0), 1,'first')];
%welim = w_e(flim(1):flim(2));
%w_eAxis = w_e(flim(1):flim(2));


%% Fourier Transform
%data = flipdim(data, 2);

fspec = fft(data2, fftsize);
fspecs = real(exp(-2*pi*shapeSet.phaseOffset*1i)*fspec);
%fspecMax = max(max(fspec(flim,:)));
%fspecs = CohereData(fspec,-0.9*fspecMax, 0.9*fspecMax);


%% Display results
%figure(65);
%contourf(w_e(flim(1):flim(2)), w_d(:), fspecs(flim(1):flim(2),:)',10); 
%disp('Contour plot rendered')
%imagesc(fspecs(flim(1):flim(2),:)')
figure(41)
%contourf(fspecs(flim(1):flim(2),:)',10);
imagesc(fftshift(fspecs,1))
%figure(74)
%plot(w_d, sum(fspecs,1))
%figure(75)
%plot(w_e(flim(1):flim(2)),sum(fspecs( flim(1):flim(2),:),2))
  
clear data2 fspecs