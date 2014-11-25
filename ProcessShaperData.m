%function data = ProcessShaperData(shapeSet)


dispDat.filename  = '';
dispDat.tTime     = 0;
dispDat.tLast     = 0;
dispDat.curNum    = 0;
dispDat.tNum      = 0;
dispDat.noiseLast = 0;
dispDat.ampLast   = 0;
dispDat.gui       = 0;
dispDat.status    = 'Starting';
dispDat = consoleUpdateShaper2D(dispDat, []);

%% Read the calibration file
dat.allWL = readWL( [shapeSet.basepath 'calibration.txt'] );

% setup storage structure
All2D.data  = [];
All2D.intensities = [];
All2D.t2    = [];
All2D.w_e   = [];
All2D.w_d   = [];
All2D.dim   = 0;
All2D.fspec = [];
All2D.fspecComp = 0;
All2D.normalizeSplices = shapeSet.normalizeSplices;


% Range within which falls the Spectra to analyze
wd_Low = 200;
wd_High = 600;

%wd_Low = 450;
%wd_High = 1150;


 %% Setup Axies
wnwidth = shapeSet.wnwidth;
lambdas = dat.allWL(wd_Low:wd_High);

% Spectrometer / detect axis
if shapeSet.CPU_correction
    %w_d = lambdas;  % 12478 is an estimate
    w_d = 10000000*lambdas.^-1 - 12478 + shapeSet.wnOff_d;
else
    w_d = 10000000*lambdas.^-1 - 12478 + shapeSet.wnOff_d;  % 12478 is an estimate 
end

% fft / excitation axis

w_e = 33.3564095*1/(.001*shapeSet.deltaT*2)* ...
       linspace(0,1,shapeSet.fftSize/2)+ shapeSet.rotFrame + shapeSet.wnOff_e; % indata2 THz

try % find the first point that is higher than the low bound
    pixlow = find( (w_d-shapeSet.w_center-wnwidth < 0), 1,'first');
catch %#ok<*CTCH>
    % if we don't find a point, use the first point
    pixlow = 1;
end 
 
if isempty(pixlow)
    pixlow = 1;
end

try
    pixhigh = find((w_d-shapeSet.w_center+wnwidth < 0), 1,'first');
catch
    pixhigh = length(w_d); 
end

if isempty(pixhigh)
    pixhigh = length(w_d);
end

pixlim = [pixlow, pixhigh];

%All2D.w_d  = w_d;
w_dAxis = w_d(pixlim(1):pixlim(2));
wdlim = [min(w_d) max(w_d)];
All2D.w_dAxis = w_dAxis;

% Zoom in on Excitation frequencies within XX cm^-1 of the transition
try % find the first point that is higher than the low bound
    pixlow = find( (w_e-shapeSet.w_center+wnwidth > 0), 1,'first');
catch
    % if we don't find a point, use the first point
    pixlow = 1;
end 
if isempty(pixlow)
    pixlow = 1;
end
 
try
    pixhigh = find((w_e-shapeSet.w_center-wnwidth > 0), 1,'first');
catch
   pixhigh = length(w_d); 
end

if isempty(pixhigh)
    pixhigh = length(w_d);
end

flim = [pixlow, pixhigh];
welim = w_e(flim(1):flim(2));
w_eAxis = w_e(flim(1):flim(2));

All2D.w_eAxis  = w_eAxis;
  
  
 dims = [length(flim(1):flim(2)), length(pixlim(1):pixlim(2))];

%% Process
% Single run, with splicing

dataold = [];
intensities = 0;
savedfiles = {};
spliceCount = size(shapeSet.splices,1);

% pre-automation style splices cell converted to automation style cell
if (size(shapeSet.splices,2) == 3)  
   shapeSet.splices(spliceCount,4) = {'-'};
    for sCount = 1:size(shapeSet.splices,1)
       shapeSet.splices{sCount,4} = [shapeSet.molPrefix sprintf('%d',shapeSet.splices{sCount,1}) ...
           '_' sprintf('%d', floor(shapeSet.deltaT)) ...
        '_' sprintf('%d',shapeSet.splices{sCount,2}) shapeSet.extPostfix];
    end
    clear t1_l t1_h
end
shapeSet.t1_low    = shapeSet.splices{1,1};

shapeSet.t1Splits = zeros(1,spliceCount);
%shapeSet.phaseOffset = (shapeSet.t1_low/1000)*(shapeSet.w_trans/33.333);
%shapeSet.phaseOffset = shapeSet.phaseOffset-floor(shapeSet.phaseOffset);

lastt2 = shapeSet.splices{1,3};

dispDat.tNum   = spliceCount;
dispDat.alltic = tic;
dispDat.ltic = tic;

for sCount = 1:spliceCount   
   
    
    t1_low      = shapeSet.splices{sCount,1};
    t1_high     = shapeSet.splices{sCount,2};
    t2          = shapeSet.splices{sCount,3};
    shapeSet.cur_t2 = t2;
    currentFile = shapeSet.splices{sCount,4};
    
    dispDat.curNum   = sCount;
    dispDat.filename = [shapeSet.input currentFile];
    
    if ~(lastt2 == t2)
        % clear dataold
        dataold = [];
        
        % extra background subtraction
        if shapeSet.remBackground
            data = data - removeAverage(data',shapeSet.remBackground)';
        end
        
        % save file
        %savedfiles(lenght(savedfiles)+1) = Save2Ddata(shapeSet, data); 
        All2D.dim = All2D.dim+1;
        
        dispDat.status   = 'Computing FFT'; dispDat = consoleUpdateShaper2D(dispDat, []);
        
        All2D.t2(All2D.dim)          = lastt2;
        All2D.data(:,:,All2D.dim)    = data;
        All2D.intensities(All2D.dim) = intensities;
        fspecs = fft(data', shapeSet.fftSize);
        fspecs = fspecs(flim(1):flim(2),pixlim(1):pixlim(2));
        
        All2D.zoom2D(1:(dims(1)),1:(dims(2)),All2D.dim) = fspecs;
        clear fspecs
        
        dispDat.tTime    = toc(dispDat.alltic);
        dispDat.tLast    = toc(dispDat.ltic);
        dispDat.status   = ['Finished t_2 = ' num2str(lastt2)];
        dispDat = consoleUpdateShaper2D(dispDat, ...
                   All2D.zoom2D(:,:,All2D.dim), w_eAxis, w_dAxis);
        dispDat.ltic = tic;

        
    end
    
    
    
    % extract extra per-slice data (not currently used)
    if (size(shapeSet.splices,2) == 7)
        shapeSet.expParam1 = shapeSet.splices{sCount,5};
        shapeSet.expParam2 = shapeSet.splices{sCount,6};
        shapeSet.expParam3 = shapeSet.splices{sCount,7};
    end
    
    shapeSet.steps     =  length(t1_low:shapeSet.deltaT:t1_high);
    shapeSet.shotpstep = ( (shapeSet.pcBGdiff == true) +1)*( (shapeSet.pcTAdiff == true)+1);  
    shapeSet.slice     = shapeSet.steps * shapeSet.shotpstep + shapeSet.cleanSpecnumber;
    
    shapeSet.t1Splits(sCount) = shapeSet.steps;
    
    %shapeSet.file         = [shapeSet.molPrefix sprintf('%d',shapeSet.t1_low) '_' sprintf('%d', floor(shapeSet.deltaT)) ...
    %    '_' sprintf('%d',shapeSet.t1_high) shapeSet.extPostfix];
    %shapeSet.pp_savefilename   = ['dataWCO_' sprintf('%d',t1_low) '_' sprintf('%d', deltaT) '.mat'];
%    shapeSet.pp_base           = [shapeSet.basepath shapeSet.input currentFile];

    %% read spectral data from file
    %disp(['Processing ' shapeSet.input currentFile])
    dispDat.status   = 'Reading file and filtering'; dispDat = consoleUpdateShaper2D(dispDat, []);
    
    [dat.allData, dat.sys.NumPixels, sys.NumY, sys.NumFrames, settings.hiCap] =...
        read_winspec([shapeSet.basepath shapeSet.input currentFile],'l',false);
    
    data = reshape(dat.allData,1340,[]);
    
    clear dat  % clear dat to save memory
    
    if shapeSet.verbose
       %try
           figure(8007)
           imagesc( data(250:1050, 1:(2*shapeSet.slice+2)) )
           title(' First Two slices, no processing ')
       %end
    end
    
    % cut data to region of interest and number of full slices
    
    ttshots = max(size(data));
    nSlices = floor(ttshots/shapeSet.slice);
    nSlices = CohereData(shapeSet.reduceingData, 1, nSlices);
    tshots  = shapeSet.slice * nSlices;
    
    data = data(wd_Low:wd_High, 2:(tshots+1));
    
    % Remove dropped frames through interpolation
    if shapeSet.interpDropFrm
       data = interpDroppedFrames(data); 
    end
    
    % Reduce Shot Noise by scaling each spectra
    if shapeSet.rmShotNoise || shapeSet.normalizeSplices
       [data intensities dispDat.noiseLast] = ReduceShotNoise(data);
       
       dispDat.ampLast = intensities;
       
       if shapeSet.normalizeSplices
            data = (1/intensities)*data;
       end % normalize
    end
    
  
    data(:,1:(10)) = 0;  % zero out some of the overflows from the spectrometer. 
                         %   only usefull for displaying pre-processed data

    dispDat.status   = 'Reducing Averages'; dispDat = consoleUpdateShaper2D(dispDat, []);
                         
  if shapeSet.reduceingData>1  
    reducedData = zeros(size(data,1),shapeSet.slice);
    indexArray  = 0:(tshots-1);
    modArray    = mod(indexArray,shapeSet.slice)+1;

    for counter = 1:(shapeSet.slice)
       indexArray = ( counter == modArray);
       reducedData(:,counter) = sum(data(:, indexArray),2)/nSlices;  
    
        if shapeSet.verbose
            if counter == 1
                tempArray = find(indexArray);
                
                disp(['First slice starts at ' int2str(tempArray(1)) ', Second slice starts at ' int2str(tempArray(2))...
                    ', Predicted speration:' int2str(shapeSet.slice)...
                    ' Measured: ' int2str(tempArray(2)-tempArray(1))])
                clear tempArray
            elseif counter == 25 &&( shapeSet.SpecInterferometry)
                % 25 is just a random spectra larger than the typical cleanSpectrometer number
                figure(50);  plot(data(150:800,indexArray));
                title('Waveforms averaging to spectra 25')
                
                figure(51);  plot(data(150:800, (shapeSet.slice-1):(shapeSet.slice+2)));
                title('Last 2 spectra slice n, first 2 slice n+1')
            end
        end % end if verbose
    end 
     
    data = reducedData;
    clear reducedData
    
    if shapeSet.verbose &&( shapeSet.SpecInterferometry)
           % 25 is just a random spectra larger than the typical cleanSpectrometer number
           figure(52);  plot(data(150:800,25));
            title('Waveform of averaged spectra 25')
                
           figure(53);
            plot(data(150:800,19:22));
            title('Waveforms of averaged spectra 19-22')
    end
    if shapeSet.verbose   
        figure(8011);
         imagesc( data )
         title('Averaged Waveforms, pre-corrections')
         addWaveformImageStuff
         
         tmp_size = size(data);
         disp(['data size is: ' int2str( tmp_size(1) ) ' by ' int2str( tmp_size(2) ) ' after averages' ])
         clear tmp_size
    end
   
  else
    % No averaging of slices, use only slice 1
    data = data(:,1:(shapeSet.slice));  
  end
      
  data(:,1:(shapeSet.cleanSpecnumber)) = [];  % remove spectra clean waveforms
  
  %% Corrections %%%%%%%%%%%%%%%%%%%%%%%%%  
  
  % Chirped Pulse Upconversion Correction
  if shapeSet.CPU_correction
      dispDat.status   = 'Correcting for Chirped Pulse'; dispDat = consoleUpdateShaper2D(dispDat, []);
      [data, axis_y] = nutoHz_interp_pp(data', lambdas, shapeSet.wnOff_d);   
	  data = chirp_corr_pp(data, axis_y, 8.32, shapeSet.highcapacity); %8.32
      data = data';
      data = real(data);
  end
  
  if shapeSet.verbose  
    figure(99)
    plot(lambdas, data(:,220:224))
     title('Spectras 220-224')
  end
  
  %% Pump Probe and Phase Cycling subtraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  sSize      = size(data);
  sSize(2)   = floor(sSize(2)/shapeSet.shotpstep);
  shiftCount = 2;
     
     dispDat.status = 'Phase Cycling Subtraction'; dispDat = consoleUpdateShaper2D(dispDat, []);
     s00 = downsample(data',shapeSet.shotpstep);
     
     if shapeSet.pcTAdiff
        data_shift = data(:,shiftCount:(size(data,2)));
        shiftCount = shiftCount+1;
        s0pi = downsample(data_shift',shapeSet.shotpstep);
     else 
        s0pi = zeros(sSize)';
     end
     
     if shapeSet.pcBGdiff
        data_shift = data(:,shiftCount:(size(data,2)));
        shiftCount = shiftCount+1;
        spipi = downsample(data_shift',shapeSet.shotpstep);
     else 
        spipi = zeros(sSize)';
     end
     
     if shapeSet.pcBGdiff && shapeSet.pcTAdiff
        data_shift = data(:,shiftCount:(size(data,2)));
        spi0 = downsample(data_shift',shapeSet.shotpstep);
     else 
        spi0 = zeros(sSize)';
     end

     data = s00' - s0pi' - spi0' + spipi';
     %data = s00' - s0pi' - spi0' + spipi';

   if shapeSet.verbose &&( shapeSet.SpecInterferometry)
       figure(54);  
         plot(lambdas, s00(2,:), lambdas, -s0pi(2,:),lambdas, -spi0(2,:),lambdas, spipi(2,:));
         title('Waveforms of spectra 2 with phase cycling')
   end   
   
   if shapeSet.verbose   
       figure(55);
        plot(lambdas, data(:,2));
        title('Waveform of averaged spectra 2 after phase cycling')
               
       figure(8012);
        imagesc( data )
        title('Averaged Waveforms post process, one slice')
   end
   
   % clean up some memory
   clear s00 s0pi spipi spi0
   clear data_shift ttshots reduceingData shiftCount indexArray modArray

    dataold = [dataold data];
   
   lastt2 = t2;
   data = dataold; % not sure if this is needed
   
end % for number of splice files

if shapeSet.verbose
    for sCount = 2:spliceCount
     shapeSet.t1Splits(sCount) = shapeSet.t1Splits(sCount-1)+shapeSet.t1Splits(sCount);
    end
    shapeSet.t1Splits'; %#ok<VUNUS>
    
    figure(8013);
     imagesc( data )
     title('All Averaged Waveforms Post-process and Splice')
end % verbose

% Store last data set
   % extra background subtraction
    if shapeSet.remBackground
       data = data - removeAverage(data',shapeSet.remBackground)';
    end
  All2D.dim                   = All2D.dim + 1;
  All2D.t2(All2D.dim)         = lastt2;
  All2D.data(:,:,  All2D.dim) = data;
  fspecs = fft(data', shapeSet.fftSize);
  fspecs = fspecs(flim(1):flim(2),pixlim(1):pixlim(2));
        
  All2D.zoom2D(1:(dims(1)),1:(dims(2)),All2D.dim) = fspecs;
  clear fspecs
  
  dispDat.filename = currentFile;
  dispDat.tTime    = toc(dispDat.alltic);
  dispDat.tLast    = toc(dispDat.ltic);
  dispDat.curNum   = sCount;
        
  dispDat = consoleUpdateShaper2D(dispDat, ...
                   All2D.zoom2D(:,:,All2D.dim), w_eAxis, w_dAxis);

%disp('ProcessShaperData: Specs Zoomed')
  
  
% save last file
%shapeSet.savedfiles(size(shapeSet.savedfiles,1)+1,1) = {Save2Ddata(shapeSet, data, 0)};




if ~isunix
  beep;
end

clear counter
disp('Processing Complete')

%eAll2D.fspec(:,:, cntr) = fft(data2', shapeSet.fftSize);nd Function ProcessShaperData()