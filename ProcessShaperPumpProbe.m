function [All2D shapeSet] = ProcessShaperPumpProbe(shapeSet)


dispDat.filename = '';
dispDat.tTime    = 0;
dispDat.tLast    = 0;
dispDat.curNum   = 0;
dispDat.tNum     = 0;
dispDat.noiseLast = 0;
dispDat.ampLast   = 0;
dispDat.gui = 0;
dispDat.status = 'Starting';
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

wd_Low  = shapeSet.wd_Low; % bounds of pixels to use from spectrometer
wd_High = shapeSet.wd_High;

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

try % find the first point that is higher than the low bound
    pixlow = find( (w_d-shapeSet.w_center-wnwidth < 0), 1,'first');
catch %#ok<*CTCH>
    % if we don't find a point, use the first point
    pixlow = 1;
end 
if isempty(pixlow);
    pixlow = 1;
end

try
    pixhigh = find((w_d-shapeSet.w_center+wnwidth < 0), 1,'first');
catch
    pixhigh = length(w_d); 
end

if isempty(pixhigh);
    pixhigh = length(w_d);
end

pixlim = [pixlow, pixhigh];

%All2D.w_d  = w_d;
w_dAxis = w_d(pixlim(1):pixlim(2));
wdlim = [min(w_d) max(w_d)];

All2D.w_dAxis = w_dAxis;
All2D.w_d = w_d;
  
  
 dims = length(pixlim(1):pixlim(2));

%% Process
% Single run, with splicing

dataold = [];
intensities = 0;
savedfiles = {};
spliceCount = size(shapeSet.splices,1);


shapeSet.t1Splits = zeros(1,spliceCount);

lastt2 = shapeSet.splices{1,3};

dispDat.tNum     = spliceCount;
dispDat.alltic   = tic;
dispDat.ltic     = tic;

for sCount = 1:spliceCount   
   
    t2          = shapeSet.splices{sCount,3};
    shapeSet.cur_t2 = t2;
    currentFile = shapeSet.splices{sCount,4};
    dispDat.filename = currentFile;
    dispDat.curNum   = sCount;
    dispDat.status   = 'Reading files';
    dispDat = consoleUpdateShaper2D(dispDat, []);
    
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
        
        All2D.t2(All2D.dim)          = lastt2;
        All2D.data(:,All2D.dim)      = ppData;
        All2D.intensities(All2D.dim) = intensities;
        All2d.QPTmode                = shapeSet.QPTmode;
        
        dispDat.filename = [shapeSet.input currentFile];
        dispDat.tTime    = toc(dispDat.alltic);
        dispDat.tLast    = toc(dispDat.ltic);
        dispDat.curNum   = sCount;
        
        dispDat = consoleUpdateShaper2D(dispDat,[]);
        figure(996874565)  
         plot(w_d, ppData)
        dispDat.ltic = tic;

        
    end
    
    
    % extract extra per-slice data
    if (size(shapeSet.splices,2) == 7)
        shapeSet.QPTmode   = shapeSet.splices{sCount,5};
        shapeSet.expParam2 = shapeSet.splices{sCount,6};
        shapeSet.expParam3 = shapeSet.splices{sCount,7};
    end
    
    shapeSet.shotpstep = 2;  
    shapeSet.slice     = 2;


    %% read spectral data from file
    %disp(['Processing ' shapeSet.input currentFile])
    [dat.allData, dat.sys.NumPixels, sys.NumY, sys.NumFrames, settings.hiCap] =...
        read_winspec([shapeSet.basepath shapeSet.input currentFile],'l',false);
    
    data = reshape(dat.allData,1340,[]);
    
    clear dat  % clear dat to save memory
    
    if shapeSet.verbose
       %try
           figure(8007)
           imagesc( data(:, 1:(50)) )
           title(' First 50 spectra, no processing ')
       %end
    end
    
    % cut data to region of interest and number of full slices
        data(:,1:(20)) = [];  % zero out some of the overflows from the spectrometer. 
                         %   only usefull for displaying pre-processed data
    
    
    ttshots = max(size(data));
    nSlices = floor(ttshots/shapeSet.slice);
    nSlices = CohereData(shapeSet.reduceingData, 1, nSlices);
    tshots  = shapeSet.slice * nSlices;
    
    data = data(wd_Low:wd_High, 1:(tshots));
    
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
    
  


    dispDat.status   = 'Averaging and Phase Cycling';
    dispDat = consoleUpdateShaper2D(dispDat, []);
    
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
      
 % data(:,1:(shapeSet.cleanSpecnumber)) = [];  % remove spectra clean waveforms
  %% Corrections %%%%%%%%%%%%%%%%%%%%%%%%%  
  
  % Chirped Pulse Upconversion Correction
  if shapeSet.CPU_correction
      dispDat.status   = 'Chirped Pulse Correction';
      dispDat = consoleUpdateShaper2D(dispDat, []);
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
  
  %% Pump Probe subtraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ppData = data(:,1) - data(:,2);
   
   clear data_shift ttshots reduceingData shiftCount indexArray modArray

    dataold = [dataold ppData];
   
   lastt2 = t2;
   ppData = dataold; % not sure if this is needed
   
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

  All2D.dim = All2D.dim+1;
        
  All2D.t2(All2D.dim)          = lastt2;
  All2D.data(:,All2D.dim)      = ppData;
  All2D.intensities(All2D.dim) = intensities;
  All2d.QPTmode                = shapeSet.QPTmode;
 
  dispDat.filename = currentFile;
  dispDat.tTime    = toc(dispDat.alltic);
  dispDat.tLast    = toc(dispDat.ltic);
  dispDat.curNum   = sCount;
  dispDat.status   = 'Processing Complete';
  dispDat = consoleUpdateShaper2D(dispDat, []);

% save last file
%shapeSet.savedfiles(size(shapeSet.savedfiles,1)+1,1) = {Save2Ddata(shapeSet, data, 0)};


if ~isunix
  beep;
end

clear counter
disp('Processing Complete')


%eAll2D.fspec(:,:, cntr) = fft(data2', shapeSet.fftSize);nd Function ProcessShaperData()


end % end of function