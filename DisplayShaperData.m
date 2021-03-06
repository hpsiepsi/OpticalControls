% function DisplayShaperData(ShaperSet, All2D)
runall = 0;
shapeSet.plotFFTAll = 0;
fftsize = 2^13;% shapeSet.fftSize;
%shapeSet.contourStyle = 't-shizel';
shapeSet.countourMax = 1;
shapeSet.phaseOffset = 0;

disp('Displaying Shaper 2D data')

 try %#ok<TRYNC>
     if ~(isempty(dispDat))
       close(dispDat.gui);
       clear dispDat
     end
 end

 shapeSet.runPlot = [1];  % which spectra to Plot
% shapeSet.runPlot = [];

w_dAxis = All2D.w_dAxis;
w_eAxis = All2D.w_eAxis;



if  length(shapeSet.runPlot)<1
    shapeSet.runPlot = 1:(All2D.dim);
    runall  = 1;
end

 for cntr = 1:length(shapeSet.runPlot)
 
     runPlot = shapeSet.runPlot(cntr);
     
     if All2D.dim ==1
        cur_t2 = All2D.t2;
        fspec = All2D.zoom2D;
     else
        cur_t2 = All2D.t2(runPlot);
        fspec = All2D.zoom2D(:,:,runPlot);
     end
     
 %% plots of Data    
%  data2  = All2D.data(:,:,runPlot);
%  data2 = data2(pixlim(1):pixlim(2),:);

%  figure(20)
%  plot(data2(200:500,:)')

%  data2 = flipdim(data2,1);
%  lambdas2  = flipdim(lambdas2,1);

%  figure(21)
%  plot(data2(200:500,:)')

%% Fourier Transformed data

% Phase tweak
fspecs = -real(exp(-2*pi*shapeSet.phaseOffset*1i) * fspec);

% Filtering in Time Domain


% Smoothing kernel F to get rid of HF noise
F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
fspecs = conv2(fspecs,F,'same');
%fspecs = conv2(fspecs,F,'same');
%fspecs = conv2(fspecs,F,'same');

% Recovered PP Signal
figure(74)
   plot(w_dAxis, 1/1e6 * real((sum(fspecs,1))));
    title('Excite Sum vs Detection (Recovered Pump-Probe Signal)')

% Set Contour plot levels
fspecsMax = max(max(fspecs));
fspecsMin = min(min(fspecs));

disp(['Spectrum intensity: Max: ' num2str(fspecsMax) '  Min: ' num2str(fspecsMin)])

fspecsMax = max( abs([fspecsMax, fspecsMin]) );
    
% Display Contour Plot
if runall == 0
 hd_fig = figure(1000 + runPlot);
else
 hd_fig = figure(1000);
end

set(hd_fig, 'Name', [shapeSet.input shapeSet.expContFile]);

   switch shapeSet.contourStyle
       case 'zanni'
         cLevels = [ -fspecsMax.*[linspace(shapeSet.countourMax, .07, 12)], ...
             fspecsMax.*linspace(.07, shapeSet.countourMax, 12)];
         contour(w_eAxis, w_dAxis, -fspecs',cLevels, '.', 'LineWidth',2); 
       case 'plain'
         cLevels = [ -fspecsMax.*[linspace(shapeSet.countourMax, .0, 40)], ...
             fspecsMax.*linspace(.0, shapeSet.countourMax, 40)] ;
         contour(w_eAxis, w_dAxis, -fspecs',cLevels);
       case 'filledlines'
         cLevels = [ -fspecsMax.*[linspace(shapeSet.countourMax, .01, 15)], ...
             fspecsMax.*linspace(.01, shapeSet.countourMax, 15)];
         contourf(w_eAxis, w_dAxis, -fspecs',cLevels, '.', 'LineWidth',.5);
       otherwise
         cLevels = [ -fspecsMax.*[linspace(shapeSet.countourMax, 0.03, 25)], ...
             fspecsMax.*linspace(0.03, shapeSet.countourMax, 25)];
         contourf(w_eAxis, w_dAxis, -fspecs'+2e-5,cLevels, 'LineStyle', 'none');         
   end
    %title([ '2D spectra of ' strrep(shapeSet.expContFile, '_', '\_') ',  t_2 = ' num2str(cur_t2) 'ps '], ...
    title([ shapeSet.molPrefix ' in some solvent,  t_2 = ' num2str(cur_t2) 'ps '], ...
        'FontWeight', 'bold', 'FontSize', 16)
    disp(['Contour plot rendered for t_2 = ' num2str(cur_t2) ' ps'])
    %axis([1850 2150 1850 2150])
    axis square;
    
    
    axis(shapeSet.w_center+[-shapeSet.wnwidth shapeSet.wnwidth -shapeSet.wnwidth shapeSet.wnwidth])
    
    min_w = min(w_eAxis);
    max_w = max(w_eAxis);
    
    if shapeSet.plotGuides
    
    for tcntr = 1:length(shapeSet.w_trans)
     % add guidelines to plot
     line([min_w max_w], [shapeSet.w_trans(tcntr) shapeSet.w_trans(tcntr)] ,...
         'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 2)
     line([shapeSet.w_trans(tcntr) shapeSet.w_trans(tcntr)], [min_w max_w], ...
         'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 2)
    end
     line([min_w max_w],[min_w max_w], 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 2)
    end % plot guides
end %for runPlot List


%% additional plots and diagnostics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image plot of entire FFT

     if All2D.dim == 1
        fspec = All2D.data;
     else
        fspec = All2D.data(:,:,1);
     end
   fspecs = fft(fspec', shapeSet.fftSize);

   %fspecs = fspecs((1:floor(shapeSet.fftSize/2)),:);
   % Smoothing kernel F to get rid of HF noise
F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
fspecs = conv2(fspecs,F,'same');
fspecs = conv2(fspecs,F,'same');
fspecs = conv2(fspecs,F,'same');

   figure(43)
     imagesc(abs( fspecs(:,:) )') ;  
     title('Entire FFT of data')
     line([flim(1) flim(2)], [1 1],       'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 2)
     line([flim(1) flim(2)], [1340 1340], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 2)
     line([flim(1) flim(1)], [1 1340],    'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 2)
     line([flim(2) flim(2)], [1 1340],    'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 2)
% end % if plotting all of FFT


    
%figure(832) 
%   contour(w_e(flim(1):flim(2)), w_d(:), fspecs(flim(1):flim(2),:)',20, '.');

%figure(75)
%    plot(w_e, sum(fspecs(1:4096,:),2)); title('Sum vs excitation axis')
    

%figure(75)
%plot(w_e(flim(1):flim(2)),sum(fspecs( flim(1):flim(2),:),2))
%clear w_d w_e t1_low w_eAxis wdlim wnwidth sSize sCount welim
%clear data2 fspecs F tshots pixlow pixhigh