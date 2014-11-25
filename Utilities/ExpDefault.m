%  Default Experimental Parameters
% Comments :
% 

shapeSet.t1_low = 0.000000;  % fs
shapeSet.deltaT = 2.000000;  % fs
shapeSet.t1_high = 500.000000;  % fs
shapeSet.pcTAdiff = true;
shapeSet.pcBGdiff = true;
shapeSet.rotFrame =  0.000;  % cm^-1
shapeSet.cleanSpecnumber = 20;
shapeSet.expType = 105;

% Manually added Settings
%shapeSet.wnOff_d     = 80; %The amount to shift the wnaxis
%shapeSet.wnOff_e     = 170;
shapeSet.phaseOffset = .0;     


%shapeSet.splices = {   0 2000    0};
shapeSet.splices = {   0  600   0;% /WCO6_splice_5/
                     602 1200   0};
      %              1202 1800   0;
      %              1802 2400   0;
      %              2402 3000   0;
      %              3002 3600   0;
      %              3602 4200   0;
      %              4202 4800   0;
      %              4802 5400   0};

shapeSet.molPrefix       = 'DMDC_';
shapeSet.PPdiff          = false;
shapeSet.SpecInterferometry = false;
shapeSet.file        = '';
shapeSet.highcapacity     = false;    % Default if the data is taken on High Capacity Spectrometer setting


%% Rare Settings %%%%%%
shapeSet.pixeltobeneg       = 737;   % The pp spec at this pixel should always be negative (if not -1*pp)
ppsettings.PPsign           = -1;
ppsettings.pphighCapacity   = false; % Set to true if data was taken on high capacity
ppsettings.remBackground    = false; % FFT out the background noise
ppsettings.filterpt400noise = false; % Applies filter to remove noise at pseduotime 400