%% This analyzes Pump-Probe geometry spectra. Specificly, with Shaped Pumps

%clear path settings sys
clear
clc

%cd '/Users/mrross/Dropbox/Research/CurrentWorkingAnalysis'  % path on Mrross Laptop
%addpath functions pseudofunctions GA_stuff Pump_Probe %Pump_Probe/functions

%% File SETTINGS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shapeSet.basepath   = '/Users/mrross/Documents/Research/Data/';             % Matt's laptop
%shapeSet.basepath = 'C:\Documents and Settings\kgroup\My Documents\Data\'; % kdata setting
shapeSet.input    = '2012-01-04/';

shapeSet.extPostfix = '.SPE';
shapeSet.file       = '';
%shapeSet.runPlot   = 1;   % plot the nth spectra

shapeSet.expContFile = 'DMDC_a';% 'exp1' will run experiment file 'exp1.m'
                   %                 in 'input' directory. Leave empty to use settings below
%shapeSet.input       = 'Referance/DMDCSlices/';  %'2011-02-25/DMDC/';    % referance OK spectra
%shapeSet.expContFile = 'DMDCExp';             % referance OK spectra
shapeSet.input       = '2012-01-13/';    % referance OK spectra
shapeSet.expContFile = 'Exp_DMDC_4'; 

%% General Settings  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shapeSet.interpDropFrm    = true;     % Use interpolation to remove dropped frames
shapeSet.remBackground    = 120;       % 0 => no BG removal, otherwise smaller number is more aggressive removal
shapeSet.reduceingData    = 1000;     % Max number of scans to average in a dataset.
shapeSet.rmShotNoise      = true;
shapeSet.normalizeSplices = false;    % Normalize splices to the same probe intensity
shapeSet.divide0fs        = false;    % Divide the first signal by two for better FFTs (Suggested by Zanni group, but I see no improvement)
shapeSet.CPU_correction   = false;    % Correct for the quad phase from CP (not currently implimented)        
shapeSet.highcapacity     = false;    % Default if the data is taken on High Capacity Spectrometer setting
shapeSet.verbose          = false;    % Extra diagnostic output and plots generated

%% Display settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%shapeSet.phaseOffset = 0.0;
shapeSet.fftSize     = 2^12;
shapeSet.wnOff_d     = -35; %The amount to shift the wnaxis
shapeSet.wnOff_e     = -20;
shapeSet.contourStyle = 'thick'; % 'plain', 'thick', 'filled', 'image'

shapeSet.countourMax = .99;
shapeSet.plotFFTAll  = 1;
%% Default Experimental Settings (EXP file overrites these)
ExpDefault;

%% Path Settings %%%%%%%%%%%%%%%%%%%%%%%%%% 
shapeSet.inputbasepath  = [shapeSet.basepath shapeSet.input];  % Where the SPE files are saved
shapeSet.outputbasepath = [shapeSet.basepath shapeSet.input];  % Where to save the output.mat
shapeSet.savefilename   = [shapeSet.molPrefix 'e2_All2D.mat']; % Filename for saved pumpprobe specs
%shapeSet.savePostfix    = '.m';
shapeSet.saveSpecs      = false;

%% Load Experiment file 
%    previously called trigger files
%    Settings above may be overridden by the experiment file !!!

cd([shapeSet.basepath shapeSet.input]);
eval(shapeSet.expContFile)
%    Note:  If there is a syntax error in the Experiment file, the error 
%           will show in the above line, not in that file.

%% Automatic Parameters set based on above settings
% First t1 times 
shapeSet.t1_low  = shapeSet.splices{1,1};
shapeSet.t1_high = shapeSet.splices{1,2};

% Display bounds
shapeSet.w_trans = [];
if strncmp(shapeSet.molPrefix, 'WCO_', 3)
        shapeSet.w_center = 1980;
        shapeSet.w_trans  = 1980;
        shapeSet.wnwidth  = 250;
        
elseif strncmp(shapeSet.molPrefix, 'DMDC_', 3)
        shapeSet.w_center   = 2014;
        shapeSet.w_trans(1) = 2014;
        shapeSet.w_trans(2) = 1983;
        shapeSet.w_trans(3) = 2044;
        shapeSet.wnwidth    = 200; 
        
elseif strncmp(shapeSet.molPrefix, 'RDC_', 3)
        shapeSet.w_center = 2052;
        shapeSet.w_trans  = 2015;
        shapeSet.w_trans  = 2084;
        shapeSet.wnwidth  = 275;
        
elseif strncmp(shapeSet.molPrefix, 'FeCO_', 3)
        shapeSet.w_center = 2000;    
        shapeSet.w_trans  = 2000;
        shapeSet.wnwidth  = 100;
        
elseif strncmp(shapeSet.molPrefix, 'FeMono_', 3)
        shapeSet.w_center = 2002;    
        shapeSet.w_trans  = 2002;
        shapeSet.wnwidth  = 40;
        
elseif strncmp(shapeSet.molPrefix, 'Tpep_', 3)
        shapeSet.w_center = 2063;    
        shapeSet.w_trans  = 2063;
        shapeSet.wnwidth  = 50;
        
        
elseif strncmp(shapeSet.molPrefix, 'SpecInf_' , 3)
        shapeSet.SpecInterferometry = true;
        shapeSet.w_center = 2000;
        shapeSet.w_trans  = 2000;
        shapeSet.wnwidth  = 150;
        shapeSet.verbose  = 1;   % assume we want more detail if we have run SpecInf
        
else
        shapeSet.w_center = 2000;
        shapeSet.w_trans  = 2000;
        shapeSet.wnwidth  = 50;
        
end

% sort splices in order of acending t2 times
shapeSet.splices = sortCellArray(shapeSet.splices,3);

%%  Process Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('Processing Initated - hold on...')
ProcessShaperData

%% Display data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~shapeSet.SpecInterferometry
    DisplayShaperData
else
    DisplaySIData
end


%% Save Data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if shapeSet.saveSpecs;
  filename = [shapeSet.outputbasepath shapeSet.savefilename];
  save(filename, 'All2D')
end