%% This analyzes Pump-Probe geometry spectra. Specificly, with Shaped Pumps

%clear path settings sys
clear
clc

shapeSet.localAnalysisPath = '/Users/mrross/Dropbox/Research/CurrentWorkingAnalysis/';  % path on Mrross Laptop
%shapeSet.localAnalysisPath = '/media/disk/kgroup/DATAANALYSIS/CurrentWorkingAnalysis/';  % DataAnalysis path for kserve

% cd to current analysis, and add all local folders to Matlab path
cd(shapeSet.localAnalysisPath)
addpath(genpath(cd))

%% File SETTINGS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shapeSet.basepath   = '/Users/mrross/Documents/Research/Data/';              % matt's laptop data path
%shapeSet.basepath   = 'C:\Documents and Settings\kgroup\MyDocuments\Data\';  % kdata  data path
%shapeSet.basepath   = '/media/disk/kgroup/Data/';                            % kserve data path

shapeSet.input      = '2012-03-17/';


shapeSet.extPostfix = '.SPE';
shapeSet.file       = '';
%shapeSet.runPlot    = 1;   % plot the nth spectra

shapeSet.expContFile = 'Exp_RDC_PPtests';% 'exp1' will run experiment file 'exp1.m'
                   %                 in 'input' directory. Leave empty to use settings below
%shapeSet.input       = 'Referance/DMDCSlices/';    % referance OK spectra
%shapeSet.expContFile = 'DMDCExp';             % referance OK spectra
%shapeSet.input       = '2012-01-13/';    % referance OK spectra
%shapeSet.expContFile = 'Exp_DMDC_4'; 

%% General Settings  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shapeSet.interpDropFrm    = true;     % Use interpolation to remove dropped frames
shapeSet.remBackground    = 50;       % 0 => no BG removal, otherwise smaller number is more aggressive removal
shapeSet.reduceingData    = 1000;     % Max number of scans to average in a dataset.
shapeSet.rmShotNoise      = true;
shapeSet.normalizeSplices = true;    % Normalize splices to the same probe intensity
shapeSet.divide0fs        = false;    % Divide the first signal by two for better FFTs (Suggested by Zanni group, but I see no improvement)
shapeSet.CPU_correction   = false;    % Correct for the quad phase from CP (not currently implimented)        
shapeSet.highcapacity     = false;    % Default if the data is taken on High Capacity Spectrometer setting
shapeSet.verbose          = false;    % Extra diagnostic output and plots generated

%% Display settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%shapeSet.phaseOffset = 0.0;
shapeSet.fftSize     = 2^13;
shapeSet.wnOff_d     = -90; %The amount to shift the wnaxis
shapeSet.wnOff_e     = 0;
shapeSet.search_w    = true; % increases the plot range to help find the signal
shapeSet.contourStyle = 'fille'; % 'plain', 'zanni', 'filled', 'image'
shapeSet.wd_Low       = 425;      % bounds of pixels to use from spectrometer
shapeSet.wd_High      = 825;
%% Default Experimental Settings (EXP file overrites these)
ExpDefault;

%% Path Settings %%%%%%%%%%%%%%%%%%%%%%%%%% 
shapeSet.inputbasepath  = [shapeSet.basepath shapeSet.input];  % Where the SPE files are saved
shapeSet.outputbasepath = [shapeSet.basepath shapeSet.input];  % Where to save the output.mat
shapeSet.savefilename   = ['All2D' shapeSet.expContFile '.mat']; % Filename for saved pumpprobe specs
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
        shapeSet.w_center   = 2006;
        shapeSet.w_trans(1) = 2006;
        shapeSet.w_trans(2) = 1980;
        shapeSet.w_trans(3) = 2020;
        shapeSet.wnwidth = 70; 
        
elseif strncmp(shapeSet.molPrefix,'RDC_', 3)
        shapeSet.w_center = 2040;
        shapeSet.w_trans  = 2040;
        shapeSet.wnwidth  = 200;
        
elseif strncmp(shapeSet.molPrefix,'FeCO_', 3)
        shapeSet.w_center   = 2000;    
        shapeSet.w_trans(1) = 2000;
        shapeSet.w_trans(1) = 2020;
        shapeSet.wnwidth    = 100;
        
elseif strncmp(shapeSet.molPrefix,'Tpep_', 3)
        shapeSet.w_center = 2063;    
        shapeSet.w_trans  = 2063;
        shapeSet.wnwidth = 50;
        
elseif strncmp(shapeSet.molPrefix,'SpecInf_' , 3)
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

if shapeSet.search_w
    shapeSet.wnwidth = 3*shapeSet.wnwidth;
end
% sort splices in order of acending t2 times
shapeSet.splices = sortCellArray(shapeSet.splices,3);

%%  Process Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('Processing Initated - hold on...')
ProcessShaperData_Jan

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