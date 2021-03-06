function [shapeSet All2D] = SetShaperAxes(shapeSet, dat, All2D firstRun)

    
wd_Low  = shapeSet.wd_Low; % bounds of pixels to use from spectrometer
wd_High = shapeSet.wd_High;

 %% Setup Axies

%% stuff for the fist run
if firstRun 
    shapeSet.wnwidth;
    lambdas = dat.allWL;
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


% Zoom in on Excitation frequencies within XX cm^-1 of the transition
try % find the first point that is higher than the low bound
    pixlow = find( (w_e-shapeSet.w_center+wnwidth > 0), 1,'first');
catch
    % if we don't find a point, use the first point
    pixlow = 1;
end 
 
try
    pixhigh = find((w_e-shapeSet.w_center-wnwidth > 0), 1,'first');
catch
   pixhigh = length(w_d); 
end
flim = [pixlow, pixhigh];
welim = w_e(flim(1):flim(2));
w_eAxis = w_e(flim(1):flim(2));

All2D.w_eAxis  = w_eAxis;
All2D.w_dAxis = w_dAxis;
shapeSet.dims = [length(flim(1):flim(2)), length(pixlim(1):pixlim(2))];

else % stuff to run teh first time only 


% trim data to resonable dimentions
data = data(:s);% trim to only needed region
All2D.w_d = All2D.w_d(w_Low:w_High);

end





end % function