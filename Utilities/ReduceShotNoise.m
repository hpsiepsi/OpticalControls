function [ data, intensity, var ] = ReduceShotNoise( data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
sSize = size(data);

% noise floor removal - average of first and last 10 pixels
spec_bg = 1/20*( mean(data( (1:10),:)) + mean(data((end-10):end,:)));

data = data - mean(spec_bg);
clear spec_bg

% integrate allong one axis

intensity = sum(data, 1);
avg_int   = mean(intensity);
std_int   = std(intensity);

 one_mat = avg_int * ones( sSize(1) ,1);
 correctionMat = one_mat *(intensity.^-1);
%correctionMat = (avg_int * ones( sSize(1) ,1)) *(intensity.^-1);

intensity = mean(intensity);
var = std_int/avg_int*100;

data = data .* correctionMat;

%disp(['Intensity variations: ' num2str(std_int/avg_int*100) '%,  Average Intensity : ' num2str(avg_int, '%10.4e')])

end

 