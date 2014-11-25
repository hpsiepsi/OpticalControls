 function input = RemoveOutliers(input, n)

%clc
%n=2;
%test input = [ 1 2 1 4 2 1 -1 -3 1 2 1 0 ;
%          1 2 1 4 2 1 -1 -3 1 2 1 0 ;
%          1 2 1 9 2 1 -9 -3 1 2 1 0 ;
%          1 2 1 2 2 1 -1 -1 1 2 1 0 ]
meanIn = mean(input')'*ones(1,size(input,2));
stdIn = std(input')'*ones(1,size(input,2));

indexArray = (abs(input-meanIn)./stdIn > n);

input(indexArray) = 1/2*( input(circshift(indexArray,[0 -1])) +...
    input(circshift(indexArray,[0 1])));

%averaged  = removeAverge(input, 1)
% input(indexArray) = 0;
% input = input + averaged.*indexArray;
 end