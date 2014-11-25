function average = removePolyFit(shapeSet,data,n)
% data(t1-steps, lambdas)

% remove transitions +/- 10 wn - may need to change later

indexWgt = ones(size(data));

for cntr = 1:length(shapeSet.w_trans)

indexWgt = ((abs(shapeSet.w_d - shapeSet.w_trans(cntr))>10).*ones(size(data,1),1) & indexWeight;

average = zeros(size(data));

len = size(data,2);

indexArray = 1:len;

for counter = 1:len
    
    indexPick = (abs(indexArray - counter) <= n);
    
    tsum = sum(data(:,indexPick),2)/sum(indexPick);
    average(:,counter) = tsum;
    
end

%data = data - average;
end