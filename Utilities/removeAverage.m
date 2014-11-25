function average = removeAverage(data,n)

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