function data = CohereData(data, lowLimit, highLimit)



indexArray = (data < lowLimit);
data(indexArray) = lowLimit;
indexArray = (data > highLimit);
data(indexArray) = highLimit;

%out = data;

%if data > highLimit
%    out = highLimit;
%end

%if data < lowLimit
%    out = lowLimit;
%end


end