function output = interpDroppedFrames( data)

d_size = size(data);

% find zeros, add colums, select rows below cut off 
rowindex = find( sum( data == 0 ) );

for rows = rowindex

    % make that row = average of nearby
    
    if (rows > 1) && ( rows < d_size(2) )
        data(:,rows) = 1/2*(data(:, rows-1)+ data(:, rows+1));
    
    elseif rows == 1
        data(:, 1) = data(:,2);
        
    elseif rows == d_size(2)
        data(:,rows) = data(:, rows-1);
    end
    
end

output = data;

end

