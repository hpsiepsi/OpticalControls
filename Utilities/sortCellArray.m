function sA = sortCellArray(cA, colm)

%isolate the column to sort with
sortCol = cA(:,colm);

% retrieve nummerical values  (Strings may not work)
sortNumCol = zeros(size(sortCol));
for count = 1:max(size(sortCol))
    sortNumCol(count) = sortCol{count};
end

% sort numbers
[sortNumCol indexArray] = sort(sortNumCol);

% sort the Cell array
for iCount = 1:size(cA,1)
    for jCount = 1:size(cA,2)
        sA(iCount,jCount) = cA(indexArray(iCount),jCount);
    end
end

end