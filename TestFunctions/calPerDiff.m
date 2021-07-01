function result = calPerDiff(dataA,dataB)
%   Calculate pixelwise Percentage difference between two datasets
result = (dataA - dataB)./dataA*100;
end

