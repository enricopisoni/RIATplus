function [ flag ] = neuralNet_isAggregatedEmis(refData, periodIndex, aqiIndex)

 %flag1=strcmp(refData.nn.Class,'Delta')==1;
 flag=strcmp(refData.nn.ArPt,'Aggregated')==1;

end