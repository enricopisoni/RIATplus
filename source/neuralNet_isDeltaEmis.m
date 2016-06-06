function [ flag ] = neuralNet_isDeltaEmis(refData, periodIndex, aqiIndex)

 flag=strcmp(refData.nn.Class,'Delta')==1;
 %flag2=strcmp(refData.nn.ArPt,'Aggregated')==1;

end