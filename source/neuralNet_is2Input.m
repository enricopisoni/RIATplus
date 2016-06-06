function [ flag ] = neuralNet_is2Input(refData, periodIndex, aqiIndex)

 NN=refData.nnSuperSet(periodIndex).nnSet(aqiIndex);
 flag=(size(NN.net,1)==48 || (size(NN.net,1)==1 && size(NN.net.inputs{1}.range,1)==16));
 %flag=(size(NN.net,1)==16 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==16));

end