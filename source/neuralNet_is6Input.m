function [ flag ] = neuralNet_is6Input(refData, periodIndex, aqiIndex)

 NN=refData.nnSuperSet(periodIndex).nnSet(aqiIndex);
 flag=(size(NN.net,1)==48 || (size(NN.net,1)==1 && size(NN.net.inputs{1}.range,1)==48));
 %flag=(size(NN.net,1)==48 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==48));

end