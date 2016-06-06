function [ nn ]= neuralNet_get_nnSuperSet_indexed(refData, index1, index2)

 nn = refData.nnSuperSet(index1).nnSet(index2);

end
