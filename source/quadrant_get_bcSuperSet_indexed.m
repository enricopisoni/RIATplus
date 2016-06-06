function [ bc ]= quadrant_get_DSuperSet_indexed(refData, index1, index2)

 bc = refData.bcSuperSet(index1).bcSet(index2);

end
