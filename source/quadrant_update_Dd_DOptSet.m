function [res]=quadrant_update_Dd_DOptSet(D, d, refData, index1)

 % check if dp or D eq 0 and update only the right one
 refData.DOptSet(index1).D = D;
 refData.DOptSet(index1).d = d;
 res=refData;

end
