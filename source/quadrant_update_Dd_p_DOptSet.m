function [res]= quadrant_update_Dd_p_DOptSet(Dp, dp, refData, index1)

 % check if dp or D eq 0 and update only the right one
 refData.DOptSet(index1).Dp = Dp;
 refData.DOptSet(index1).dp = dp;
 res=refData;

end
