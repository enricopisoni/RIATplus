function [firstguessStruct]=firstguess_Prepare_init(refInfo, common, refData)

 % interface
 prefix='firstguess_';
 %neuralNetStruct.getDFunction=strcat(prefix, 'get_D');
 firstguessStruct.get_buildEmissionFunction=strcat(prefix, 'buildEmission');

 %firstguessStruct.get_DSuperSetFunction=strcat(prefix, 'get_DSuperSet');
 firstguessStruct.get_DSuperSetFunction=strcat(prefix, 'get_aqi_D');
 firstguessStruct.get_DSuperSetFunction_indexed=strcat(prefix, 'get_DSuperSet_indexed');
 firstguessStruct.get_DOptSetFunction=strcat(prefix, 'get_DOptSet');
 firstguessStruct.get_DOptSet_indexedFunction=strcat(prefix, 'get_DOptSet_indexed');
 firstguessStruct.get_Dd_DOptSetFunction=strcat(prefix, 'get_Dd_DOptSet');
 firstguessStruct.get_Dd_p_DOptSetFunction=strcat(prefix, 'get_Dd_p_DOptSet');
 firstguessStruct.get_build_solution_D_indexedFunction=strcat(prefix, 'get_DSuperSet_indexed');
 firstguessStruct.get_rebuildOrderedEmisFunction=strcat(prefix, 'rebuildOrderedEmis');
 firstguessStruct.get_expectedOutputFunction=strcat(prefix, 'get_expectedOutput');
 
 firstguessStruct.update_Dd_p_DOptSetFunction=strcat(prefix, 'update_Dd_DOptSet');
 firstguessStruct.update_Dd_DOptSetFunction=strcat(prefix, 'update_Dd_p_DOptSet');

end
