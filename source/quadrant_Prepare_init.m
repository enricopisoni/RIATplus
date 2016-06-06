function [quadrantStruct]=quadrant_Prepare_init(refInfo, common, refData)

%compute quadrant emissions

%domain dimensions and other variables

quadrantStruct.geomFactor=4;
%quadrantStruct.radius=commonData.radius;
%quadrantStruct.radius=4;

%  DSuperSet = struct('DSet', {});
%  for h=1:3,
%     DSetHorizon = struct('D', {}, 'd', {}, 'Dp', {}, 'dp', {});
%     for i=1:commonData.AQINum,
%         if (isequal(strtrim(commonData.pathANN(h).ANNs(i,:)),'-999')==0)
%             if commonData.optim_flags.mode_ce_mo==3 %in case aggregated scenario analysis, do not load Dd
%                 D=struct('D', {-999}, 'd', {-999}, 'Dp', {-999}, 'dp', {-999});
%             else
%                 D=load(strtrim(commonData.pathDd(h).Dd(i,:)));
%             end
%         else
%             D=struct('D', {-999}, 'd', {-999}, 'Dp', {-999}, 'dp', {-999});
%         end
%         DSetHorizon = [DSetHorizon, D];
%      end
%      DSuperSet = [DSuperSet, struct('DSet', DSetHorizon)];
%  end
% 
%  % DOptSet and nnOptSet are the data structures to be used in evaluating the
%  % target AQIs
%  DOptSet = struct('D', {}, 'd', {}, 'Dp', {}, 'dp', {});
%  for o=1:commonData.optim_flags.optAQINum,
%      DOptSet = [DOptSet, DSuperSet(commonData.aqi_horizon(o)).DSet(commonData.aqi_obj(o)+1)];
%  end
%  
%  % data
%  quadrantStruct.DSetHorizon=DSetHorizon;
%  quadrantStruct.DSuperSet=DSuperSet;
%  quadrantStruct.DOptSet=DOptSet;
 
 % interface
 prefix='quadrant_';
 %neuralNetStruct.getDFunction=strcat(prefix, 'get_D');
 quadrantStruct.get_buildEmissionFunction=strcat(prefix, 'buildEmission');

 quadrantStruct.get_DSuperSetFunction=strcat(prefix, 'get_DSuperSet');
 quadrantStruct.get_DSuperSetFunction_indexed=strcat(prefix, 'get_DSuperSet_indexed');
 quadrantStruct.get_DOptSetFunction=strcat(prefix, 'get_DOptSet');
 quadrantStruct.get_DOptSet_indexedFunction=strcat(prefix, 'get_DOptSet_indexed');
 quadrantStruct.get_Dd_DOptSetFunction=strcat(prefix, 'get_Dd_DOptSet');
 quadrantStruct.get_Dd_p_DOptSetFunction=strcat(prefix, 'get_Dd_p_DOptSet');
 quadrantStruct.get_build_solution_D_indexedFunction=strcat(prefix, 'get_DSuperSet_indexed');
 quadrantStruct.get_rebuildOrderedEmisFunction=strcat(prefix, 'rebuildOrderedEmis');
 quadrantStruct.get_expectedOutputFunction=strcat(prefix, 'get_expectedOutput');
 
 quadrantStruct.update_Dd_p_DOptSetFunction=strcat(prefix, 'update_Dd_DOptSet');
 quadrantStruct.update_Dd_DOptSetFunction=strcat(prefix, 'update_Dd_p_DOptSet');

end
