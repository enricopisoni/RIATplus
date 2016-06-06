function [firstguessStruct]=firstguess_Prepare_main(refInfo, commonData)

 DSuperSet = struct('DSet', {});
 for h=1:3,
    DSetHorizon = struct('D', {}, 'd', {}, 'Dp', {}, 'dp', {});
    for i=1:commonData.AQINum,
        if (isequal(strtrim(commonData.pathANN(h).ANNs(i,:)),'-999')==0)
            if commonData.optim_flags.mode_ce_mo==3 %in case aggregated scenario analysis, do not load Dd
                D=struct('D', {-999}, 'd', {-999}, 'Dp', {-999}, 'dp', {-999});
            else
                D=load(strcat('.', filesep, strtrim(commonData.pathDd(h).Dd(i,:))));
            end
        else
            D=struct('D', {-999}, 'd', {-999}, 'Dp', {-999}, 'dp', {-999});
        end
        DSetHorizon = [DSetHorizon, D];
     end
     DSuperSet = [DSuperSet, struct('DSet', DSetHorizon)];
 end

 % DOptSet and nnOptSet are the data structures to be used in evaluating the
 % target AQIs
 DOptSet = struct('D', {}, 'd', {}, 'Dp', {}, 'dp', {});
 DOptSet1 = struct('D', {}, 'd', {}, 'Dp', {}, 'dp', {});
 optNum=0;

 for o=1:commonData.optim_flags.optAQINum,
     DOptSet = [DOptSet, DSuperSet(commonData.aqi_horizon(o)).DSet(commonData.aqi_obj(o)+1)];
 end
 
%  for o=1:commonData.AQINum
%      if (isequal(strtrim(commonData.pathANN(1).ANNs(o,:)),'-999')==0)
%          DOptSet1 = [DOptSet1, DSuperSet(commonData.aqi_horizon(o)).DSet(commonData.aqi_obj(o)+1)];
%          optNum = optNum+1
%          if commonData.optim_flags.optAQINum == optNum
%              break
%          end
%      end
%  end
 
 % data
 firstguessStruct.DSetHorizon=DSetHorizon;
 firstguessStruct.DSuperSet=DSuperSet;
 firstguessStruct.DOptSet=DOptSet;
 
 % interface
 prefix='firstguess_';
 %neuralNetStruct.getDFunction=strcat(prefix, 'get_D');
 firstguessStruct.get_buildEmissionFunction=strcat(prefix, 'buildEmission');

 firstguessStruct.get_DSuperSetFunction=strcat(prefix, 'get_DSuperSet');
 firstguessStruct.get_DSuperSetFunction_indexed=strcat(prefix, 'get_DSuperSet_indexed');
 firstguessStruct.get_DOptSetFunction=strcat(prefix, 'get_DOptSet');
 firstguessStruct.get_DOptSet_indexedFunction=strcat(prefix, 'get_DOptSet_indexed');
 firstguessStruct.get_Dd_DOptSetFunction=strcat(prefix, 'get_Dd_DOptSet');
 firstguessStruct.get_Dd_p_DOptSetFunction=strcat(prefix, 'get_Dd_p_DOptSet');
 firstguessStruct.get_build_solution_D_indexedFunction=strcat(prefix, 'get_DSuperSet_indexed');
 firstguessStruct.get_aqi_DFunction=strcat(prefix, 'get_aqi_D');
 firstguessStruct.get_aqi_bcFunction=strcat(prefix, 'get_aqi_bc');
 firstguessStruct.get_aqi_nnFunction=strcat(prefix, 'get_aqi_nn');
 firstguessStruct.get_aqi_nnOptSetFunction=strcat(prefix, 'get_aqi_nnOptSet');
 
 % new MM 20160428
 firstguessStruct.get_fillcomputesolFunction=strcat(prefix, 'fillcomputesol');
 firstguessStruct.get_buildEmissionFunction=strcat(prefix, 'buildEmission');

 firstguessStruct.get_rebuildOrderedEmisFunction=strcat(prefix, 'rebuildOrderedEmis');
 firstguessStruct.get_bcSuperSetFunction=strcat(prefix, 'get_bcSuperSet');
 
 firstguessStruct.update_Dd_p_DOptSetFunction=strcat(prefix, 'update_Dd_DOptSet');
 firstguessStruct.update_Dd_DOptSetFunction=strcat(prefix, 'update_Dd_p_DOptSet');

end
