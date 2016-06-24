function [commonDataInfo]= load_CommonDataInfo(f1, type)

 % move to prepare CommonDataInfo
 [dirs, info]=init_CommonVar(f1);%, type);
 [pathANN, pathDd, pathBc]=init_AQIDefinition(dirs.pathDEF, info.AQINum, type);
 commonDataInfo.dirs=dirs;
 commonDataInfo.info=info;
 commonDataInfo.pathANN=pathANN;
 commonDataInfo.pathDd=pathDd;
 commonDataInfo.pathBc=pathBc;
 commonDataInfo.AQINum=info.AQINum;
 %[configurationFiles]=init_ConfigurationVar(dirs.pathCONF);
% [pathstr,name,ext] = fileparts(dirInfo.pathFOO); 
%dirInfo.pathCONF=strcat(pathstr, filesep, 'aggregated_configuration_', callTypeTag, '.txt');         % path for scenario analysis file, containing ms-pollutant % emission reductions

 %[type,configurationFiles]=init_ConfigurationVar(f2);

 [optim_flags, cell_threshold_set, tmp_thres_cost, aqi_weights_init, aqi_weights, aqi_horizon, aqi_obj_function_type, aqi_obj, MNum, MPBudget]=init_OptimDefinition(dirs.pathFOO, info.AQINum);

 commonDataInfo.optim_flags=optim_flags;
 commonDataInfo.aqi_horizon=aqi_horizon;
 commonDataInfo.aqi_obj=aqi_obj;

% commonDataInfo.aggType=type;
% commonDataInfo.aggConfigurationFiles=configurationFiles;
 
 commonDataInfo.cell_threshold_set=cell_threshold_set;
 commonDataInfo.tmp_thres_cost=tmp_thres_cost;
 commonDataInfo.aqi_weights_init=aqi_weights_init;
 commonDataInfo.aqi_weights=aqi_weights;
 commonDataInfo.aqi_horizon=aqi_horizon;
 commonDataInfo.aqi_obj_function_type=aqi_obj_function_type;
 commonDataInfo.aqi_obj=aqi_obj;
 commonDataInfo.MNum=MNum;
 %commonDataInfo.optimizerCondition='(optimizerValues==1 | optimizerValues==2)';
 commonDataInfo.MPBudget=MPBudget;
 
 [commonDataInfo]=commonDataInfo;
 
end
