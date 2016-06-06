function [flags, cell_threshold_set, tmp_thres_cost, ...
   aqi_weights_init, aqi_weights, aqi_horizon, aqi_obj_function_type, aqi_obj, ...
    MNum, MPBudget]=init_OptimDefinition(pathFOO, AQINum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       flag_optim_oth.txt
%                  (optimization configuration)

flag_optim_oth = load(pathFOO);

flags.coordinateFlip       =flag_optim_oth(1); % 0 means no data flipping is required, 1 means data flipping needed (as in RIAT).
flags.reg_net         =flag_optim_oth(2); % DUMMY VARIABLE - NOT READ HERE. IN RIAT, THIS WAS 0 means you consider linear model, 1 neural network model, for all AQIs.
flags.constraints     =flag_optim_oth(3); % flag_constraints=0 means all techs are not replaceable (keep LB as they are).
flags.mode_ce_mo      =flag_optim_oth(4); % 0 means multi-objective, 1 means cost-effectiveness.
flags.paretopoints    =flag_optim_oth(5); % If 0, compute 5 points of pareto curve. If greataer than 0, it is the execution time to be dedicated to pareto points computation.
cell_threshold_set = zeros(AQINum,1);
cell_threshold_set(1,1) = flag_optim_oth(6);
cell_threshold_set(2,1) = flag_optim_oth(7);
cell_threshold_set(3,1) = flag_optim_oth(8);
cell_threshold_set(4,1) = flag_optim_oth(9);
cell_threshold_set(5,1) = flag_optim_oth(10);
cell_threshold_set(6,1) = flag_optim_oth(11);
cell_threshold_set(7,1) = flag_optim_oth(12);
flags.conv_value           =flag_optim_oth(13); % Convergence value for the optimization algorithm.
flags.areal_point          =flag_optim_oth(14); % 0 means areal and point emissions summed, 1 means areal and point emissions separated.
flags.abs_perc             =flag_optim_oth(15); % Specify if the cost constraints expressed in absolute (0) or MRR% values (1). If cost-effectiveness, put 0.
tmp_thres_cost(1)    =flag_optim_oth(16); % Costs (Meuro over CLE) for Pareto curve (CLE and MRR automatically created). If cost-effectiveness, only this value is considered.
tmp_thres_cost(2)    =flag_optim_oth(17); % Costs (Meuro over CLE) for Pareto curve (CLE and MRR automatically created). If cost-effectiveness, this value is not considered.
tmp_thres_cost(3)    =flag_optim_oth(18); % Costs (Meuro over CLE) for Pareto curve (CLE and MRR automatically created). If cost-effectiveness, this value is not considered.
optAQINum            =flag_optim_oth(19); % Number of AQIs to consider during the optimization.

offset = 19; % The number of the last row read.

% Read other infos for optimization, as:
% - number of AQIs to be optimized;
% - type of aggregation;
% - time horizon.
aqi_obj = zeros(optAQINum,1);
for o=1:optAQINum,
    aqi_obj(o) = flag_optim_oth(offset+o);
end
aqi_obj_function_type = zeros(optAQINum,1);
for o=1:optAQINum,
    aqi_obj_function_type(o) = flag_optim_oth(offset+optAQINum+o);
end
aqi_horizon = zeros(optAQINum,1);
for o=1:optAQINum,
    aqi_horizon(o) = flag_optim_oth(offset+optAQINum+optAQINum+o);
end
aqi_weights = zeros(optAQINum,1);
for o=1:optAQINum,
    aqi_weights(o) = flag_optim_oth(offset+optAQINum+optAQINum+optAQINum+o);
end

% If aqi_weights(1)==1, we are in the "fairness" approach, and all the
% weights are automatically recomputed by the system.
% Otherwise the values of the weights do not change.
aqi_weights_init=aqi_weights(1);
if aqi_weights(1)==2
    normalizing_factor = 0;
    alpha = 0.1;
    for i=0:optAQINum-1,
        normalizing_factor = normalizing_factor + alpha^i;
    end
    for i=0:optAQINum-1,
        aqi_weights(i+1) = (alpha^i) / normalizing_factor;
    end
end

% Read the macrosector budget constraints' rhs.
index = offset+optAQINum+optAQINum+optAQINum+optAQINum+1;
MNum = flag_optim_oth(index);
MId = zeros(MNum,1);
MPBudget = zeros(MNum,1);
for m=1:MNum,
    MId(m) = flag_optim_oth(index+m);
end
for m=1:MNum,
    MPBudget(m) = flag_optim_oth((index+MNum)+m);
end

%aqi.weights_init=aqi_weights_init;
%aqi.weights=aqi_weights;
%aqi.horizon=aqi_horizon;
%aqi.obj_function_type=aqi_obj_function_type;
%aqi.obj=aqi_obj;
flags.optAQINum = optAQINum;

end