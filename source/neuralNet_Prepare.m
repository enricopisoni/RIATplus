%function [neuralNetStruct]= neuralNet_Prepare(precursors, indicators, x, y, nx, ny, totalCells, refInfo, ...
%    prepareDataInfo, optimizerValues, optimizerCondition)
function [neuralNetStruct]= neuralNet_Prepare( refInfo, commonDataInfo)

%prepare Neural Net

DSuperSet = struct('DSet', {});
nnSuperSet = struct('nnSet', {});
bcSuperSet = struct('bcSet', {});
for h=1:3,
    DSetHorizon = struct('D', {}, 'd', {}, 'Dp', {}, 'dp', {});
    % ET 20140331 - New classes for nnSetHorizon
    nnSetHorizon = struct('ps_input', {}, 'ps_target', {}, 'net', {}, 'icells', {},'Class', {}, 'PRECs', {},'ps_pca',{},'ArPt',{});
    % ET 20140331 - set bcSetHorizon - structure for basecases
    bcSetHorizon= struct('emi_bc',{}, 'conc_bc', {});
    
    for i=1:commonDataInfo.AQINum,
        if (isequal(strtrim(commonDataInfo.pathANN(h).ANNs(i,:)),'-999')==0)
            if commonDataInfo.optim_flags.mode_ce_mo==3 %in case aggregated scenario analysis, do not load Dd
                D=struct('D', {-999}, 'd', {-999}, 'Dp', {-999}, 'dp', {-999});
            else
                D=load(strtrim(commonDataInfo.pathDd(h).Dd(i,:)));
            end
            %nn=load(strtrim((pathANN(h).ANNs(i,:))));
            [nn]=net_read(strtrim((commonDataInfo.pathANN(h).ANNs(i,:)))); %20140617 ET - read model from txt file
            % ET 20140331 load basecases for all nets
            if strcmp(nn.Class,'Delta')==1
                nnp=strtrim((commonDataInfo.pathANN(h).ANNs(i,:))); %net path
                idx1= strfind(nnp,'/'); % modify net path to createBC path
                idx2= strfind(nnp,'_');
                nnp1=nnp(idx2(2):end); %save _TP#
                nnp(idx2(2):end)=[]; %erase _TP#
                nnp(idx1(3)+1:idx2(1))=[]; %erase net_
                bc_concentrations.BC_conc=importdata(strcat(nnp,'_concBC',nnp1,'.txt')); %%20140617 ET - read basecase from txt file
                bc_emissions.BC_emi=importdata(strcat(nnp,'_emiBC',nnp1,'.txt')); %%20140617 ET - read basecase from txt file
                %bc_emissions=load(strcat(nnp,'_emiBC',nnp1));
                %bc_concentrations=load(strcat(nnp,'_concBC',nnp1));
                bc = struct('emi_bc',{bc_emissions.BC_emi}, 'conc_bc', {bc_concentrations.BC_conc});
            else
                bc = struct('emi_bc',{-999}, 'conc_bc', {-999});
            end
        else
            D=struct('D', {-999}, 'd', {-999}, 'Dp', {-999}, 'dp', {-999});
            % ET 20140331 - New classes for nnSetHorizon
            nn = struct('ps_input', {-999}, 'ps_target', {-999}, 'net', {-999}, 'icells', {-999},'Class', {-999}, 'PRECs', {-999},'ps_pca',{-999},'ArPt',{-999});
            bc = struct('emi_bc',{-999}, 'conc_bc', {-999});
        end
        DSetHorizon = [DSetHorizon, D];
        nnSetHorizon= [nnSetHorizon, nn];
        bcSetHorizon= [bcSetHorizon, bc]; % ET 20140331 - New basecase set
    end
    DSuperSet = [DSuperSet, struct('DSet', DSetHorizon)];
    nnSuperSet = [nnSuperSet, struct('nnSet', nnSetHorizon)];
    bcSuperSet = [bcSuperSet, struct('bcSet', bcSetHorizon)]; % ET 20140331 - New basecase Superset
end

% DOptSet and nnOptSet are the data structures to be used in evaluating the
% target AQIs
DOptSet = struct('D', {}, 'd', {}, 'Dp', {}, 'dp', {});
% ET 20140331 - New classes for nnSetHorizon and bcOptSet added
nnOptSet = struct('ps_input', {}, 'ps_target', {}, 'net', {}, 'icells', {},'Class', {}, 'PRECs', {},'ps_pca',{},'ArPt',{});
bcOptSet =struct('emi_bc',{}, 'conc_bc', {});
for o=1:commonDataInfo.optim_flags.optAQINum,
    DOptSet = [DOptSet, DSuperSet(commonDataInfo.aqi_horizon(o)).DSet(commonDataInfo.aqi_obj(o)+1)];
    nnOptSet = [nnOptSet, nnSuperSet(commonDataInfo.aqi_horizon(o)).nnSet(commonDataInfo.aqi_obj(o)+1)];
    bcOptSet = [bcOptSet, bcSuperSet(commonDataInfo.aqi_horizon(o)).bcSet(commonDataInfo.aqi_obj(o)+1)];
end
% data
neuralNetStruct.nnSetHorizon=nnSetHorizon;
neuralNetStruct.nnSuperSet=nnSuperSet;
neuralNetStruct.bcSuperSet=bcSuperSet;
neuralNetStruct.nnOptSet=nnOptSet;
neuralNetStruct.bcOptSet=bcOptSet;

% interface
prefix='neuralNet_';
%neuralNetStruct.getDFunction=strcat(prefix, 'get_D');
neuralNetStruct.get_buildEmissionFunction=strcat(prefix, 'buildEmission');

neuralNetStruct.get_nnSuperSetFunction=strcat(prefix, 'get_nnSuperSet');
neuralNetStruct.get_nnOptSetFunction=strcat(prefix, 'get_nnOptSet');
neuralNetStruct.get_aqipercellFunction=strcat(prefix, 'get_aqipercell');
%neuralNetStruct.get_aqi_per_cellFunction=strcat(prefix, 'get_aqi_per_cell');
%neuralNetStruct.get_do_aqi_per_cellFunction=strcat(prefix, 'do_aqi_per_cell');
neuralNetStruct.get_aqi_per_cellFunction=strcat(prefix, 'do_aqi_per_cell');
neuralNetStruct.get_aqi_bcFunction=strcat(prefix, 'get_aqi_bc');
neuralNetStruct.get_aqi_nnFunction=strcat(prefix, 'get_aqi_nn');
neuralNetStruct.get_aqi_nnOptSetFunction=strcat(prefix, 'get_aqi_nnOptSet');
neuralNetStruct.get_bcSuperSetFunction=strcat(prefix, 'get_bcSuperSet');
neuralNetStruct.get_isDeltaEmisFunction=strcat(prefix, 'get_isDeltaEmis');
neuralNetStruct.get_isAggregatedEmisFunction=strcat(prefix, 'get_isAggregatedEmis');
neuralNetStruct.get_fillcomputesolFunction=strcat(prefix, 'fillcomputesol');
 
neuralNetStruct.is2InputFunction=strcat(prefix, 'is2Input');
neuralNetStruct.is6InputFunction=strcat(prefix, 'is6Input');

%neuralNetStruct.isDelta=strcmp(nn.Class,'Delta')==1;
%neuralNetStruct.isAggregated=strcmp(nn.ArPt,'Aggregated')==1;

end
