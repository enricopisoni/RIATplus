function []=old_MAINmain_61c(f1,f2,f3)

%MM 2015
%%Issue(s)
% Done
% _added a line in flag_optim_path file
% _added a folder (configuration)
% _implemented interface approach (you implement physical code for your
% specific solution, but the abstract layer is always the same)
% ToDo Both
% _test code on different input specific input files (needed for each style of run)
% _check means of variable Dd_(p)_DOptSet (neuralNet_get_Dd_DOptSet function)
% _check means of variable Dd_DOptSet (neuralNet_get_Dd_DOptSet function)
% _check common input files (needed for each style of run)
% _check specific input files (needed for specific style of run)
% _check rightness of input->output variables for emissionAggregationType (agg shape->weight to assign)
% _check rightness of input->output variables for emissionModelUseType (each cell or all the same)
% _check rightness of input->output variables for emissionModelType (nn or regression)
% _check rightness of input->output variables for emissionAbsoluteDeltaValues (concentration or emission)
% ToDo Mirko
% _integrate new "client side" version (domain info loading from text file)
% _link (files, variable...) between init and main
% _gaussian weight assignements AND regression to compute single cell value
% _multiple Ring weight assignements AND ?? to compute single cell value

%End

%20130830-EP
% - aggregated scenario analysis: ANNs bounds can be exceeded of -/+ 30%
% (this is required to apply aggregated scenario analysis to Alsace)
% - approximation of results: if results less than 10^-6, put results to 0
% - statusLOg.txt and command windows messages modified to contain, at the
% beginning of the row, "PROGRESSION: "
%TO CREATE MATLAB EXE, USE MATLAB2012B 32BIT, AND DEPLOYTOOL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%main program
%to run program write:
% MAINmain('inputOPERAaspa/flag_optim_path.txt','inputOPERAaspa/flag_optim_path_postproc.txt')

%load configuration files
% f1='inputOPERAaspa/flag_optim_path.txt';
% f2='inputOPERAaspa/flag_optim_path_postproc.txt';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format short

% 20130618
% with aggregated Scenario Analysis, POSTmain is not used
% (maps of emissions and AQIs are generated in MAINmain)

%20130402
% found bug on "detailed approach" for scenario analysis

%20130320
% now each AQI model can be expressed as linear or anns....in the routine
% that run the model, if command manages the two cases (linear or anns) and
% also the case for with 2 inputs (nox and voc, as for O3) and 6 inputs
% (for PM and NOx)

%20130319
% in the routine "UPDATE" two changes have been made:
% - bug fixed, when selecting techs that do not affect emission reductions (there was a NaN problem)
% - when UB-LB<10^(-4) techs are kept fixed to CLE

%20120620 - 1.4.5
%main_compute_aqi_real_objective now computes also aggregation type 2 and 3

% now the philippe thunis relation has been implemented
% y = 4.0399 x - 74.172
% y=number of exceedances
% x=meanPM10
% relation is derived from chimere/traffic station cells, at EU level
% relation is used with aggregation type 3 (see flag_optim_oth) and only for PM10 ANNs

%20120619 - 1.4.5
%GHG
% It is important to add a new file: out_clemfr_ghg.txt
% that contains, from ms-sec-act-tec-low/high (columns 1-5)
% and then other 8 columns with
%RE(CO2, N20, CH4, Fgas) and unEF CO2, N20, CH4, Fgas [kton/actlev]
% then GHG remaining (after reduction) emissions are computed from actlev,
% applying optimal AR computed for air quality

%traffic constraints
%- nuovo file "traffic_duplication" letto in flag_optim_path
%- nel file si specifica se c'� "traffic duplication, e ID sector per HIG,
%EXT, URB
%- poi programma PERL creare file con nuovi vincoli.

%number of daily threshold exceedances
% - waiting for philippe relation between meanPM10 and number of daily threshold exceedances
% - in optimization, we minimize the total number of exceedances over the domain, while at the end
%     we show the total number of cells with exceedances

% %PAD (policy application domain) VS ACD (aqi computation domain)
% added in flag_optim_path the files
% - flag_optim_dom_pop (for population average aqi)
% - flag_aqi_dom (ACD)
% - flag_optim_dom (already available) si for PAD
% in main_computa_[aqi;aqi_real_obj)] now computation of AQI is performed on ACD

% multipollutant is now managed, with
% - fairnessapproach (weights==2 --> normalization, sorting, using predefined weights)
% - user-defined (weights<1 --> normalization, using user-defined weights)
% added in flag_optim_dom these weights!

%20120528 - update in v.1.4.3
% now also working for CREER domain
%
% AQI population average implemented
% file flag_optim_dom_pop.txt is read, to allow for population averaged AQI
% steps:
% - read file population
% - treat if necessary (if RIAT case, reshape and cut)
% - aggregate using population
%
% Now also daily PM10 number of threshold exceedances is managed. This AQI
% is computed using meanPM10 ANN, implementing an aggregation that apply a
% factor (from meanPM10 to dailyPM10numberOfExceedances) and then computes
% how many cells are beyond the PM10 threshold

%20120517 - updates in v1.4.2
% INIT is modified to manage the "border cells". Now the tool read emissions
% emissions directory (usually "./input/emissions/") considering two subdirectories:
% 1) detailed (for emissions detailed as sec-act-tec)
% 2) aggregated (for emissions aggregated for macrosector)
% in case the cells are partly inside and partly outside the optimization
% domain, the emissions are read in both the directories

%read emissions are saved in emi{:,1} (aggregated emissions), and emi{:,2}
%(detailed emissions). For border cells you have both columns of data

%modified routines: MAINload_emi, MAINbase_emi (if cell is completely or
%partly in optimization domainm you sum up detailed and, if available,
%aggregated data); MAINmain and POSTmain now consider that optimization
%domain is defined by flag_optim_dom==1 and flag_optim_dom==2

%for border cells, in "MAINbase_emi" the contribution of the part of the
%cell that is outside the optimization domain, is added to the base_emi_low
%and base_emi_low

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%20120213 - EP
%added in the "flag_optim_path.txt" file the path for the file containing
%the list of techs ("techsList.csv"). This file is used to get code of
%technologies involved in the NH3 constraints creation

%20120213 - EP
%function "MAINlin_constr" has been modified, to create on-the-fly, in the
%code, the NH3 constraints. Now in the input directory you need to provide
%the file "techsList", and then the code create the list of NH3 constraints.
%I.e. the constraints for SA are: sumAR (SA+SA_LNA+SA_CSA) ... > MFR(SA)
%It means that the routine looks for combined techs with SA, and then sum
%up all the combined techs containing SA, to have sumAR < MFR.

%20120210 - EP
% - POSTmain has been "cleaned up", to remove all the useless variables
%   still used
% - "POSTflag_optim_oth" file has been remove% - inserted comments on the "Update" function (to remove not useful techs)
% - in "flag_optim_path" file now also paths for D, d and scenario mode file
%   have to be specified

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OPERA REQUIREMENTS
%20111206 - now the code works also as scenario mode, and also considering
%different domains to evaluate the results
%to use scenario mode, modify the flag in "flag_optim_oth.txt"

%to consider a different domain to evaluate results, modify
%"flag_optim_dom".txt file

%NB: TO RUN THE SCENARIO MODE MODIFYING THE DOMAIN, YOU ALSO NEED
%CORRESPONDING MODIFIED EMISSIONS (MS LEVEL OUTSIDE DOMAIN, MS-S-A-T LEVEL
%INSIDE DOMAIN)
%NB: ALSO YOU NEED TO PROVIDE THE FILE ./input/xscenario.txt(CONTAINING THE
%AR YOU WANT TO APPLY IN THE SCENARIO MODE)

%20111206 replaceable/not replaceable techs, multiobj vs cost-effectiveness cost
%constraints both in abs and perc: all these RIAT capabilities are working

%20111206 WORKING WITH: meanPM10, meanPM25, aot40, somo35

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%d%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - FLAG_USABLE_TECH: subset of technologies to be used
% - FLAG_REPLACEABLE_TECH: flag_constraints=1 contains replaceable techs

%equations to compute emissions:
%E_ijkp = sum_t(Tijk) A_ijk ef_ijkp (1-eff_ijkp) AR_ijkt
%E_p = E_p(initial) - sum_ijk E_ijkp

%equations to compute costs:
%C_ijk = sum_t(ijk) C_ijkt A_ijk AR_ijkt
%C =  sum_ijk C_ijk

%constraints
% LBt < ARt < UBt - inserted with LB, UB
% sum_t(sap) ARt <= 1
% sum_t(sap) ARt >= sum_t(sap) ARcle

%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%to keep track of computing time
tic

%remove warning
%needed otherwise the compile toolbox is not working
warning off MCR:ClassInUseAtExit
warning off MATLAB:ClassInstanceExists
warning off MATLAB:MKDIR:DirectoryExists

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       READING INPUT FILES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       flag_optim_path.txt

commonDataInfo=load_CommonDataInfo(f1,f3);

aggregationInfo=load_AggregationInfo(commonDataInfo.configurationFiles);
%[dirs, info]=init_CommonVar(f1, 'main');

%new MM
%[configurationFiles]=init_ConfigurationVar(dirs.pathCONF);

%aggregationInfo=load_AggregationInfo(configurationFiles);
%new MM end

%                       aqi_definition.txt

% preparing the data structure to store the AQIs configuration
% 7 columns are the 7 AQIs
% 3 structures are: 1)annual, 2)winter, 3)summer

%[pathANN, pathDd]=init_AQIDefinition(dirs.pathDEF, info.AQINum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       flag_optim_oth.txt
%                  (optimization configuration)

%[optim_flags, cell_threshold_set, tmp_thres_cost, aqi_weights_init, aqi_weights, aqi_horizon, aqi_obj_function_type, aqi_obj, MNum, MPBudget]=init_OptimDefinition(dirs.pathFOO, info.AQINum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FROM HERE ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create directory name
mkdir('output')

if (commonDataInfo.optim_flags.mode_ce_mo==0 | commonDataInfo.optim_flags.mode_ce_mo==1 | commonDataInfo.optim_flags.mode_ce_mo==2)
    %load gains DB
    %load file "out_clemfr.txt", containing informations all quadrupes corinair
    %macrosector-sectorIiasa-activityIiasa-tech, and related informations
    % global_data=load('input/out_clemfr.txt');
    global_data=load(commonDataInfo.dirs.pathOCM);
    
    %extract part of gains DB related to emissions for which no reduction
    %technologies are available
    global_data_noc=global_data(find(global_data(:,4)==commonDataInfo.info.nocID),:);
    global_data(find(global_data(:,4)==commonDataInfo.info.nocID),:)=[];
    
    %extract arADSdet and arREF variable...
    % 20130402 - arADSdet is a variable used for the detailed scenario mode...it
    % has to contain AR values between 0 and 1
    arADSdet=global_data(:,24);
    arREF=global_data(:,25);
    global_data(:,25)=[];
    global_data(:,24)=[];
end

%COLUMNS OF OUT_CLEMFR FILE
%macrosector sector activity tecnology lowHighEmissionFlag
%removalEfficiencies (NOx,VOC,NH3,PM10,PM25,SO2)
%unabatedEmissionFactores (NOx,VOC,NH3,PM10,PM25,SO2)
%unitcost arCLE_year arPOT_year flagReplace flagOptim flagNTM_TM arADSdet
%atREFyear

%hl_flag=1 (low), =2(high);
%replace=1 (not replaceable)
%optim=1 (to be used in the optimization)
%ntm_tm=1 (means technical measure....0 means NTM)
%NOTE: AR and RE are in %

%rows are quadruples * 2 (low and high). Quadruples are related to
%GAINS activities, without NOC

%select optimization domain: %1 means cell in optimization domain, 0 means
%outside (2 is border cell, partly inside and partly outside)

[commonDataInfo]=init_Domain(commonDataInfo);
commonDataInfo.special.xutm=commonDataInfo.domainInfo.special.xutm;
commonDataInfo.special.yutm=commonDataInfo.domainInfo.special.yutm;
commonDataInfo.special.global_data=global_data;
%[domainData, domainInfo]=init_Domain(dirs.pathFOD);

%load cells emissions
%emi: emissions of the starting case (2 rows in the case of cells
%outside domain, sector-activity-techs rows if inside domain)
%if the cell is inside the domain, the correspondent file contains both
%emissions (first 6 columns) and activity level (seventh column)
%emi contains two columnes: column 1 is the aggregated data, column 2 the
%detailed data
nx=commonDataInfo.domainInfo.nx;
ny=commonDataInfo.domainInfo.nx;

pathEMI=commonDataInfo.dirs.pathEMI;
if (commonDataInfo.optim_flags.mode_ce_mo==0 | commonDataInfo.optim_flags.mode_ce_mo==1 | commonDataInfo.optim_flags.mode_ce_mo==2)
    pathEMI=strcat(pathEMI,'TP1/');
end
%[emi]=MAINload_emi(commonDataInfo.domainInfo.ncel,pathEMI,commonDataInfo.domainInfo.flag_optim_dom,commonDataInfo.optim_flags.mode_ce_mo);
%UNIBS(ET)20131001 - flag_optim_dom changed with flag_region_dom in input
[ncel, nx, ny]=calcCellNo('latlon', commonDataInfo);

%[emi,flag_ADS_tp]=MAINload_emi(commonDataInfo.domainInfo.ncel,pathEMI,commonDataInfo.domainInfo.flag_region_dom,commonDataInfo.optim_flags.mode_ce_mo);
[emi,flag_ADS_tp]=MAINload_emi(ncel,pathEMI,commonDataInfo.domainInfo.flag_region_dom,commonDataInfo.optim_flags.mode_ce_mo);
commonDataInfo.special.flag_ADS_tp=flag_ADS_tp;
%in case of RIAT, matrices and input are cut and flipped (due to an initial
%problem in the data). Otherwise no change to the data is performed
if commonDataInfo.optim_flags.coordinateFlip==1
    [newDomainData, newDomainInfo, newEmi]=gridCorrectionCase1(commonDataInfo.domainData,commonDataInfo.domainInfo)
    commonDataInfo.special.xutm=newDomainInfo.special.xutm;
    commonDataInfo.special.yutm=newDomainInfo.special.yutm;
    domainData=newDomainData;
    domainInfo=newDomainInfo;
    emi=newEmi;
end

if commonDataInfo.optim_flags.coordinateFlip==2 %ASPA CASE
    [newDomainData, newDomainInfo, newEmi]=gridCorrectionCase2(commonDataInfo.domainData,commonDataInfo.domainInfo)
    commonDataInfo.special.xutm=newDomainInfo.special.xutm;
    commonDataInfo.special.yutm=newDomainInfo.special.yutm;
    domainData=newDomainData;
    domainInfo=newDomainInfo;
    emi=newEmi;
end

%define number of optimization cells: to be used to apply D and d in the
%routine MAINcompute_aqi (1 is for cells inside optim domain, 2 for border
%cells)
%ncelopt=length(find(domainInfo.flag_optim_dom==1 | domainInfo.flag_optim_dom==2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%part in common with INIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load PM10 to PM25 relationship
pm10aveToExceed=load(commonDataInfo.dirs.pathPM);

%output file name
nomeOUTPUT=commonDataInfo.dirs.pathOUF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE LOG FILE WITH STATUS OF THE PROGRAMME
mkdir(commonDataInfo.dirs.pathOUF)
fileStatus=strcat(commonDataInfo.dirs.pathOUF,'/statusLOG.csv');
fidStatus=fopen(fileStatus,'wt');

fileExit=strcat(commonDataInfo.dirs.pathOUF,'/exitLOG.csv');
fidExit=fopen(fileExit,'wt');

%                    DATA STRUCTURES DEFINITION
%      (TO MANAGE ANN, D AND d FOR ALL THE POSSIBLE CONFIGURATIONS)

% DSuperSet(1) and nnSuperSet(1): annual
% DSuperSet(2) and nnSuperSet(2): winter
% DSuperSet(3) and nnSuperSet(3): summer

%update / merge code 6 november (nn as textual)

%one structure with n sections
%emissionAggregationType (geom)
%emissionModelUseType (each cell or all the same)
%emissionModelType (nn or regression)
%emissionAbsoluteDeltaValues (concentration or emission)

geometryIntermediateData=interface_prepare(aggregationInfo.geometryDataInfo, commonDataInfo);
aggregationInfo.geometryIntermediateData=geometryIntermediateData;

% MM 20160420 comment out
% mathIntermediateData=interface_prepare(aggregationInfo.mathDataInfo, commonDataInfo);
% aggregationInfo.mathIntermediateData=mathIntermediateData;

if size(commonDataInfo.domainData.data,2)==6
    commonDataInfo.flag_region_dom=commonDataInfo.domainData.data(:,3); %regional domain
    flag_optim_dom=commonDataInfo.domainData.data(:,4);%optimization domain
    flag_aqi_dom=commonDataInfo.domainData.data(:,5); %aqi computation domain
    pop=commonDataInfo.domainData.data(:,6); %population
else
    commonDataInfo.flag_region_dom=commonDataInfo.domainData.data(:,3);
    flag_optim_dom=commonDataInfo.domainData.data(:,3);%optimization domain
    flag_aqi_dom=commonDataInfo.domainData.data(:,4); %aqi computation domain
    pop=commonDataInfo.domainData.data(:,5); %population
end

% creating solution data structures
solutionSet = struct('BUDGET', {},...
    'SOL_LIST', struct('X', {}, 'X_free', {}, 'AQIs', {}, ...
    'COST', {}, 'BUDGET', {}, ...
    'COSTPERMACROSECTOR', {}, ...
    'BUDGETPERMACROSECTOR', {}));

solutionList = struct('X', {}, 'X_free', {}, 'AQIs', {}, ...
    'COST', {}, 'BUDGET', {}, ...
    'COSTPERMACROSECTOR', {}, ...
    'BUDGETPERMACROSECTOR', {});

%in all cases except aggregated SA
%load or compute base emissions for each cell sum of emissions, for each cell,
%separated for low and high rows extracte from emi (one row per each sec-act) and then summed
%emissions of the part of the border cell outside optimziation domain, is
%summed up to base_emi_low and base_emi_high
%actlev_final=activity level of all sector-activity-techs
%actlev_final_sum=optimization domain wide sum of activity levels
%save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint14_fix', 'solutionSet');

if (commonDataInfo.optim_flags.mode_ce_mo==0 | commonDataInfo.optim_flags.mode_ce_mo==1 | commonDataInfo.optim_flags.mode_ce_mo==2)
    [actlev_final,actlev_final_sum,base_emi_low,base_emi_high,base_emi_low_noc,...
        base_emi_high_noc]=MAINbase_emi(emi,global_data,flag_optim_dom,commonDataInfo.flag_region_dom,global_data_noc);
    
    %TO CHECK EMISSION MAPS USE THIS COMMAND
    % imagesc(flipud(reshape(base_emi_low(:,4),62,95)))
    
    %select replaceable technologies:   flag_repl_tech
    %flag_repl_tech. 0 means you can substitute the tech, 1 you can't
    flag_repl_tech=global_data(:,21);
    
    %select techs to be used: flag_usab_tech: 0=not usable; 1=usable
    flag_usab_tech=global_data(:,22);
    
    % distinguish between technical and non-technical measures:
    % technical measures: flag_tech_nontech = 1
    % non-technical measures: flag_tech_nontech = 0
    flag_tech_nontech=global_data(:,23);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %20101110 - PROCEDURE TO MANAGE REPLACEABLE TECHNOLOGIES AND SUBSET
    %OF TECHNOLOGIES TO BE USED IN OPTIMIZATION
    % If consideing the out_cle_mfr file:
    % - third column (from final one): 0=tech can be substituted; 1=tech cannot be substituted
    % - second column (from final one): 0=tech not usable; 1=tech usable
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ... TO HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %                     READING INPUT FILES - END
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %starting point for optimization
    CLE = global_data(:,19)/100;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                            CONSTRAINTS
    
    % creation of linear constraints concerning the application rate of the
    % technologies
    [A, B] = ...
        MAINcreate_constraints_matrix(global_data,CLE,commonDataInfo.optim_flags.constraints,commonDataInfo.dirs.pathLST);
    
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint1_new','A','B');
    % [A,B] = MAINlin_constr(global_data,CLE,flag_constraints,pathLST);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                        UB & LB CONSTRAINTS
    
    % LB and UB defined between 0 and 1
    LB=global_data(:,19)/100;
    UB=global_data(:,20)/100;
    
    % if there are no usable techs, modify UB to take this into account, and to
    % get that for these technologies UB=LB
    indNotUsabTec=find(flag_usab_tech==0);
    UB(indNotUsabTec)=LB(indNotUsabTec);
    
    % update LB keeping into account that there are replaceable technologies.
    if commonDataInfo.optim_flags.constraints==1
        LB(find(flag_repl_tech==0))=0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           COST CONSTRAINT
    costconstr = MAINcreate_cost_constraint();
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   COST CONSTRAINT PER MACROSECTOR
    
    % compute the minimum cost per macrosector, i.e. the cost for each
    % macrosector when technologies are applied at the CLE application rate
    minimumCostPerMacrosector = computeMinimumCostPerMacrosector();
    
    % for each macrosector build the mask (costPerMacrosectorMask) to identify
    % the coefficients associated with the technologies which apply to the
    % given macrosector
    
    [num_of_technologies, foo] = size(flag_usab_tech);
    costPerMacrosectorMask = zeros(num_of_technologies,commonDataInfo.MNum);
    for m=1:commonDataInfo.MNum,
        for t=1:num_of_technologies,
            if(global_data(t,1) == MId(m))
                costPerMacrosectorMask(t,m) = 1;
            end
        end;
    end
    
    % for each macrosector derive from the  costconstr's coefficients
    % the coefficients associated with the technologies which apply to
    % the given macrosector (costPerMacrosector)
    
    % costPerMacrosector is used to derive the constraint
    % (costconstrPerMacrosectorLHS) associated with the given macrosector
    
    % costPerMacrosector is also used to compute the cost for each
    % macrosector when technologies are pplied at a given application
    % rate (see buildSolution function)
    
    costPerMacrosector = zeros(num_of_technologies,commonDataInfo.MNum);
    for m=1:commonDataInfo.MNum,
        costPerMacrosector(:,m) = ...
            costPerMacrosectorMask(:,m) .* costconstr;
    end;
    
    % for each macrosector define the constraint associated  with the given
    % macrosector
    %
    % Ex.
    % If only x1 and x2 apply to a given macrosector, and for this macrosector
    % the internal cost increment must be at most 20% of the total internal
    % cost increment (where the increment is computed with respect to internal
    % costs imposed by the CLE), then we have:
    % c1 x1 + c2 x2 - (c1 x1^CLE + c2 x2^CLE) <=
    %   0.2 ((c1 x1 + c2 x2 + ... + cN xN) - (c1 x1^CLE + c2 x2^CLE + ... + cN xN^CLE))
    % and thus:
    % 0.8c1 x1 + 0,.8c2 x2 - 0.2c3 x3 ... - 0.2cN xN <=
    %   (c1 x1^CLE + c2 x2^CLE) - 0.2(c1 x1^CLE + c2 x2^CLE + ... + cN xN^CLE)
    
    % lhs term
    costconstrPerMacrosectorLHS = zeros(num_of_technologies,commonDataInfo.MNum);
    for m=1:commonDataInfo.MNum,
        costconstrPerMacrosectorLHS(:,m) = ...
            costPerMacrosector(:,m) - costconstr * MPBudget(m);
    end;
    
    % rhs term
    costconstrPerMacrosectorRHS = zeros(commonDataInfo.MNum,1);
    for m=1:commonDataInfo.MNum,
        costconstrPerMacrosectorRHS(m,1) = ...
            minimumCostPerMacrosector(MId(m)) - ...
            MPBudget(m) * (sum(minimumCostPerMacrosector));
    end
    
    if(commonDataInfo.MNum > 0)
        A = [A; costconstrPerMacrosectorLHS'];
        B = [B; costconstrPerMacrosectorRHS];
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       TRAFFIC CONSTRAINTS
    
    % read configuration file
    flag_optim_oth = load(commonDataInfo.dirs.pathTRD);
    flagTRD = flag_optim_oth(1); %0 means no traffic duplication, 1 means duplication
    if flagTRD==1
        hig     = flag_optim_oth(2); %hig sector ID
        urb     = flag_optim_oth(3); %urb sector ID
        ext     = flag_optim_oth(4); %ext sector ID
    end
    
    % flagTRD=0;
    if flagTRD==1
        
        trafficConstraintsMatrix = ...
            generateTrafficConstraintsMatrix(hig, urb, ext);
        
        [trafficConstraintsRowNums, trafficConstraintsColNums] = ...
            size(trafficConstraintsMatrix);
        if(trafficConstraintsRowNums > 0)
            for r=1:trafficConstraintsRowNums,
                ind1 = find(trafficConstraintsMatrix(r,:)==-1);
                ind2 = find(trafficConstraintsMatrix(r,:)==1);
                [foo,num1] = size(ind1);
                [foo,num2] = size(ind2);
                if((num1 == 1) && ...
                        (num2 == 1))
                    A = [A; trafficConstraintsMatrix(r,:)];
                    B = [B; 0];
                else
                    if(~((num1 == 0) && (num2 == 0)))
                        fprintf('Some traffic constraints cannot be imposed...\n');
                        fprintf('A check for consistency concerning the technologies'' DB is required...\n');
                    end
                end
                
            end;
        end;
        
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % the cost constraint holds only while computing intermediate pareto curve
    % points
    A1= [A; costconstr'];
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint2_new','A1');
    
    %                         CONSTRAINTS - END
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       DATA STRUCTURES UPDATING
    % update data to be considered for single or multi-pollutant optimization
    [num_of_technologies, foo] = size(flag_usab_tech);
    deltaCost = 0;
    
    % to be defined in the Update() function
    % inhibited_x is the set of techs not to be optimized
    inhibited_x = zeros(num_of_technologies,1);
    
    % to be defined in the Update() function
    % flag to define if the tech is or is not used in optimization
    % (1 means it is kept fixed)
    fixed_x = zeros(num_of_technologies,1);
    
    % 20120210 - EP - this function cut the technologies that does not have
    % impact on the problem dimension.
    % original version without parameter
    %Update();
    % MM version
    nAggregationInfo=Update_(aggregationInfo, commonDataInfo);
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint3_new','A','A1','costconstr');
    aggregationInfo=0;
    aggregationInfo=nAggregationInfo;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           OPTIMIZATION
    % 'active-set', 'trust-region-reflective', 'interior-point', 'levenberg-marquardt','sqp'.
    % opt=optimset('LargeScale','off',...'Algorithm','sqp',...
    % NOTE: SQP is the best!!!
    opt=optimset('LargeScale','off','Algorithm','sqp',...
        'Display','Iter',...
        'Diagnostic','off',...
        'FunValCheck','off',...
        'TolFun',commonDataInfo.optim_flags.conv_value);
    
    % compute CLE values
    x_CLE_free=CLE;
    
    %ENR20130415
    strStatus='PROGRESSION: computing CLE results...';
    disp(strStatus);
    fprintf(fidStatus, '%s\n',strStatus);
    
    % original version
    %sol_CLE = buildSolution(x_CLE_free, fixed_x, inhibited_x, ...
    %    costconstr'*(x_CLE_free) + deltaCost, ...
    %    costconstr'*(x_CLE_free) + deltaCost, ...
    %    );
    % end
    % MM Version (two extra parameters)
    sol_CLE = buildSolution_(x_CLE_free, fixed_x, inhibited_x, ...
        costconstr'*(x_CLE_free) + deltaCost, ...
        costconstr'*(x_CLE_free) + deltaCost, ...
        aggregationInfo, commonDataInfo);
    % MM Version
    solutionList = [solutionList, sol_CLE];
    
    solutionSet = [solutionSet, ...
        struct('BUDGET', sol_CLE.BUDGET, 'SOL_LIST', solutionList)];
        
end
%20130424 END THE PART SPECIFIC TO MO, CE, DETAILED SA

%save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint14_new', 'solutionSet');
% if multi-objective
if commonDataInfo.optim_flags.mode_ce_mo==0       %multi-objective
    
    % fix the starting point for the optimization
    x0=CLE;
    
    %ENR20130415
    strStatus='PROGRESSION: computing MFR optimal results...';
    disp(strStatus);
    fprintf(fidStatus, '%s\n',strStatus);
    
    % compute best solution
    %original code
    %x_BEST_free = MAINcompute_sol(A,B,LB,UB,x0);
    %MM new (aggregationInfo and period parameter added
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint15a_new', 'A','B','LB','UB','x0');
    fName=strtrim(commonDataInfo.pathANN(1).ANNs(1,:));
    aggregationInfo.firstguess=0;
    if (strcmp(fName,'-999') == 0)
         [alpha, omega, radius, flatWeight, pollutantList]=FG_read(fName);
         aggregationInfo.firstguess.alpha=alpha;
         aggregationInfo.firstguess.omega=omega;
         aggregationInfo.firstguess.radius=radius;
         aggregationInfo.firstguess.flatWeight=flatWeight;
         aggregationInfo.firstguess.pollutantList=pollutantList;
    end
    x_BEST_free = MAINcompute_sol_(A,B,LB,UB,x0,aggregationInfo, commonDataInfo, 1);
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint24_new', 'x_BEST_free');
    %end MM
    
    %Original version
    %     sol_BEST = buildSolution(x_BEST_free, fixed_x, inhibited_x, ...
    %         costconstr'*(x_BEST_free) + deltaCost, ...
    %         costconstr'*(x_BEST_free) + deltaCost);
    % MM Version (one extra parameter)
    sol_BEST = buildSolution_(x_BEST_free, fixed_x, inhibited_x, ...
        costconstr'*(x_BEST_free) + deltaCost, ...
        costconstr'*(x_BEST_free) + deltaCost,...
        aggregationInfo, commonDataInfo);
    
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint25_new', 'sol_BEST');
    %     startingcputime=cputime;
    
    % if (flag_paretopoints > 0) you manage the case with a constraint on
    % the total computation time, and not the generation of a fixed number
    % of pareto curve points
    if(commonDataInfo.optim_flags.paretopoints == 0)
        % 3 points are added to the pareto curve
        
        % always use CLE as starting point, at the moment
        x0=CLE;
        
        % pareto curve points definition
        for i = 1:3
            
            % different threshold (budget) if costs are in absolute or
            % percentage values
            if commonDataInfo.optim_flags.abs_perc==0
                thres_cost = sol_CLE.COST + commonDataInfo.tmp_thres_cost(i);
            else
                thres_cost=(sol_CLE.COST + ...
                    (sol_BEST.COST - sol_CLE.COST) * commonDataInfo.tmp_thres_cost(i));
            end
            
            % insert the cost constraint
            B1 = [B;thres_cost - deltaCost];
            
            %ENR20130415
            strStatus = sprintf('PROGRESSION: computing point %d optimal results...',i+1);
            disp(strStatus)
            fprintf(fidStatus, '%s\n',strStatus);
            
            % compute best solution
            % MM Version (two extra parameters)
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint15b_new', 'A1','B1','LB','UB','x0');
            x_sol_free = MAINcompute_sol_(A1,B1,LB,UB,x0,aggregationInfo, commonDataInfo,i);
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint26a_new', 'x_sol_free');
            
            %             sol = buildSolution(x_sol_free, fixed_x, inhibited_x, ...
            %                 costconstr'*(x_sol_free) + deltaCost, ...
            %                 thres_cost);
            % MM Version (one extra parameter)
            sol = buildSolution_(x_sol_free, fixed_x, inhibited_x, ...
                costconstr'*(x_sol_free) + deltaCost, ...
                thres_cost, aggregationInfo, commonDataInfo);
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint26b_new', 'sol');
            
            solutionList = struct('X', {}, 'X_free', {}, 'AQIs', {}, ...
                'COST', {}, 'BUDGET', {}, ...
                'COSTPERMACROSECTOR', {}, ...
                'BUDGETPERMACROSECTOR', {});
            
            solutionList = [solutionList, sol];
            
            solutionSet = [solutionSet, ...
                struct('BUDGET', thres_cost, 'SOL_LIST', solutionList)];
            
        end;
        
        solutionList = struct('X', {}, 'X_free', {}, 'AQIs', {}, ...
            'COST', {}, 'BUDGET', {}, ...
            'COSTPERMACROSECTOR', {}, ...
            'BUDGETPERMACROSECTOR', {});
        
        solutionList = [solutionList, sol_BEST];
        
        solutionSet = [solutionSet, ...
            struct('BUDGET', sol_BEST.BUDGET, 'SOL_LIST', solutionList)];
        
    end;
    
    % in case you want to create a pareto curve with a constraint on
    % optimization time (fixed in maxTime)
    if(commonDataInfo.optim_flags.paretopoints > 0 )
        
        time = toc;
        % computation time to be dedicated to pareto curve
        maxTime=commonDataInfo.optim_flags.paretopoints;
        
        solutionList = struct('X', {}, 'X_free', {}, 'AQIs', {}, ...
            'COST', {}, 'BUDGET', {}, ...
            'COSTPERMACROSECTOR', {}, ...
            'BUDGETPERMACROSECTOR', {});
        
        solutionList = [solutionList, sol_BEST];
        
        % possibily find multiple solutions to insert in
        % solutionList by applying multistart
        % ...
        % possibily find multiple solutions to insert in
        % solutionList by applying multistart - end
        
        solutionSet = [solutionSet, ...
            struct('BUDGET', sol_BEST.BUDGET, 'SOL_LIST', solutionList)];
        
        while(time <= maxTime)
            
            generatedSolutionSet = struct('BUDGET', {},...
                'SOL_LIST', struct('X', {}, 'X_free', {}, 'AQIs', {}, ...
                'COST', {}, 'BUDGET', {}, ...
                'COSTPERMACROSECTOR', {}, ...
                'BUDGETPERMACROSECTOR', {}));
            
            [foo, setSize] = size(solutionSet);
            for i = 1:setSize - 1,
                
                time = toc;
                if (time > maxTime) break; end;
                
                thres_cost = solutionSet(i).BUDGET + ...
                    ((solutionSet(i+1).BUDGET - solutionSet(i).BUDGET)/2);
                
                B1 = [B;thres_cost - deltaCost];
                
                solutionList = struct('X', {}, 'X_free', {}, 'AQIs', {}, ...
                    'COST', {}, 'BUDGET', {}, ...
                    'COSTPERMACROSECTOR', {}, ...
                    'BUDGETPERMACROSECTOR', {});
                
                
                % possibily find multiple solutions to insert in
                % solutionList by applying multistart
                
                % fix the starting point for the optimization
                x0 = solutionSet(i).SOL_LIST(1).X_free + ...
                    ((solutionSet(i+1).SOL_LIST(1).X_free - ...
                    solutionSet(i).SOL_LIST(1).X_free)/2) * 1/2;
                
                indexes = find(LB > x0);
                x0(indexes) = LB(indexes);
                
                % compute solution
                %ENR20130415
                strStatus = sprintf('PROGRESSION: computing optimal results, time (seconds) is at %d...',time);
                disp(strStatus);
                fprintf(fidStatus, '%s\n',strStatus);
                
                %x_sol_free = MAINcompute_sol(A1,B1,LB,UB,x0);
                % MM Version (one extra parameter)
                %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint15c_new', 'A1','B1','LB','UB','x0');
                x_sol_free = MAINcompute_sol_(A1,B1,LB,UB,x0,aggregationInfo, commonDataInfo,1);
                %                 sol = buildSolution(x_sol_free, fixed_x, inhibited_x, ...
                %                     costconstr'*(x_sol_free) + deltaCost, ...
                %                     thres_cost);
                % MM Version (one extra parameter)
                %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint27a_new', 'x_sol_free');
                sol = buildSolution_(x_sol_free, fixed_x, inhibited_x, ...
                    costconstr'*(x_sol_free) + deltaCost, ...
                    thres_cost, aggregationInfo, commonDataInfo);
                %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint27b_new', 'sol');
                
                solutionList = [solutionList, sol];
                
                % possibily find multiple solutions to insert in
                % solutionList by applying multistart - end
                
                
                generatedSolutionSet = [generatedSolutionSet, ...
                    struct('BUDGET', thres_cost, 'SOL_LIST', solutionList)];
                
            end;
            
            
            % insert the new solution in solutionSet keeping it ordered
            % input: solutionSet = [S1, S2, S3, ..., Sn]
            % output: solutionSet = [S1, newS1, S2, newS2, S3, ..., newSn-1, Sn]
            [foo, genSetSize] = size(generatedSolutionSet);
            
            [foo, setSize] = size(solutionSet);
            solutionSetIndex = 1;
            for i=1:genSetSize,
                solutionSet = [solutionSet(1:solutionSetIndex),...
                    generatedSolutionSet(i:i),...
                    solutionSet(solutionSetIndex+1:setSize)];
                solutionSetIndex = solutionSetIndex + 2;
                [foo, setSize] = size(solutionSet);
            end;
            % insert the new solution in solutionSet keeping it ordered - end
            time = toc;
            
        end;
        
    end;
    
elseif commonDataInfo.optim_flags.mode_ce_mo==1   %cost-effectiveness
    
    % fix the starting point for the optimization
    x0=CLE;
    
    % insert the cost constraint
    thres_cost = sol_CLE.COST + commonDataInfo.tmp_thres_cost(1);
    B1 = [B; thres_cost - deltaCost];
    
    strStatus='PROGRESSION: computing the optimal cost-effectiveness results...';
    disp(strStatus);
    fprintf(fidStatus, '%s\n',strStatus);
    
    %compute solution
    %x_sol_free = MAINcompute_sol(A1,B1,LB,UB,x0);
    % MM Version (one extra parameter)
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint15d_new', 'A1','B1','LB','UB','x0');
    x_sol_free = MAINcompute_sol_(A1,B1,LB,UB,x0, aggregationInfo, commonDataInfo,1);
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint28a_new', 'x_sol_free');
    
    %     sol = buildSolution(x_sol_free, fixed_x,
    %     inhibited_x, ...
    %         costconstr'*(x_sol_free) + deltaCost, ...
    %         thres_cost);
    % MM Version (one extra parameter)
    sol = buildSolution_(x_sol_free, fixed_x, ...
        inhibited_x, ...
        costconstr'*(x_sol_free) + deltaCost, ...
        thres_cost, aggregationInfo, commonDataInfo);
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint28b_new', 'sol');
    
    solutionList = struct('X', {}, 'X_free', {}, 'AQIs', {}, ...
        'COST', {}, 'BUDGET', {}, ...
        'COSTPERMACROSECTOR', {}, ...
        'BUDGETPERMACROSECTOR', {});
    
    solutionList = [solutionList, sol];
    
    solutionSet = [solutionSet, ...
        struct('BUDGET', thres_cost, 'SOL_LIST', solutionList)];
    
elseif commonDataInfo.optim_flags.mode_ce_mo==2   % scenario mode
    
    % 20130402 - arADSdet is a variable used for the detailed scenario mode...it
    % has to contain AR values between 0 and 1
    x_sol=arADSdet/100;
    % in case of scenario mode, use arADSdet AR, from out_clemfr (last) column
    %routine "update" is not applied to data
    fixed_x_scenariomode = zeros(num_of_technologies,1);
    originalcostconstr = MAINcreate_cost_constraint();
    
    strStatus='PROGRESSION: computing the detailed scenario analysis results...';
    disp(strStatus);
    fprintf(fidStatus, '%s\n',strStatus);
    
    %     sol = buildSolution(x_sol, fixed_x_scenariomode, [], ...
    %         originalcostconstr'*(x_sol), ...
    %         originalcostconstr'*(x_sol));
    % MM Version (one extra parameter)
    sol = buildSolution_(x_sol, fixed_x_scenariomode, [], ...
        originalcostconstr'*(x_sol), ...
        originalcostconstr'*(x_sol), aggregationInfo, commonDataInfo);
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint29a_new', 'sol');
    
    solutionList = struct('X', {}, 'X_free', {}, 'AQIs', {}, ...
        'COST', {}, 'BUDGET', {}, ...
        'COSTPERMACROSECTOR', {}, ...
        'BUDGETPERMACROSECTOR', {});
    
    solutionList = [solutionList, sol];
    
    %20130403 - fixed bug on cost of detailed scenario mode
    solutionSet = [solutionSet, ...
        struct('BUDGET', sol.COST, 'SOL_LIST', solutionList)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %AGGREGATED SCENARIO MODE
elseif commonDataInfo.optim_flags.mode_ce_mo==3   %aggregated scenario analysis
    %% MM: no quadrants needed
    
    strStatus='PROGRESSION: computing the aggregated scenario analysis results...';
    disp(strStatus);
    fprintf(fidStatus, '%s\n',strStatus);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LOAD EMISSION REDUCTION PERCENTAGES
    %     aggsa=importdata(strcat(pathEMI,'aggregated_ads_perc_reductions.txt'));
    aggsa=importdata(commonDataInfo.dirs.pathADS);
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\sa1_fix', 'aggsa');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %COMPUTE REDUCED EMISSIONS (BOTH AREAL AND POINT...WE EXPECT BOTH TYPE
    %OF FILES)
    %conisder the seasons - yea, win, sum
    for sea=1:3
        for i=1:commonDataInfo.domainInfo.ncel
            %emiCleRed contains both CLE {sea,1} and REDUCED {sea,2}
            %emission fields
            %emi init
            emiCleRed{sea,1}(i,:)=sum(emi{i,1}(:,1+12*(sea-1):12+12*(sea-1))); %CLE contains emissions
            emiCleRed{sea,2}(i,:)=sum(emi{i,1}(:,1+12*(sea-1):12+12*(sea-1))); %emi red contains initially the same CLE emissions
            
            %update PAD emissions, considering the emission reduction
            %application...NB: only consider cells completely inside
            %regione
            if flag_optim_dom(i)==1
                emiCLEtmp{i,1}=emi{i,1}(:,1+12*(sea-1):12+12*(sea-1))-aggsa.data/100.*emi{i,1}(:,1+12*(sea-1):12+12*(sea-1));
                %emi red
                emiCleRed{sea,2}(i,:)=sum(emiCLEtmp{i,1});
            end
        end
    end
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\sa2_fix', 'emiCleRed');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %SOME LABELS
    nomeaqi={'PM10','PM25','AOT40','SOMO35','MAX8H','NO2','PM10dailyExceed'};
    aqiUM={'[microg/m3]','[microg/m3]','[microg/m3*h]','[microg/m3*d]','[microg/m3]','[microg/m3]','[# days]'};
    nomeemi={'NOX','VOC','NH3','PM10','PM25','SO2'};
    emiUM={'[ton/year]','[ton/year]','[ton/year]','[ton/year]','[ton/year]','[ton/year]'};
    nomesea={'TP1','TP2','TP3'};
    delim=',';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %CREATE DIRS
    mkdir(strcat(commonDataInfo.dirs.pathOUF,'\maps_aqi'));
    mkdir(strcat(commonDataInfo.dirs.pathOUF,'\maps_emi'));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %SAVING EMISSIONS, OF CLE AND OF REDUCED EMISSIONS
    for indpar=1:2 %2 points to be produced
        %UNIBS(ET)20130927
        %loop over seasons introduced only if data for that season exists
        for sea2=1:3
            if commonDataInfo.special.flag_ADS_tp(sea2)==2
                emiADS=emiCleRed{sea,indpar};
                for emiind=1:6
                    emiADSFinal=emiADS(:,emiind)+emiADS(:,emiind+6); %sum areal and point
                    emiADSFinal(flag_optim_dom==0)=-999;
                    file=strcat(commonDataInfo.dirs.pathOUF,'/maps_emi/',nomeemi{emiind},nomesea{sea2},'_point',int2str(indpar),'.csv');
                    fid=fopen(file,'wt');
                    fprintf(fid,'%s %c %s %c %s \n','xutm[km]', delim, 'yutm[km]', delim, strcat(nomeemi{emiind},emiUM{emiind}));
                    dlmwrite(file,[commonDataInfo.special.xutm commonDataInfo.special.yutm emiADSFinal],'-append','roffset', 0,'precision','%6.2f');
                    fclose(fid);
                end
            end
        end
    end
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\sa3_fix', 'emiADSFinal');
    
    %LOOP OVER AQIS - CONSIDER ONLY YEARLY EMISSIONS AND AQIS
    for ii = 1:3, %season
        for jj = 1:commonDataInfo.info.AQINum, %aqi
            %UNIBS(ET)20130927
            %Compute AQI only if data for that season exists
            if commonDataInfo.special.flag_ADS_tp(ii)==2
                if (isequal(strtrim(commonDataInfo.pathANN(ii).ANNs(jj,:)),'-999')==0)
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %COMPUTE AQI PER CELL, RELATED TO INITIAL AND REDUCED CASES
                    %CONSTRAINTS ON ANNS IDENTIFICATION BOUNDS ARE CHECKED HERE
                    NN=neuralNet_get_nnSuperSet_indexed(aggregationInfo.mathIntermediateData, ii, jj);
                    [aqiCLEini]=MAINaggregated_scenario_mode__(emiCleRed{ii,1}, NN,ii,jj, commonDataInfo, aggregationInfo);
                    [aqiCLEred]=MAINaggregated_scenario_mode__(emiCleRed{ii,2}, NN,ii,jj, commonDataInfo, aggregationInfo);
                    
                    %if AQI=#exceed > threshold
                    if(jj == 7)
                        % EU level relation, fron Philippe Thunis, from Chimere/Urban stations
                        aqiCLEini = pm10aveToExceed(1).*aqiCLEini-pm10aveToExceed(2);
                        aqiCLEini(aqiCLEini<0)=0;
                        
                        aqiCLEred = pm10aveToExceed(1).*aqiCLEred-pm10aveToExceed(2);
                        aqiCLEred(aqiCLEred<0)=0;
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %put infos in one single structure
                    aqiTOT{1}=aqiCLEini;
                    aqiTOT{2}=aqiCLEred;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %SAVING REDUCED AQI
                    for indpar=1:2 %2 points to be produced
                        aqiADS=aqiTOT{indpar};
                        
                        aqiFinal=zeros(length(flag_optim_dom),1);
                        aqiFinal(flag_optim_dom==0)=-999;
                        aqiFinal(find(flag_optim_dom==1 | flag_optim_dom==2))=aqiADS;
                        
                        file=strcat(commonDataInfo.dirs.pathOUF,'/maps_aqi/',nomeaqi{jj},nomesea{ii},'_point',int2str(indpar),'.csv');
                        fid=fopen(file,'wt');
                        fprintf(fid,'%s %c %s %c %s \n','xutm[km]', delim, 'yutm[km]', delim, strcat(nomeaqi{jj},aqiUM{jj}));
                        dlmwrite(file,[commonDataInfo.special.xutm commonDataInfo.special.yutm aqiFinal],'-append','roffset', 0,'precision','%6.2f');
                        fclose(fid);
                    end
                end
            end
        end
    end
    
    Error=0;
    fprintf(fidExit, int2str(Error));
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       DATA STRUCTURES RESTORING
%                  (DATA STRUCTURES REQUIRED IN POSTmain)

paretoSolutions = struct('X', {}, 'X_free', {}, 'AQIs', {}, ...
    'COST', {}, 'BUDGET', {}, ...
    'COSTPERMACROSECTOR', {}, ...
    'BUDGETPERMACROSECTOR', {});

[foo, num_of_budgets] = size(solutionSet);
for b=1:num_of_budgets,
    [foo, num_of_points] = size(solutionSet(b).SOL_LIST);
    for p=1:num_of_points,
        paretoSolutions = [paretoSolutions, solutionSet(b).SOL_LIST(p)];
    end;
end;

if (commonDataInfo.optim_flags.mode_ce_mo==0 | commonDataInfo.optim_flags.mode_ce_mo==1 | commonDataInfo.optim_flags.mode_ce_mo==2)
    CLE = global_data(:,19)/100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save results in mat file
% 20160421 MM FG Vs quadrant/NN
% Quadrant/NN
%localnnSuperSet=interface_get_nnSuperSet(aggregationInfo.mathIntermediateData);
%localDSuperSet=interface_get_DSuperSet(aggregationInfo.geometryIntermediateData);
%localbcSuperSet=interface_get_bcSuperSet(aggregationInfo.mathIntermediateData);
% FirstGuess
localnnSuperSet=0;
localDSuperSet=FG_get_aqi_D(aggregationInfo.geometryIntermediateData, 1);
localbcSuperSet=0;
%end
save(nomeOUTPUT,'paretoSolutions','localDSuperSet','localnnSuperSet','solutionSet');
toc

strStatus='PROGRESSION: Computation part is finished.';
disp(strStatus);
fprintf(fidStatus, '%s\n',strStatus);

% fprintf(fidStatus, '%s\n','Computation part is finished.');

if commonDataInfo.optim_flags.mode_ce_mo<3
    % %postprocessing of results
    % 20160421 MM FG Vs quadrant/NN
    % Quadrant/NN
    %localnnSuperSet=interface_get_nnSuperSet(aggregationInfo.mathIntermediateData);
    %localDSuperSet=interface_get_DSuperSet(aggregationInfo.geometryIntermediateData);
    %localbcSuperSet=interface_get_bcSuperSet(aggregationInfo.mathIntermediateData);
    % FirstGuess
    localnnSuperSet=0;
    localDSuperSet=FG_get_DSuperSet(aggregationInfo.geometryIntermediateData, 1);
    localbcSuperSet=0;
    %end
    % MM Version load locally needed structures...
    %POSTmain(f2,nnSuperSet,DSuperSet,pathANN,domainInfo.flag_optim_dom, domainInfo.flag_aqi_dom,...
    % 20160422 : MM Version Quadrant/NN Version
    %POSTmain(f2,localnnSuperSet,localDSuperSet,localbcSuperSet,commonDataInfo.pathANN, commonDataInfo.flag_region_dom, flag_optim_dom, flag_aqi_dom,...
    %    commonDataInfo.domainData.coordinate,emi,global_data,base_emi_low,base_emi_high,CLE, ...
    %    base_emi_low_noc,base_emi_high_noc,commonDataInfo.dirs.pathOCM,paretoSolutions,...
    %    commonDataInfo.optim_flags.areal_point,commonDataInfo.optim_flags.reg_net,commonDataInfo.dirs.pathOUF,pop,commonDataInfo.info.nocID,commonDataInfo.dirs.pathAR,pm10aveToExceed,...
    %    commonDataInfo.cell_threshold_set,fidStatus,fidExit)
    % 20160422 : MM Version First Guess/SR Version
    global_data=commonDataInfo.special.global_data;
    POSTmain_FG(f2,localDSuperSet,aggregationInfo,commonDataInfo.pathANN, commonDataInfo.flag_region_dom, flag_optim_dom, flag_aqi_dom,...
        commonDataInfo.domainData.coordinate,emi,global_data,base_emi_low,base_emi_high,CLE, ...
        base_emi_low_noc,base_emi_high_noc,commonDataInfo.dirs.pathOCM,paretoSolutions,...
        commonDataInfo.optim_flags.areal_point,commonDataInfo.optim_flags.reg_net,commonDataInfo.dirs.pathOUF,pop,commonDataInfo.info.nocID,commonDataInfo.dirs.pathAR,pm10aveToExceed,...
        commonDataInfo.cell_threshold_set,fidStatus,fidExit, commonDataInfo)
    %     POSTmain(f2,localnnSuperSet,localDSuperSet,pathANN,domainInfo.flag_optim_dom, domainInfo.flag_aqi_dom,...
    %         domainData.coordinate,emi,global_data,base_emi_low,base_emi_high,CLE, ...
    %         base_emi_low_noc,base_emi_high_noc,dirs.pathOCM,paretoSolutions,...
    %         optim_flags.areal_point,1,dirs.pathOUF,domainInfo.pop,info.nocID,dirs.pathAR,pm10aveToExceed,...
    %         cell_threshold_set,fidStatus,fidExit);
    %optim_flags.areal_point,flag_optim_flags.reg_net,dirs.pathOUF,domainInfo.pop,info.nocID,dirs.pathAR,pm10aveToExceed,...
end

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%                          NESTED FUNCTIONS                              %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function costconstr = MAINcreate_cost_constraint()
        %create cost constraints
        
        %sn1g	secg	actg	tecg	hl_flag		reNOx		reCOV		reNH3
        %rePM10 rePM2.5 	reSO2  emfNOx  emfCOV  emfNH3  emfPM10  emfPM2.5  emfSO2
        %unitcost   arCLE_year   arPOT_year  replace  optim
        %removal efficiencies
        uc = global_data(:,18);
        
        %emission order: NOX, COV, NH3, PM10, PM25, SO2
        
        %remove the last rows of the actlev_final_sum file (the noc-related rows)
        % actlev_final_sum(length(actlev_final_sum)-12:length(actlev_final_sum),:)=[];
        actlev_final_sum(size(global_data,1)+1:end,:) = [];
        costconstr = uc.*actlev_final_sum;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function minimumCostPerMacrosector = computeMinimumCostPerMacrosector()
        
        costconstr=MAINcreate_cost_constraint();
        
        % for each macrosector build the mask (costPerMacrosectorMask) to
        % identify the coefficients associated with the technologies which
        % apply to the given macrosector
        
        [num_of_technologies, foo] = size(flag_usab_tech);
        internalCostPerMacrosectorMask = zeros(num_of_technologies,11);
        for mm=1:11,
            for tt=1:num_of_technologies,
                if(global_data(tt,1) == mm)
                    internalCostPerMacrosectorMask(tt,mm) = 1;
                end
            end;
        end
        
        % for each macrosector derive from the  costconstr's coefficients
        % the coefficients associated with the technologies which apply to
        % the given macrosector (internalCostPerMacrosector)
        
        internalCostPerMacrosector = zeros(num_of_technologies,11);
        for mm=1:11,
            internalCostPerMacrosector(:,mm) = ...
                internalCostPerMacrosectorMask(:,mm) .* costconstr;
        end;
        
        % compute the cost for each macrosector when technologies are
        % applied at the CLE application rate
        minimumCostPerMacrosector = zeros(11,1);
        for mm=1:11,
            minimumCostPerMacrosector(mm,1) = ...
                internalCostPerMacrosector(:,mm)'*CLE;
        end;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function defaultTrafficConstraintsSubmatrix = ...
            generateTrafficConstraints(TRAFFIC_MATRIX_1, TRAFFIC_MATRIX_2, num_of_technologies_tc)
        
        defaultTrafficConstraintsSubmatrix = [];
        
        [TM1_rows_num, foo] = size(TRAFFIC_MATRIX_1);
        [TM2_rows_num, foo] = size(TRAFFIC_MATRIX_2);
        
        for i_TM1=1:TM1_rows_num,
            TM1_index = TRAFFIC_MATRIX_1(i_TM1,1);
            TM1_activity = TRAFFIC_MATRIX_1(i_TM1,2);
            TM1_technology = TRAFFIC_MATRIX_1(i_TM1,3);
            TM1_technology_H_L = TRAFFIC_MATRIX_1(i_TM1,4);
            
            for i_TM2=1:TM2_rows_num,
                TM2_index = TRAFFIC_MATRIX_2(i_TM2,1);
                TM2_activity = TRAFFIC_MATRIX_2(i_TM2,2);
                TM2_technology = TRAFFIC_MATRIX_2(i_TM2,3);
                TM2_technology_H_L = TRAFFIC_MATRIX_2(i_TM2,4);
                
                if((TM1_activity == TM2_activity) && ...
                        (TM1_technology == TM2_technology) && ...
                        (TM1_technology_H_L == TM2_technology_H_L))
                    
                    min_index = -1;
                    max_index = -1;
                    if(TM1_index < TM2_index)
                        min_index = TM1_index; max_index = TM2_index;
                    else
                        max_index = TM1_index; min_index = TM2_index;
                    end;
                    
                    defaultTrafficConstraint = [];
                    
                    for ii=1:min_index-1
                        defaultTrafficConstraint = ...
                            [defaultTrafficConstraint, 0];
                    end;
                    defaultTrafficConstraint = ...
                        [defaultTrafficConstraint, 1];
                    for ii=min_index+1:max_index-1
                        defaultTrafficConstraint = ...
                            [defaultTrafficConstraint, 0];
                    end;
                    defaultTrafficConstraint = ...
                        [defaultTrafficConstraint, -1];
                    for ii=max_index+1:num_of_technologies_tc
                        defaultTrafficConstraint = ...
                            [defaultTrafficConstraint, 0];
                    end;
                    
                    defaultTrafficConstraintsSubmatrix = ...
                        [defaultTrafficConstraintsSubmatrix; defaultTrafficConstraint];
                    
                    break;
                    
                end;
                
            end;
            
        end;
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function trafficConstraintsMatrix = ...
            generateTrafficConstraintsMatrix(hig, urb, ext)
        
        % find in the DB of technologies the entries associated with
        % highway, urban and extra-urban activities
        global_data_tc=load(dirs.pathOCM);
        [num_of_technologies_tc, foo] = size(global_data_tc);
        
        % distinguish between technical and non-technical measures:
        % technical measures: flag_tech_nontech = 1
        % non-technical measures: flag_tech_nontech = 0
        flag_tech_nontech_tc = global_data_tc(:,23);
        
        hig_matrix = [];
        urb_matrix = [];
        ext_matrix = [];
        for tt=1:num_of_technologies_tc,
            if((global_data_tc(tt,2) == hig) && ...
                    (flag_tech_nontech_tc(tt) == 1))
                hig_matrix = [hig_matrix;
                    tt, global_data_tc(tt,3), global_data_tc(tt,4), global_data_tc(tt,5)];
            end
            if((global_data_tc(tt,2) == urb) && ...
                    (flag_tech_nontech_tc(tt) == 1))
                urb_matrix = [urb_matrix;
                    tt, global_data_tc(tt,3), global_data_tc(tt,4), global_data_tc(tt,5)];
            end
            if((global_data_tc(tt,2) == ext) && ...
                    (flag_tech_nontech_tc(tt) == 1))
                ext_matrix = [ext_matrix;
                    tt, global_data_tc(tt,3), global_data_tc(tt,4), global_data_tc(tt,5)];
            end
        end;
        
        % define the constraints matching the technologies associated with
        % highway and urban activities
        
        defaultTrafficConstraintsSubmatrixA = ...
            generateTrafficConstraints(hig_matrix, urb_matrix, ...
            num_of_technologies_tc);
        
        % remove from the constraints the coefficients associated with
        % nocID technologies
        
        defaultTrafficConstraintsSubmatrixA = defaultTrafficConstraintsSubmatrixA';
        
        trafficConstraintsSubmatrixA = ...
            defaultTrafficConstraintsSubmatrixA(find(~(global_data_tc(:,4)==info.nocID)),:);
        
        trafficConstraintsSubmatrixA = trafficConstraintsSubmatrixA';
        
        % define the constraints matching the technologies associated with
        % urban and extra-urban activities
        
        defaultTrafficConstraintsSubmatrixB = ...
            generateTrafficConstraints(urb_matrix, ext_matrix,...
            num_of_technologies_tc);
        
        % remove from the constraints the coefficients associated with
        % nocID technologies
        
        defaultTrafficConstraintsSubmatrixB = defaultTrafficConstraintsSubmatrixB';
        
        trafficConstraintsSubmatrixB = ...
            defaultTrafficConstraintsSubmatrixB(find(~(global_data_tc(:,4)==info.nocID)),:);
        
        trafficConstraintsSubmatrixB = trafficConstraintsSubmatrixB';
        
        
        % define the final matrix
        trafficConstraintsMatrix = ...
            [trafficConstraintsSubmatrixA;...
            trafficConstraintsSubmatrixB];
        
        
        clear defaultTrafficConstraintsMatrix;
        clear global_data_tc;
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MM Version (one extra parameter)
    function newAgg = Update_(aggregationInfo, commonDataInfo)
        
        %if techs change emissions less than threshold (sum for all
        %pollutants) the techs are disregarded
        DOptSet=neuralNet_get_FullDOptSet(aggregationInfo.geometryIntermediateData);
        %DOptSet=interfacet_get_FullDOptSet(aggregationInfo.geometryIntermediateData);
        %if techs change emissions less than threshold (sum for all
        %pollutants) the techs are disregarded
        threshold = 1e-6;
        
        final_active_technologies = zeros(num_of_technologies,1);
        
        for ii = 1:commonDataInfo.optim_flags.optAQINum,
            
            %in case of areal+point summed, D and d contain infos both or areal
            %and point, together
            %in case of areal and point separated, D and d are only related to
            %areal emissions (the point infos Dp and dp will be used later on)
            D = DOptSet(ii).D;
            d = DOptSet(ii).d;
            
            totalReductionFactor = zeros(num_of_technologies,1);
            
            for tt = 1:num_of_technologies,
                %check for each tech the % of reduced emission (D*UB) over
                %initial emissions (d)
                ReductionFactors = (D(:,tt) .* (UB(tt))) ./ d;
                ReductionFactors = abs(ReductionFactors);
                ReductionFactors = ReductionFactors .* 100.0;
                %ENR 20130318 - error when d=0...this create NaN elements,
                %and sum(ReductionFactors)=NaN
                ReductionFactors(isnan(ReductionFactors))=0;
                %sum over pollutants-quadrants-cells...is the tech reducing
                %emissions?
                totalReductionFactor(tt,1) = sum(ReductionFactors);
                
            end;
            
            %check also for point emissions
            if (commonDataInfo.optim_flags.areal_point==1)
                D = DOptSet(ii).Dp;
                d = DOptSet(ii).dp;
                
                for tt = 1:num_of_technologies,
                    %the same as for areal techs
                    ReductionFactors = (D(:,tt) .* (UB(tt))) ./ d;
                    ReductionFactors = abs(ReductionFactors);
                    ReductionFactors = ReductionFactors .* 100.0;
                    %ENR 20130318 - error when d=0...this create NaN elements,
                    %and sum(ReductionFactors)=NaN
                    ReductionFactors(isnan(ReductionFactors))=0;
                    %cumulate the effect of both areal and point emissions
                    totalReductionFactor(tt,1) = ...
                        totalReductionFactor(tt,1) + sum(ReductionFactors);
                    
                end;
                
            end
            
            %these are techs that will be used in optimization
            %in the beginning all are used, except ...
            active_technologies = ones(num_of_technologies,1);
            
            for tt = 1:num_of_technologies,
                %table of the cases to be managed in the frame of smartupdate.
                %In these cases, put active_technologies(tt,1) = 0;
                
                %no substitution configuration:
                %    OPTIMIZABLE     REPLACEABLE    EFFICIENT   UB=LB
                %         0
                %         1                             0
                %         1                             0         X
                %substitution configuration
                %    OPTIMIZABLE     REPLACEABLE    EFFICIENT   UB=LB
                %         0               1
                %         1               1             0         X
                % (OPTIMIZABLE=1 means optimizable
                %  REPLACEABLE=1 means not replaceable
                %  EFFICIENT=1 means efficient)
                %  UB=LB
                
                if (((commonDataInfo.optim_flags.constraints == 0) && (flag_usab_tech(tt,1) == 0)) || ...
                        ((commonDataInfo.optim_flags.constraints == 0) && (flag_usab_tech(tt,1) == 1) && (totalReductionFactor(tt,1) <= threshold)) || ...
                        ((commonDataInfo.optim_flags.constraints == 0) && (flag_usab_tech(tt,1) == 1) && (abs(UB(tt)-LB(tt))<=1e-4)) || ...
                        ((commonDataInfo.optim_flags.constraints == 1) && (flag_usab_tech(tt,1) == 0) && (flag_repl_tech(tt,1) == 1)) || ...
                        ((commonDataInfo.optim_flags.constraints == 1) && (flag_usab_tech(tt,1) == 1) && (flag_repl_tech(tt,1) == 1) && (totalReductionFactor(tt,1) <= threshold)) || ...
                        ((commonDataInfo.optim_flags.constraints == 1) && (flag_usab_tech(tt,1) == 1) && (flag_repl_tech(tt,1) == 1) && (abs(UB(tt)-LB(tt))<=1e-4)))
                    active_technologies(tt,1) = 0;
                end
                
            end;
            
            final_active_technologies = final_active_technologies | active_technologies;
            
        end;
        
        %subset of control variable to fix (not used in optimization)
        inhibited_x = zeros(num_of_technologies,1);
        
        for tt = 1:num_of_technologies,
            if(final_active_technologies(tt,1)==0)
                %this tech is not used in optimization....keep track of its
                %CLE value
                inhibited_x(tt,1) = CLE(tt,1);
                fixed_x(tt,1) = 1;
            end;
        end;
        
        for ii = 1:commonDataInfo.optim_flags.optAQINum,
            %in case of areal+point summed, D and d contain infos both or areal
            %and point, together
            %in case of areal and point separated, D and d are only related to
            %areal emissions (the point infos Dp and dp will be used later on)
            DOptSet(ii).d = DOptSet(ii).d - (DOptSet(ii).D * inhibited_x);
            if (commonDataInfo.optim_flags.areal_point==1)
                DOptSet(ii).dp = DOptSet(ii).dp - (DOptSet(ii).Dp * inhibited_x);
            end
            
        end
        
        %remove from constraints the part of the problem related to the
        %techs not to be used in optimization
        B = B - (A * inhibited_x);
        
        
        %part of cost related to techs at CLE (not to be optimizized)
        %also techs not to be optimized and with with cost <0 are kept at
        %CLE, to allow for perfect computation of CLE
        deltaCost = sum(costconstr.*(inhibited_x));
        
        % update of UB - keep only technologies to be optimized
        UB = UB(final_active_technologies==1);
        
        % update of LB - keep only technologies to be optimized
        LB = LB(final_active_technologies==1);
        
        % update of CLE - keep only technologies to be optimized
        CLE = CLE(final_active_technologies==1);
        
        % update of Dx - keep only technologies to be optimized
        for ii = 1:commonDataInfo.optim_flags.optAQINum,
            DOptSet(ii).D = DOptSet(ii).D(:,(final_active_technologies==1));
            if (commonDataInfo.optim_flags.areal_point==1)
                DOptSet(ii).Dp = DOptSet(ii).Dp(:,(final_active_technologies==1));
            end
        end
        
        % keep only technologies to be optimized
        % update of A;
        A = A(:,(final_active_technologies==1));
        % update of A1;
        A1 = A1(:,(final_active_technologies==1));
        % update of costconstr;
        costconstr = costconstr(final_active_technologies==1);
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\DOptSet_before_new', 'DOptSet');
        newAgg=neuralNet_set_FullDOptSet(aggregationInfo, DOptSet);
        DOptSetWrong=neuralNet_get_FullDOptSet(aggregationInfo.geometryIntermediateData);
        DOptSetRight=neuralNet_get_FullDOptSet(newAgg.geometryIntermediateData);
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\DOptSet_after_new', 'DOptSetWrong', 'DOptSetRight');
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MM Version (one extra parameter)
    function sol = buildSolution_(x_free, fixed_x, inhibited_x, ...
            cost, budget, aggInfo, commonDataInfo)
        
        % solution restoring - go back to full set of control variables
        % (both optimized and kept fixed)
        x = zeros(num_of_technologies,1);
        index = 1;
        for tt = 1:num_of_technologies,
            % if tech has remained fixed
            if(fixed_x(tt,1)==1)
                %retrieve the CLE
                x(tt,1)=inhibited_x(tt,1);
            else
                % else if tech has been used in optimization, retrieve
                % optimal solution
                x(tt,1)=x_free(index,1);
                index = index + 1;
            end;
        end;
        
        cost_Per_Macrosector = zeros(commonDataInfo.MNum,1);
        for mm=1:commonDataInfo.MNum,
            cost_Per_Macrosector(mm,1) = ...
                costPerMacrosector(:,mm)'*x;
        end;
        
        %computing the AQIS for the 3 temporal periods, 6 aqis, 3
        %aggregation types
        aqi = ones(3,commonDataInfo.info.AQINum,3) * -999;
        %isDelta=interface_get_isDelta(aggInfo);
        %isAggregated=interface_get_isAggregated(aggInfo);
        
        for ii = 1:3,
            for jj = 1:commonDataInfo.info.AQINum,
                if (isequal(strtrim(commonDataInfo.pathANN(ii).ANNs(jj,:)),'-999')==0)
                    % MM Version (added D get from function and some parameters in the call)
                 fName=strtrim(commonDataInfo.pathANN(ii).ANNs(jj,:));
                 aggInfo.firstguess=0;
                 if (strcmp(fName,'-999') == 0)
                    [alpha, omega, radius, flatWeight, pollutantList]=FG_read(fName);
                    aggInfo.firstguess.alpha=alpha;
                    aggInfo.firstguess.omega=omega;
                    aggInfo.firstguess.radius=radius;
                    aggInfo.firstguess.flatWeight=flatWeight;
                    aggInfo.firstguess.pollutantList=pollutantList;
                 end
                    
                    %cycleD=interface_get_DSuperSetFunction_indexed(aggInfo.mathIntermediateData, ii, jj);
                    cycleD=quadrant_get_DSuperSet_indexed(aggInfo.geometryIntermediateData, ii, jj);
                    %cyclebc=quadrant_get_bcSuperSet_indexed(aggInfo.mathIntermediateData, ii, jj);
                    % 20160420 MM temporary setting
                    aggInfo.mathIntermediateData=0;
                    %cycleNN=neuralNet_get_nnSuperSet_indexed(aggInfo.mathIntermediateData, ii, jj);
                    % 20160420 MM temporary setting
                    cycleNN=0;
                    cyclebc=0;
                    % FG version
                    aqi(ii,jj,1) = MAINcompute_aqi_(x, cycleNN, cycleD, 0,commonDataInfo.cell_threshold_set(jj),jj-1,1, cyclebc, aggInfo, commonDataInfo, ii, jj);
                    % FG version only year period (third index=1)...
                    aqi(ii,jj,2) = MAINcompute_aqi_(x, cycleNN, cycleD, 1,commonDataInfo.cell_threshold_set(jj),jj-1,1, cyclebc, aggInfo, commonDataInfo, ii, jj);
                    % FG version only year period (third index=1)...
                    aqi(ii,jj,3) = MAINcompute_aqi_(x, cycleNN, cycleD, 2,commonDataInfo.cell_threshold_set(jj),jj-1,1, cyclebc, aggInfo, commonDataInfo, ii, jj);
                    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint12_new', 'aqi');
                end
                
            end
        end
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint13_new', 'aqi');
        
        sol = struct('X', x, 'X_free', x_free, 'AQIs', aqi, ...
            'COST', cost, 'BUDGET', budget, ...
            'COSTPERMACROSECTOR', cost_Per_Macrosector, ...
            'BUDGETPERMACROSECTOR', budget * commonDataInfo.MPBudget);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MM Version (added parameters)
%x,D,BCset, aggregationInfo, commonDataInfo, periodIndex, aqiIndex
    %function aqi_per_cell=MAINcompute_aqi_per_cell_(x, DD, bcSet, aggregationInfo,commonDataInfo,periodIndex,aqiIndex)
    function aqi_per_cell=MAINcompute_aqi_per_cell_(x, NN, DD, bcSet,  aggregationInfo,commonDataInfo,periodIndex,aqiIndex)
        
        % in case of areal+point summed, D and d contain infos both or areal
        % and point, together
        % in case of areal and point separated, D and d are only related to
        % areal emissions (the point infos Dp and dp will be used later on)
        D=DD.D;
        d=DD.d;
        
        restoreDataInfo.areal_point=commonDataInfo.optim_flags.areal_point;
        
        % E is the reduced emissions, depending on application rates
        % e are the initial emissions
        
        % calculation starts from sparse matrices
        % delta or absolute (x = 0 or CLE)
        %TODO implement delta or absolute
        E = D * sparse(x);
        E = d - E;
        E_full = full(E);
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint5_new','E_full');
        
        %for quadrant case
        %from global matrix, extract quadrant informations
        %MM check split in 4 quadrant case
        
        %!!!!Load radius 4 or 16
        % 20160420 fg version
        emis=FG_rebuildOrderedEmis(E_full, aggregationInfo.geometryIntermediateData, commonDataInfo);
        % 20160420 nn/quadrant version
        %emis=interface_rebuildOrderedEmis(E_full, aggregationInfo.geometryIntermediateData, commonDataInfo);
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint6_new','emis');
        
        %CASE OF AREAL+POINT, SUMMED
        if commonDataInfo.optim_flags.areal_point==0
            
            emisP=0;
            %model indepent NO
            
            % 20160420 nn/quadrant version
            bcSet=0;
            NN=0;
            aggregationInfo.mathInfo=0;
            % 20160420 nn/quadrant version
            %emissioni=interface_buildEmission(emis, emisP, bcSet, NN, aggregationInfo.mathInfo,'A',periodIndex,aqiIndex);
            % 20160420 fg version
            emissioni=FG_buildEmission(emis, emisP, NN, bcSet, aggregationInfo.mathInfo, aggregationInfo,'A',periodIndex,aqiIndex);
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPointEmissioni_a_new','emissioni');
            
            
            %CASE FOR POINT EMISSIONS (AREAL INFOS HAVE BEEN PROCESSED IN
            %THE LINES ABOVE)
        elseif commonDataInfo.optim_flags.areal_point==1
            Dp=DD.Dp;
            dp=DD.dp;
            
            % check here absolute or Delta emission...
            Ep = Dp * sparse(x);
            Ep = dp - Ep;
            Ep_full = full(Ep);
            
            % 20160420 fg version
            emisP=FG_rebuildOrderedEmis(Ep_full, aggregationInfo.geometryIntermediateData, commonDataInfo);
            % 20160420 nn/quadrant version
            %emisP=interface_rebuildOrderedEmis(Ep_full, aggregationInfo.geometryIntermediateData, commonDataInfo);
            emissioni=interface_buildEmission(emis, emisP, NN, bcSet, aggregationInfo.mathDataInfo,aggregationInfo.mathIntermediateData,'P',periodIndex,aqiIndex);
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPointEmissioni_p_new','emissioni');
        end
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint7_new','emissioni');
            % 20160420 fg version
         pathANN=commonDataInfo.pathANN;
         fName=strtrim(pathANN(periodIndex).ANNs(aqiIndex,:));
         %strtrim(pathANN(k).ANNs(indaqi,:))
         % nn case
         % 20160418: MM & EP net/quadrants input
         %[net]=net_read(fName);
         % commonDataInfo.radius=net.icells;
         % 20160418: MM & EP regression input
         
         aqi_per_cell=FG_do_aqi_per_cell(emissioni, NN, aggregationInfo, commonDataInfo, periodIndex, aqiIndex);
            % 20160420 nn/quadrant version
        %aqi_per_cell=do_interface_get_aqi_per_cell(emissioni, NN, aggregationInfo, commonDataInfo, periodIndex, aqiIndex);
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function aqi_val = ...
%        MAINcompute_aqi_(x,nn,D,obj_function_type,cell_threshold,AQIId,REAL,prepareDataInfo,finalizeJobInfo,restoreDataInfo)
%function aqi_val = ...
%        MAINcompute_aqi_(x,D,obj_function_type,cell_threshold,AQIId,REAL,prepareDataInfo,finalizeJobInfo,restoreDataInfo)
%obj_function_type function applied
%obj_function_type function applied
%one single object split by section (inner structures, for example)
    function aqi_val = ...
            MAINcompute_aqi_(x, NN, D, obj_function_type, cell_threshold, AQIId, REAL, BCset, aggregationInfo, commonDataInfo, periodIndex, aqiIndex)
        
        % compute aqi per cell
        %30140403 ET - Added BCset input
        %aqi_per_cell = MAINcompute_aqi_per_cell_(x,nn,D,BCset);
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint4_new',...
        %'x', 'D', 'obj_function_type', 'cell_threshold', 'AQIId', 'REAL','BCset');
        
        aqi_per_cell = MAINcompute_aqi_per_cell_(x, NN,D,BCset, aggregationInfo, commonDataInfo, periodIndex, aqiIndex);
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint10_new', 'aqi_per_cell');
        
        % managing the fact that PAD and ACD are different
        % define indexes for pad domain
        indpad=(flag_optim_dom==1 | flag_optim_dom==2);
        % define updated indaqi
        indacd=flag_aqi_dom(indpad);
        % remove from aqi_per_cell what is not in the acd
        aqi_per_cell(indacd==0)=[];
        
        if(AQIId == 6)
            % EU level relation, fron Philippe Thunis, from Chimere/Urban stations
            aqi_per_cell = ...
                pm10aveToExceed(1).*aqi_per_cell-pm10aveToExceed(2);
            %             thresDailyPm10=50;
            %             aqi_per_cell_T = ones(size(aqi_per_cell_mod)) * thresDailyPm10;
            %             aqi_per_cell_T = aqi_per_cell_mod - aqi_per_cell_T;
            aqi_per_cell (aqi_per_cell < 0) = 0;
            
            %             aqi_per_cell = aqi_per_cell_T;
        end
        
        switch obj_function_type
            case 0
                % average
                aqi_val = mean(aqi_per_cell);
                
            case 1
                % number of cells over thresold
                T = ones(size(aqi_per_cell)) * cell_threshold;
                
                if(REAL == 1)
                    aqi_val = sum(aqi_per_cell > T);
                else
                    T = aqi_per_cell - T;
                    T(T < 0) = 0;
                    % Theta'
                    aqi_val = sum(T);
                    
                    % Theta''
                    %T = T.^(1/4);
                    %aqi_val = sum(T);
                end
                
            case 2
                % population average
                indpop=commonDataInfo.domainInfo.pop;
                indpop(commonDataInfo.domainInfo.flag_aqi_dom==0)=[];
                aqi_val=sum(aqi_per_cell'.*indpop)/sum(indpop);
        end
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint11_new', 'aqi_val');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     function aqi_val=MAINcompute_aqi_mp(x,optAQIBestValues)
%
%         optAQIValues = zeros(optim_flags.optAQINum,1);
%
%         for ii = 1:optim_flags.optAQINum,
%
%             optAQIValues(ii,1) = MAINcompute_aqi(x,nnOptSet(ii),DOptSet(ii),...
%                 aqi_obj_function_type(ii),...
%                 cell_threshold_set(aqi_obj(ii)+1),...
%                 aqi_obj(ii),...
%                 0);
%
%             optAQIcle(ii,1)=sol_CLE.AQIs(aqi_horizon(ii),aqi_obj(ii)+1,aqi_obj_function_type(ii)+1);
%
%         end
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %NEW VERSION
%         %20130430 - testing the claudio proposal to improve objective
%         %function computation
%         delta_init=optAQIcle-optAQIBestValues;
%         delta_prod=ones(optim_flags.optAQINum,1);
%         for ii=1:optim_flags.optAQINum
%             for kk=1:optim_flags.optAQINum
%                 if kk~=ii
%                     delta_prod(ii,1)=delta_prod(ii,1).*delta_init(kk,1);
%                 end
%             end
%         end
%         if aqi_weights_init <= 1
%             % user defined weights
%             normalizedOptAQIValues = (optAQIValues-optAQIBestValues) .* delta_prod;
%             aqi_val = normalizedOptAQIValues' * aqi_weights;
%             % "fairness" approach weights
%         elseif aqi_weights_init == 2
%             normalizedOptAQIValues = (optAQIValues-optAQIBestValues) .* delta_prod;
%             sortedNormalizedOptAQIValues = sort(normalizedOptAQIValues,'descend');
%             aqi_val = sortedNormalizedOptAQIValues' * aqi_weights;
%         end
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %         %OLD VERSION
%         %         %multi-pollutant aggregation depends on the assumptions
%         %         %(user-defined or automatic-fairness weights)
%         %         if aqi_weights_init < 1
%         %             % user defined weights
%         %             normalizedOptAQIValues = (optAQIValues-optAQIBestValues) ./ (optAQIcle-optAQIBestValues);
%         %             aqi_val = normalizedOptAQIValues' * aqi_weights;
%         %             % "fairness" approach weights
%         %         elseif aqi_weights_init == 1
%         %             normalizedOptAQIValues = (optAQIValues-optAQIBestValues) ./ (optAQIcle-optAQIBestValues);
%         %             sortedNormalizedOptAQIValues = sort(normalizedOptAQIValues,'descend');
%         %             aqi_val = sortedNormalizedOptAQIValues' * aqi_weights;
%         %         end
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function aqi_val=MAINcompute_aqi_mp_(x,optAQIBestValues,aggregationInfo,commonDataInfo, periodIndex)
        
        optAQIValues = zeros(commonDataInfo.optim_flags.optAQINum,1);
        
        for ii = 1:commonDataInfo.optim_flags.optAQINum,
            % 20140403 ET - Added bcSuperSet(ii) input
            D=interface_get_aqi_D(aggregationInfo.geometryIntermediateData, ii);
            bc=interface_get_aqi_bc(aggregationInfo.mathIntermediateData, ii);
            nn=interface_get_aqi_nnOptSet(aggregationInfo.mathIntermediateData, ii);
            
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint21_new', 'D','bc');

            optAQIValues(ii,1) = MAINcompute_aqi_(x,...
                nn,...
                D,...
                commonDataInfo.aqi_obj_function_type(ii),...
                commonDataInfo.cell_threshold_set(commonDataInfo.aqi_obj(ii)+1),...
                commonDataInfo.aqi_obj(ii),...
                0,...
                bc, ...
                aggregationInfo, ...
                commonDataInfo, ...
                periodIndex,...
                ii);
            %%%%%%%%%%%%
            %optAQIValues(ii,1) = MAINcompute_aqi(x,nnOptSet(ii),DOptSet(ii),...
            %    commonDataInfo.aqi_obj_function_type(ii),...
            %    commonDataInfo.cell_threshold_set(commonDataInfo.aqi_obj(ii)+1),...
            %    aqi_obj(ii),...
            %    0,bcOptSet(ii));
            
            optAQIcle(ii,1)=sol_CLE.AQIs(commonDataInfo.aqi_horizon(ii),commonDataInfo.aqi_obj(ii)+1,commonDataInfo.aqi_obj_function_type(ii)+1);
            
        end
        
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint22_new', 'optAQIcle', 'optAQIBestValues');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %NEW VERSION
        %20130430 - testing the claudio proposal to improve objective
        %function computation
        delta_init=optAQIcle-optAQIBestValues;
        delta_prod=ones(commonDataInfo.optim_flags.optAQINum,1);
        for ii=1:commonDataInfo.optim_flags.optAQINum
            for kk=1:commonDataInfo.optim_flags.optAQINum
                if kk~=ii
                    delta_prod(ii,1)=delta_prod(ii,1).*delta_init(kk,1);
                end
            end
        end
 
        aqi_weights_init=commonDataInfo.aqi_weights_init;
        aqi_weights=commonDataInfo.aqi_weights;
        if aqi_weights_init <= 1
            % user defined weights
            normalizedOptAQIValues = (optAQIValues-optAQIBestValues) .* delta_prod;
            aqi_val = normalizedOptAQIValues' * aqi_weights;
            % "fairness" approach weights
        elseif aqi_weights_init == 2
            normalizedOptAQIValues = (optAQIValues-optAQIBestValues) .* delta_prod;
            sortedNormalizedOptAQIValues = sort(normalizedOptAQIValues,'descend');
            aqi_val = sortedNormalizedOptAQIValues' * aqi_weights;
        end
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint_23_new', 'aqi_weights', 'delta_prod', 'aqi_weights_init', 'aqi_val', 'delta_prod');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         %OLD VERSION
        %         %multi-pollutant aggregation depends on the assumptions
        %         %(user-defined or automatic-fairness weights)
        %         if aqi_weights_init < 1
        %             % user defined weights
        %             normalizedOptAQIValues = (optAQIValues-optAQIBestValues) ./ (optAQIcle-optAQIBestValues);
        %             aqi_val = normalizedOptAQIValues' * aqi_weights;
        %             % "fairness" approach weights
        %         elseif aqi_weights_init == 1
        %             normalizedOptAQIValues = (optAQIValues-optAQIBestValues) ./ (optAQIcle-optAQIBestValues);
        %             sortedNormalizedOptAQIValues = sort(normalizedOptAQIValues,'descend');
        %             aqi_val = sortedNormalizedOptAQIValues' * aqi_weights;
        %         end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function x_sol = MAINcompute_sol_(A,B,LB,UB,x0, aggregationInfo, periodIndex)
    function x_sol = MAINcompute_sol_(A,B,LB,UB,x0, aggregationInfo, commonDataInfo, periodIndex)
        
        %nnOptSet=interface_get_nnOptSet(mathInfo);
        %DOptSet=interface_get_DOptSet(mathInfo);
        if(commonDataInfo.optim_flags.optAQINum == 1)
            %             a=MAINcompute_aqi_(x,D,...
            %                 aqi_obj_function_type(1),...
            %                 cell_threshold_set(aqi_obj(1)+1),...
            %                 aqi_obj(1),...
            %                 0,...
            %                 aggregationInfo, ...
            %                 periodIndex,...
            %                 1);
            % 20160421 First guess setting
            D=FG_get_aqi_D(aggregationInfo.geometryIntermediateData, 1);
            % 20160421 quadrant/NN setting
            %D=interface_get_aqi_D(aggregationInfo.geometryIntermediateData, 1);
            %cyclebc, aggInfo, commonDataInfo, ii, jj
            %D=interface_get_aqi_D(aggregationInfo.geometryIntermediateData, 1);
            %bc: delta or abs values
            % 20160421 First guess setting
            bc=1;
            %bc: delta or abs values
            % 20160421 quadrant/NN setting
            %bc=interface_get_aqi_bc(aggregationInfo.mathIntermediateData, 1);
            % 20160421 First guess setting
            nn=1;
            % 20160421 quadrant/NN setting
            %nn=interface_get_aqi_nnOptSet(aggregationInfo.mathIntermediateData, 1);

            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint18a_new', 'D','bc');
            x_sol=fmincon(@(x) MAINcompute_aqi_(x,nn,D,...
                commonDataInfo.aqi_obj_function_type(1),...
                commonDataInfo.cell_threshold_set(commonDataInfo.aqi_obj(1)+1),...
                commonDataInfo.aqi_obj(1),...
                0,...
                bc, ...
                aggregationInfo, ...
                commonDataInfo, ...
                periodIndex,...
                1),...
                x0, A, B, [], [], LB, UB, [], opt);
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\var_x_sol_new', 'x_sol');
        else
            
            optAQIBestValues = zeros(commonDataInfo.optim_flags.optAQINum,1);
            
            for ii = 1:commonDataInfo.optim_flags.optAQINum,
                
                %D=aggregationInfo.mathIntermediateData.DOptSet(ii);
                % 20160421 quadrant/NN setting
                D=interface_get_aqi_D(aggregationInfo.geometryIntermediateData, ii);
                % 20160421 First guess setting
                bc=1;
                % 20160421 quadrant/NN setting
                %bc=interface_get_aqi_bc(aggregationInfo.mathIntermediateData, ii);
                % 20160421 First guess setting
                nn=1;
                % 20160421 quadrant/NN setting
                %nn=interface_get_aqi_nnOptSet(aggregationInfo.mathIntermediateData, ii);
                %nn=interface_get_aqi_nn(aggregationInfo.mathIntermediateData, ii);
                %[Da, Dpa]=interface_get_Dd_DOptSet(aggregationInfo.mathIntermediateData, ii)
                %disp(size(D));
                %disp(size(Da));
                %disp(size(Dpa));
                
                %x, D, obj_function_type, cell_threshold, AQIId, REAL, BCset, aggregationInfo, commonDataInfo, periodIndex, aqiIndex
                %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint18b_new', 'D','bc');
                x_aqi=fmincon(@(x) MAINcompute_aqi_(x,...
                    nn,...
                    D,...
                    commonDataInfo.aqi_obj_function_type(ii),...
                    commonDataInfo.cell_threshold_set(commonDataInfo.aqi_obj(ii)+1),...
                    commonDataInfo.aqi_obj(ii),...
                    0,...
                    bc, ...
                    aggregationInfo, ...
                    commonDataInfo, ...
                    periodIndex,...
                    ii),...
                    x0, A, B, [], [], LB, UB, [], opt);
                %
                %aqi_obj(ii),...
                %    0,bcOptSet(ii)
                %
                %MM check missing parameters here!
                %                 optAQIBestValues(ii,1) = MAINcompute_aqi_(x_aqi,nnOptSet(ii),DOptSet(ii),...
                %                     aqi_obj_function_type(ii),...
                %                     cell_threshold_set(aqi_obj(ii)+1),...
                %                     aqi_obj(ii),....
                %                     0);
                %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\var_x_aqi_new', 'x_aqi');
                optAQIBestValues(ii,1) = MAINcompute_aqi_(x_aqi,...
                    nn,...
                    D,...
                    commonDataInfo.aqi_obj_function_type(ii),...
                    commonDataInfo.cell_threshold_set(commonDataInfo.aqi_obj(ii)+1),...
                    commonDataInfo.aqi_obj(ii),....
                    0,...
                    bc, ...
                    aggregationInfo, ...
                    commonDataInfo, ...
                    periodIndex,...
                    ii);
            end
            
            opt=optimset(opt,'TolFun',1e-6,'MaxIter',70);
            
            %x_sol=fmincon(@(x) MAINcompute_aqi_mp(x,optAQIBestValues),...
            %    x0, A, B, [], [], LB, UB, [], opt);
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\var_optAQIBestValues_new','optAQIBestValues');
            x_sol=fmincon(@(x) MAINcompute_aqi_mp_(x,optAQIBestValues, aggregationInfo, commonDataInfo,...
                periodIndex),...
                x0, A, B, [], [], LB, UB, [], opt);
            
            opt=optimset(opt,'TolFun',commonDataInfo.optim_flags.conv_value,'MaxIter',400);
            
        end
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varpoint23_new', 'x_sol','opt');
        
    end


%     function [spi1,spi2,spi3,spi4]=MAINquadrant_emi(valori_emi,r,nx,ny)
%
%         %compute quadrant emissions
%
%         %domain dimensions and other variables
%         dimx=nx;
%         dimy=ny;
%         prova_emi=valori_emi(:,1);
%         emi=reshape(prova_emi,dimy,dimx);
%
%         %enlarged emission, with increased zeros all around (used to be able to
%         %perform simple products among emi and factor matrix
%         emi_enlarged=zeros(size(emi,1)+2*(r-2),size(emi,2)+2*(r-2));
%         emi_enlarged(r-2+1:r-2+1+size(emi,1)-1,r-2+1:r-2+1+size(emi,2)-1)=emi;
%         [sa sb]=size(emi_enlarged);
%
%         %
%         %left quadrant
%         %factors to perform products of cell emissions
%         fattore=tril(ones(r+1,r+1))-diag(repmat(0.5,r+1,1));
%         sotto=flipud(fattore);
%         fattore=[fattore;sotto(2:end,:)];
%         fattore(r+1,r+1)=0.25;
%
%         %initialize variable, and set dimensions
%         spicchio3=zeros(dimy,dimx);
%
%         %down quadrant
%         %factors to perform products of cell emissions
%         fattoreD=rot90(fattore,3);
%
%         %initialize variable, and set dimensions
%         spicchio1=zeros(dimy,dimx);
%
%         %right quadrant
%         %factors to perform products of cell emissions
%         fattoreR=rot90(fattore,2);
%
%         %initialize variable, and set dimensions
%         spicchio4=zeros(dimy,dimx);
%
%         %up quadrant
%         %factors to perform products of cell emissions
%         fattoreU=rot90(fattore,1);
%
%         %initialize variable, and set dimensions
%         spicchio2=zeros(dimy,dimx);
%
%         %loop to aggregate emissions
%         for i=r+1:sa-r
%             for j=r+1:sb-r
%                 %left quadrant
%                 spicchio3(i-r+2,j-r+2)=sum(sum(emi_enlarged(i-r:i+r,j-r:j).*fattore));
%                 %up quadrant (down as matrix, up geographically speaking)
%                 spicchio1(i-r+2,j-r+2)=sum(sum(emi_enlarged(i-r:i,j-r:j+r).*fattoreD));
%                 %right quadrant
%                 spicchio4(i-r+2,j-r+2)=sum(sum(emi_enlarged(i-r:i+r,j:j+r).*fattoreR));
%                 %down quadrant (up as matrix, down geographically speaking)
%                 spicchio2(i-r+2,j-r+2)=sum(sum(emi_enlarged(i:i+r,j-r:j+r).*fattoreU));
%             end
%         end
%
%         spi1=reshape(spicchio1,dimx*dimy,1);
%         spi2=reshape(spicchio2,dimx*dimy,1);
%         spi3=reshape(spicchio3,dimx*dimy,1);
%         spi4=reshape(spicchio4,dimx*dimy,1);
%
%     end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTING AGGREGATED SCENARIO MODE SOLUTIONS
%     function [aqi_per_cell]=MAINaggregated_scenario_mode(emiTMP,aggregationInfo,ii,jj)
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %ANNS CONSIDERED
%         %NN=nnSuperSet(ii).nnSet(jj);
%         %icells=NN.icells;
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %AREAL CASE
%         % original version
%         [s1_NOX,s2_NOX,s3_NOX,s4_NOX]=MAINquadrant_emi(emiTMP(:,1),icells,nx,ny);
%         [s1_VOC,s2_VOC,s3_VOC,s4_VOC]=MAINquadrant_emi(emiTMP(:,2),icells,nx,ny);
%         [s1_NH3,s2_NH3,s3_NH3,s4_NH3]=MAINquadrant_emi(emiTMP(:,3),icells,nx,ny);
%         [s1_PM10,s2_PM10,s3_PM10,s4_PM10]=MAINquadrant_emi(emiTMP(:,4),icells,nx,ny);
%         [s1_PM25,s2_PM25,s3_PM25,s4_PM25]=MAINquadrant_emi(emiTMP(:,5),icells,nx,ny);
%         [s1_SO2,s2_SO2,s3_SO2,s4_SO2]=MAINquadrant_emi(emiTMP(:,6),icells,nx,ny);
%         NH3_allOld=[s1_NH3,s2_NH3,s3_NH3,s4_NH3];
%         NOX_allOld=[s1_NOX,s2_NOX,s3_NOX,s4_NOX];
%         PM10_allOld=[s1_PM10,s2_PM10,s3_PM10,s4_PM10];
%         PM25_allOld=[s1_PM25,s2_PM25,s3_PM25,s4_PM25];
%         SO2_allOld=[s1_SO2,s2_SO2,s3_SO2,s4_SO2];
%         VOC_allOld=[s1_VOC,s2_VOC,s3_VOC,s4_VOC];
%         % original version end
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %CONSIDER CASE WITH 6 OR 2 PERCURSOR EMISSIONS AS INPUT,
%         %AREAL CASE. ALWAYS THERE ARE 4 QUADRANTS TO CONSIDER WIND
%         %DIRECTIONS
%         if optim_flags.areal_point==0
%             if (size(NN.net,1)==24 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==24))% ANNlinear, 6 input
%                 emissioni=[NH3_all,NOX_all,PM10_all,PM25_all,SO2_all,VOC_all];
%             elseif (size(NN.net,1)==8 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==8))%ANNlinear, 2 input
%                 emissioni=[NOX_all,VOC_all];
%             end
%
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %IF ALSO POINT SOURCES
%         elseif optim_flags.areal_point==1
%             %point
%             %original version
%             [s1_NOXp,s2_NOXp,s3_NOXp,s4_NOXp]=MAINquadrant_emi(emiTMP(:,7),icells,nx,ny);
%             [s1_VOCp,s2_VOCp,s3_VOCp,s4_VOCp]=MAINquadrant_emi(emiTMP(:,8),icells,nx,ny);
%             [s1_NH3p,s2_NH3p,s3_NH3p,s4_NH3p]=MAINquadrant_emi(emiTMP(:,9),icells,nx,ny);
%             [s1_PM10p,s2_PM10p,s3_PM10p,s4_PM10p]=MAINquadrant_emi(emiTMP(:,10),icells,nx,ny);
%             [s1_PM25p,s2_PM25p,s3_PM25p,s4_PM25p]=MAINquadrant_emi(emiTMP(:,11),icells,nx,ny);
%             [s1_SO2p,s2_SO2p,s3_SO2p,s4_SO2p]=MAINquadrant_emi(emiTMP(:,12),icells,nx,ny);
%
%             NH3_allpOld=[s1_NH3p,s2_NH3p,s3_NH3p,s4_NH3p];
%             NOX_allpOld=[s1_NOXp,s2_NOXp,s3_NOXp,s4_NOXp];
%             PM10_allpOld=[s1_PM10p,s2_PM10p,s3_PM10p,s4_PM10p];
%             PM25_allpOld=[s1_PM25p,s2_PM25p,s3_PM25p,s4_PM25p];
%             SO2_allpOld=[s1_SO2p,s2_SO2p,s3_SO2p,s4_SO2p];
%             VOC_allpOld=[s1_VOCp,s2_VOCp,s3_VOCp,s4_VOCp];
%             % end
%
%             if isequal(NH3_allpOld, NH3_allp)
%                 disp('Right Yo!!!!');
%             else disp('Mamma miaaaa!!!!');
%             end
%             if isequal(PM25_allpOld, PM25_allp)
%                 disp('Right Yo!!!!');
%             else disp('Mamma miaaaa!!!!');
%             end
%
%             %create input structure, to be used in the ANNs
%             if (size(NN.net,1)==48 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==48))% ANN/linear, 6 input
%                 emissioni=[NH3_all,NOX_all,PM10_all,PM25_all,SO2_all,VOC_all,...
%                     NH3_allp,NOX_allp,PM10_allp,PM25_allp,SO2_allp,VOC_allp];
%             elseif (size(NN.net,1)==16 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==16))%ANNlinear, 2 input
%                 emissioni=[NOX_all,VOC_all,NOX_allp,VOC_allp];
%             end
%
%         end
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %ANN INPUT DATA
%         %keep only optim domain
%         emissioni(find(domainInfo.flag_optim_dom==0),:)=[];
%         %20130820 - consider only if cell completerly in PAD
%         %         emissioni(find(flag_optim_dom==0 | flag_optim_dom==2),:)=[];
%         input_rete2=emissioni';
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %IN CASE QUADRANT EMISSIONS TOO CLOSE TO CTM BOUNDARY
%         if strcmp(dirs.pathAR,'-1')==0
%             AR=load(dirs.pathAR);                             %load not in optim
%             AR=AR.Ratio;                                 %rename variable
%             AR(domainInfo.flag_optim_dom==0,:)=[];                  %remove cells outside flag_optim_dom
%             if (size(NN.net,1)==48 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==48))% ANN/linear, 6 input
%                 ARreordrepm=repmat(AR,1,12);            %repmat
%             elseif (size(NN.net,1)==16 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==16))%ANNlinear, 2 input
%                 ARreordrepm=repmat(AR,1,4);            %repmat
%             end
%             input_rete2=emissioni'./ARreordrepm';        %rewrite input_rete2 dividing emissions by area_ratio
%         end
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %CHECK IF EMISSIONS ARE INSIDE ANNS IDENTIFICATION
%         %BOUNDS...IF NOT, STOP RIAT+
%         %compare ANNs input and ps_input, to see if we are in the
%         %ANNs bounds
%         %check in the PAD
%         %         for indcel=1:ncel
%         %                 if sum(input_rete2(:,indcel)<NN.ps_input.xmin)>0 | sum(input_rete2(:,indcel)>NN.ps_input.xmax)>0
%         %                     find(input_rete2(:,indcel)<NN.ps_input.xmin)
%         %                     find(input_rete2(:,indcel)>NN.ps_input.xmax)
%         %                     figure;plot([input_rete2(:,indcel)-NN.ps_input.xmin  NN.ps_input.xmax-input_rete2(:,indcel)])
%         %                     error('SR model identification bounds not respected...RIAT+ is terminated')
%         %                 end
%         %         end
%         %matrix of min and max ANNS values
%         minVal=repmat(NN.ps_input.xmin,1,size(emissioni,1));
%         maxVal=repmat(NN.ps_input.xmax,1,size(emissioni,1));
%         %check if my scenario is out of bounds
% %         indmin=find(input_rete2<minVal);  OK
% %         indmax=find(input_rete2>maxVal);  OK
%         %20130828 - input_rete2>0 added to manage issue of the Alsace cells
%         %adjacent to the CTM grid
%         indmin=find(input_rete2<minVal & input_rete2>0);
%         indmax=find(input_rete2>maxVal & input_rete2>0);
%         %check bounds problem.....if limited problems, simulate the
%         %scenario, otherwise exit
% %         ii
% %         jj
% %         max((input_rete2(indmax)-maxVal(indmax))./(maxVal(indmax))*100)
% %         min((input_rete2(indmin)-minVal(indmin))./minVal(indmin)*100)
%
%         if length(indmax)>0
%             %             if max(input_rete2(indmax)-maxVal(indmax))>1
%             if max((input_rete2(indmax)-maxVal(indmax))./(maxVal(indmax))*100)>30
%
%                 Error=3;
%
%                 fprintf(fidExit, int2str(Error));
%                 strStatus='PROGRESSION: Aggregated scenario analysis finished.';
%                 disp(strStatus);
%                 fprintf(fidStatus, '%s\n',strStatus);
%                 error('SR model identification bounds not respected...RIAT+ is terminated')
%             end
%         end
%
%         if length(indmin)>0
%             %             if max(input_rete2(indmin)-minVal(indmin))<-1
%             if min((input_rete2(indmin)-minVal(indmin))./minVal(indmin)*100)<-30
%
%                 Error=3;
%                 fprintf(fidExit, int2str(Error));
%                 strStatus='PROGRESSION: Aggregated scenario analysis finished.';
%                 disp(strStatus);
%                 fprintf(fidStatus, '%s\n',strStatus);
%                 error('SR model identification bounds not respected...RIAT+ is terminated')
%             end
%         end
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %IF BOUNDS RESPECTED, RUN SR MODEL
%         %if constraints respected for this AQI, create results
%         %define if lin o net
%         if size(NN.net,1)==1
%             optim_flags.reg_net=1;
%         elseif size(NN.net,1)>1
%             optim_flags.reg_net=0;
%         end
%
%         switch optim_flags.reg_net
%             case 0 %linear case
%                 aqi_per_cell=(input_rete2'*NN.net)';
%             case 1 %nonlinear case
%                 %20111215 - run the network only on the optimization cells
%                 NN.ps_input.no_change=0; %to be compatible among different matlab versions
%                 [input_rete_norm2]=mapminmax('apply',input_rete2,NN.ps_input);
%
%                 %check for results outside bounds - force to be inside bounds
%                 input_rete_norm2(find(input_rete_norm2<-1))=-1;
%                 input_rete_norm2(find(input_rete_norm2>1))=1;
%
%                 output_rete_norm2=sim(NN.net,input_rete_norm2);
%
%                 %%%
%                 %check for results outside bounds - force to be inside bounds
%                 output_rete_norm2(find(output_rete_norm2<-1))=-1;
%                 output_rete_norm2(find(output_rete_norm2>1))=1;
%                 %%%
%
%                 NN.ps_target.no_change=0; %to be compatible among different matlab versions
%                 aqi_per_cell=mapminmax('reverse',output_rete_norm2,NN.ps_target);
%         end
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%keep track of computing time
toc

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTING AGGREGATED SCENARIO MODE SOLUTIONS
function [aqi_per_cell]=MAINaggregated_scenario_mode_(emiTMP,NN,periodIndex,aqiIndex, aggregationInfo, commonDataInfo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANNS CONSIDERED
%check how to pass parameter (ii and jj)
%NN=interface_get_nnSuperSet_indexed(aggregationInfo.mathIntermediateData, ii, jj);
%icells=aggregationInfo.commonDataInfo.totCells;
%NN=nnSuperSet(ii).nnSet(jj);
%icells=NN.icells;

coordinate=commonDataInfo.domainData.data(:,1:2);
stepsize=round(max(abs(coordinate(2,1)-coordinate(1,1)),abs(coordinate(2,2)-coordinate(1,2))));
nxny=round((max(coordinate)-min(coordinate))/stepsize+1);
nx=nxny(1,1);
ny=nxny(1,2);
ncel=nx*ny;

orderedPolls = {'NOX';'VOC';'NH3';'PM10';'PM25';'SO2'}
np=size(orderedPolls);
quad=4;
indicators=0;
x=0;
y=0;

%flag_region_dom=commonDataInfo.domainData.data(:,3);
indopt=flag_region_dom; %cells in optimization domain
indoptrep=repmat(indopt,quad*np,1); %index repeated for pollutants and quadrants
indoptrep1=repmat(indopt,quad,1); %index repeated for pollutants and quadrants
indiciMAT=reshape(indopt,ny,nx);
splitResult=1;

% fill details for this case (pathann, radius...)
sixInput=interface_is6Input(aggregationInfo.mathIntermediateData,periodIndex,aqiIndex);
twoInput=interface_is2Input(aggregationInfo.mathIntermediateData,periodIndex,aqiIndex);
commonDataInfo=fill_commonDataInfo(commonDataInfo, periodIndex, aqiIndex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AREAL CASE
% MM new version

thisPrecursor=emiTMP(:,1);
[thisIntermediateResult, thisResult]=do_Job(aggregationInfo.geometryDataInfo, commonDataInfo, thisPrecursor,...
    indicators, x, y, nx, ny, ...
    orderedPolls(1), indiciMAT, splitResult);
NOX_all=thisResult.resGrid;

thisPrecursor=emiTMP(:,2);
[thisIntermediateResult, thisResult]=do_Job(aggregationInfo.geometryDataInfo, commonDataInfo, thisPrecursor,...
    indicators, x, y, nx, ny, ...
    orderedPolls(2), indiciMAT, splitResult);
VOC_all=thisResult.resGrid;

thisPrecursors=emiTMP(:,3);
[thisIntermediateResult, thisResult]=do_Job(aggregationInfo.geometryDataInfo, commonDataInfo, thisPrecursor,...
    indicators, x, y, nx, ny, ...
    orderedPolls(3), indiciMAT, splitResult);
NH3_all=thisResult.resGrid;

thisPrecursors=emiTMP(:,4);
[thisIntermediateResult, thisResult]=do_Job(aggregationInfo.geometryDataInfo, commonDataInfo, thisPrecursor,...
    indicators, x, y, nx, ny, ...
    orderedPolls(4), indiciMAT, splitResult);
PM10_all=thisResult.resGrid;

thisPrecursors=emiTMP(:,5);
[thisIntermediateResult, thisResult]=do_Job(aggregationInfo.geometryDataInfo, commonDataInfo, thisPrecursor,...
    indicators, x, y, nx, ny, ...
    orderedPolls(5), indiciMAT, splitResult);
PM25_all=thisResult.resGrid;

thisPrecursors=emiTMP(:,6);
[thisIntermediateResult, thisResult]=do_Job(aggregationInfo.geometryDataInfo, commonDataInfo, thisPrecursor,...
    indicators, x, y, nx, ny, ...
    orderedPolls(6), indiciMAT, splitResult);
SO2_all=thisResult.resGrid;

% MM end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONSIDER CASE WITH 6 OR 2 PERCURSOR EMISSIONS AS INPUT,
%AREAL CASE. ALWAYS THERE ARE 4 QUADRANTS TO CONSIDER WIND
%DIRECTIONS
if commonDataInfo.optim_flags.areal_point==0
    
    if sixInput% ANN/linear, 6 input
        emissioni=[NH3_all,NOX_all,PM10_all,PM25_all,SO2_all,VOC_all];
        
        %elseif (size(NN.net,1)==16 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==16))%ANNlinear, 2 input
    elseif twoInput%ANNlinear, 2 input
        emissioni=[NOX_all,VOC_all];
        
    end
    
    %if (size(NN.net,1)==24 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==24))% ANNlinear, 6 input
    %    emissioni=[NH3_all,NOX_all,PM10_all,PM25_all,SO2_all,VOC_all];
    %elseif (size(NN.net,1)==8 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==8))%ANNlinear, 2 input
    %    emissioni=[NOX_all,VOC_all];
    %end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %IF ALSO POINT SOURCES
elseif commonDataInfo.optim_flags.areal_point==1
    %point
    % MM new version
    
    orderedPolls={'NOX_d';'VOC_d';'NH3_d';'PM10_d';'PM25_d';'SO2_d'};
    
    thisPrecursors=emiTMP(:,7);
    [thisIntermediateResult, thisResult]=do_Job(aggregationInfo.geometryDataInfo, commonDataInfo, thisPrecursor,...
        indicators, x, y, nx, ny, ...
        orderedPolls(1), indiciMAT, splitResult);
    NOX_allp=thisResult.resGrid;
    
    thisPrecursors=emiTMP(:,8);
    [thisIntermediateResult, thisResult]=do_Job(aggregationInfo.geometryDataInfo, commonDataInfo, thisPrecursor,...
        indicators, x, y, nx, ny, ...
        orderedPolls(2), indiciMAT, splitResult);
    VOC_allp=thisResult.resGrid;
    
    thisPrecursors=emiTMP(:,9);
    [thisIntermediateResult, thisResult]=do_Job(aggregationInfo.geometryDataInfo, commonDataInfo, thisPrecursor,...
        indicators, x, y, nx, ny, ...
        orderedPolls(3), indiciMAT, splitResult);
    NH3_allp=thisResult.resGrid;
    
    thisPrecursors=emiTMP(:,10);
    [thisIntermediateResult, thisResult]=do_Job(aggregationInfo.geometryDataInfo, commonDataInfo, thisPrecursor,...
        indicators, x, y, nx, ny, ...
        orderedPolls(4), indiciMAT, splitResult);
    PM10_allp=thisResult.resGrid;
    
    thisPrecursors=emiTMP(:,11);
    [thisIntermediateResult, thisResult]=do_Job(aggregationInfo.geometryDataInfo, commonDataInfo, thisPrecursor,...
        indicators, x, y, nx, ny, ...
        orderedPolls(5), indiciMAT, splitResult);
    PM25_allp=thisResult.resGrid;
    
    thisPrecursors=emiTMP(:,12);
    [thisIntermediateResult, thisResult]=do_Job(aggregationInfo.geometryDataInfo, commonDataInfo, thisPrecursor,...
        indicators, x, y, nx, ny, ...
        orderedPolls(6), indiciMAT, splitResult);
    SO2_allp=thisResult.resGrid;
    
    % MM end
    %create input structure, to be used in the ANNs
    if sixInput% ANN/linear, 6 input
        emissioni=[NH3_all,NOX_all,PM10_all,PM25_all,SO2_all,VOC_all,...
            NH3_allp,NOX_allp,PM10_allp,PM25_allp,SO2_allp,VOC_allp];
    elseif twoInput% ANN/linear, 2 input
        emissioni=[NOX_all,VOC_all,NOX_allp,VOC_allp];
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANN INPUT DATA
%keep only optim domain
% already done???
% emissioni dims = 10496X48
emissioni(find(commonDataInfo.domainInfo.flag_optim_dom==0),:)=[];
%emissioni(find(indoptrep1==0),:)=[];
%aqi_per_cell=do_aqi_per_cell(emissioni, aggregationInfo);
aqi_per_cell=do_interface_get_aqi_per_cell(emissioni, NN, aggregationInfo, commonDataInfo, periodIndex, aqiIndex);

%         %20130820 - consider only if cell completerly in PAD
%         %         emissioni(find(flag_optim_dom==0 | flag_optim_dom==2),:)=[];
%         input_rete2=emissioni';
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %IN CASE QUADRANT EMISSIONS TOO CLOSE TO CTM BOUNDARY
%         if strcmp(dirs.pathAR,'-1')==0
%             AR=load(dirs.pathAR);                             %load not in optim
%             AR=AR.Ratio;                                 %rename variable
%             AR(domainInfo.flag_optim_dom==0,:)=[];                  %remove cells outside flag_optim_dom
%             %if (size(NN.net,1)==48 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==48))% ANN/linear, 6 input
%             if sixInput
%                 ARreordrepm=repmat(AR,1,12);            %repmat
%             %elseif (size(NN.net,1)==16 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==16))%ANNlinear, 2 input
%             elseif twoInput%ANNlinear, 2 input
%                 ARreordrepm=repmat(AR,1,4);            %repmat
%             end
%             input_rete2=emissioni'./ARreordrepm';        %rewrite input_rete2 dividing emissions by area_ratio
%         end
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %CHECK IF EMISSIONS ARE INSIDE ANNS IDENTIFICATION
%         %BOUNDS...IF NOT, STOP RIAT+
%         %compare ANNs input and ps_input, to see if we are in the
%         %ANNs bounds
%         %check in the PAD
%         %         for indcel=1:ncel
%         %                 if sum(input_rete2(:,indcel)<NN.ps_input.xmin)>0 | sum(input_rete2(:,indcel)>NN.ps_input.xmax)>0
%         %                     find(input_rete2(:,indcel)<NN.ps_input.xmin)
%         %                     find(input_rete2(:,indcel)>NN.ps_input.xmax)
%         %                     figure;plot([input_rete2(:,indcel)-NN.ps_input.xmin  NN.ps_input.xmax-input_rete2(:,indcel)])
%         %                     error('SR model identification bounds not respected...RIAT+ is terminated')
%         %                 end
%         %         end
%         %matrix of min and max ANNS values
%         xMin=interface_get_xMin(aggregation, periodIndex,aqiIndex);
%         xMax=interface_get_xMax(aggregation, periodIndex,aqiIndex);
%         %minVal=repmat(NN.ps_input.xmin,1,size(emissioni,1));
%         %maxVal=repmat(NN.ps_input.xmax,1,size(emissioni,1));
%         minVal=repmat(xmin,1,size(emissioni,1));
%         maxVal=repmat(xmax,1,size(emissioni,1));
%         %check if my scenario is out of bounds
% %         indmin=find(input_rete2<minVal);  OK
% %         indmax=find(input_rete2>maxVal);  OK
%         %20130828 - input_rete2>0 added to manage issue of the Alsace cells
%         %adjacent to the CTM grid
%         indmin=find(input_rete2<minVal & input_rete2>0);
%         indmax=find(input_rete2>maxVal & input_rete2>0);
%         %check bounds problem.....if limited problems, simulate the
%         %scenario, otherwise exit
% %         ii
% %         jj
% %         max((input_rete2(indmax)-maxVal(indmax))./(maxVal(indmax))*100)
% %         min((input_rete2(indmin)-minVal(indmin))./minVal(indmin)*100)
%
%         if length(indmax)>0
%             %             if max(input_rete2(indmax)-maxVal(indmax))>1
%             if max((input_rete2(indmax)-maxVal(indmax))./(maxVal(indmax))*100)>30
%
%                 Error=3;
%
%                 fprintf(fidExit, int2str(Error));
%                 strStatus='PROGRESSION: Aggregated scenario analysis finished.';
%                 disp(strStatus);
%                 fprintf(fidStatus, '%s\n',strStatus);
%                 error('SR model identification bounds not respected...RIAT+ is terminated')
%             end
%         end
%
%         if length(indmin)>0
%             %             if max(input_rete2(indmin)-minVal(indmin))<-1
%             if min((input_rete2(indmin)-minVal(indmin))./minVal(indmin)*100)<-30
%
%                 Error=3;
%                 fprintf(fidExit, int2str(Error));
%                 strStatus='PROGRESSION: Aggregated scenario analysis finished.';
%                 disp(strStatus);
%                 fprintf(fidStatus, '%s\n',strStatus);
%                 error('SR model identification bounds not respected...RIAT+ is terminated')
%             end
%         end
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %IF BOUNDS RESPECTED, RUN SR MODEL
%         %if constraints respected for this AQI, create results
%         %define if lin o net
%
%         aqi_per_cell=interface_get_aqi_per_cell(aggregationInfo, input_rete2)
%         %if size(NN.net,1)==1
%         %    optim_flags.reg_net=1;
%         %elseif size(NN.net,1)>1
%         %    optim_flags.reg_net=0;
%         %end
%
%         %switch optim_flags.reg_net
%         %    case 0 %linear case
%         %        aqi_per_cell=(input_rete2'*NN.net)';
%         %    case 1 %nonlinear case
%         %        %20111215 - run the network only on the optimization cells
%         %        NN.ps_input.no_change=0; %to be compatible among different matlab versions
%         %        [input_rete_norm2]=mapminmax('apply',input_rete2,NN.ps_input);
%         %
%         %        %check for results outside bounds - force to be inside bounds
%         %        input_rete_norm2(find(input_rete_norm2<-1))=-1;
%         %        input_rete_norm2(find(input_rete_norm2>1))=1;
%         %
%         %        output_rete_norm2=sim(NN.net,input_rete_norm2);
%         %
%         %        %%%
%         %        %check for results outside bounds - force to be inside bounds
%         %        output_rete_norm2(find(output_rete_norm2<-1))=-1;
%         %        output_rete_norm2(find(output_rete_norm2>1))=1;
%         %        %%%
%
%         %        NN.ps_target.no_change=0; %to be compatible among different matlab versions
%         %        aqi_per_cell=mapminmax('reverse',output_rete_norm2,NN.ps_target);
%         %end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
