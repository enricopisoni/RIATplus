function []=MAINmain(f1,f2)
%#function network

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       READING INPUT FILES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       flag_optim_path.txt

fid=fopen(f1);
C1=textscan(fid,'%s');
% variables
pathEMI=cell2mat(C1{1}(1));          % path for emission files
pathFOD=cell2mat(C1{1}(2));          % path for coordinate, policy application domain, AQI computation domain, population
pathFOO=cell2mat(C1{1}(3));          % path for optimization configuration file
pathOCM=cell2mat(C1{1}(4));          % path for measures DB
pathLST=cell2mat(C1{1}(5));          % path for file with techs list (to create NH3 constraints)
pathTRD=cell2mat(C1{1}(6));          % path for traffic duplication
pathDEF=cell2mat(C1{1}(7));          % path for AQI definition
pathPM =cell2mat(C1{1}(8));          % PM10: relation to convert "yearly average" to "# daily threshold exceedances"
pathOUF=cell2mat(C1{1}(9));          % path for output file
pathAR=cell2mat(C1{1}(10));          % path for areal_ratio variable (used for ASPA, for cells close to domain boundary)
pathADS=cell2mat(C1{1}(11));         % path for scenario analysis file, containing ms-pollutant % emission reductions
nocID=cell2mat(C1{1}(12));           % no control technology ID
AQINum=str2num(cell2mat(C1{1}(13))); % number of AQIs to consider
nocID=str2num(nocID);                % nocID to be transformed in number

fclose(fid);

%                       aqi_definition.txt

% preparing the data structure to store the AQIs configuration
% 7 columns are the 7 AQIs
% 3 structures are: 1)annual, 2)winter, 3)summer
pathANN = struct ('ANNs',{});
pathDd = struct ('Dd',{});
for h=1:3;
    path = cellstr('-999');
    for p=1:AQINum-1,
        path = [path, cellstr('-999')];
    end;
    path = char(path);
    
    pathANN = [pathANN, struct('ANNs',path)];
    pathDd = [pathDd, struct('Dd',path)];
end;

fid = fopen(pathDEF);
C1 = textscan(fid,'%s','delimiter','\n');

row = 1;
entryNum = str2num(cell2mat(C1{1}(row))); row = row + 1;

row = row + 1;
for k=1:entryNum,
    id = regexp(cell2mat(C1{1}(row)),'^[0-9]+','match');   row = row + 1;
    months = regexp(cell2mat(C1{1}(row)),'^[0-9]+','match');   row = row + 1;
    %     ann = regexp(cell2mat(C1{1}(row)),'\.*(/\w+)+','match');  row = row + 1;
    %     Dd = regexp(cell2mat(C1{1}(row)),'\.*(/\w+)+','match');  row = row + 1;
    %20130731 - generalize the use of paths - both absolute and relative
    ann=sscanf(cell2mat(C1{1}(row)),'%s %'); row = row + 1;
    Dd=sscanf(cell2mat(C1{1}(row)),'%s %'); row = row + 1;
    row = row + 1;
    
    h = 0;
    if (isequal(strtrim(char(months)),'1')==1)
        h = 1;
    end
    if (isequal(strtrim(char(months)),'2')==1)
        h = 2;
    end
    if (isequal(strtrim(char(months)),'3')==1)
        h = 3;
    end
    
    jj = cell2mat(id)-char('0') + 1;
    if(jj == 1)
        pathANN(h) = struct('ANNs', char(char(ann), pathANN(h).ANNs(2:AQINum,:)));
        pathDd(h) = struct('Dd', char(char(Dd), pathDd(h).Dd(2:AQINum,:)));
    end
    if(jj == AQINum)
        pathANN(h) = struct('ANNs', char(pathANN(h).ANNs(1:AQINum-1,:), char(ann)));
        pathDd(h) = struct('Dd', char(pathDd(h).Dd(1:AQINum-1,:), char(Dd)));
    end
    if((jj > 1) && (jj < AQINum))
        pathANN(h) = ...
            struct('ANNs', char(pathANN(h).ANNs(1:jj-1,:), ...
            char(ann), ...
            pathANN(h).ANNs(jj+1:AQINum,:)));
        pathDd(h) = ...
            struct('Dd', char(pathDd(h).Dd(1:jj-1,:), ...
            char(Dd), ...
            pathDd(h).Dd(jj+1:AQINum,:)));
    end
end;

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       flag_optim_oth.txt
%                  (optimization configuration)

flag_optim_oth = load(pathFOO);
coordinateFlip       =flag_optim_oth(1); % 0 means no data flipping is required, 1 means data flipping needed (as in RIAT).
flag_reg_net         =flag_optim_oth(2); % DUMMY VAR - NOT USED 0 means you consider linear model, 1 neural network model, for all AQIs.
flag_constraints     =flag_optim_oth(3); % flag_constraints=0 means all techs are not replaceable (keep LB as they are).
flag_mode_ce_mo      =flag_optim_oth(4); % 0 means multi-objective, 1 means cost-effectiveness, 2 means detailed SA, 3 means aggregated SA.
flag_paretopoints    =flag_optim_oth(5); % If 0, compute 5 points of pareto curve. If greataer than 0, it is the execution time to be dedicated to pareto points computation.
cell_threshold_set = zeros(AQINum,1);
cell_threshold_set(1,1) = flag_optim_oth(6);
cell_threshold_set(2,1) = flag_optim_oth(7);
cell_threshold_set(3,1) = flag_optim_oth(8);
cell_threshold_set(4,1) = flag_optim_oth(9);
cell_threshold_set(5,1) = flag_optim_oth(10);
cell_threshold_set(6,1) = flag_optim_oth(11);
cell_threshold_set(7,1) = flag_optim_oth(12);
conv_value           =flag_optim_oth(13); % Convergence value for the optimization algorithm.
areal_point          =flag_optim_oth(14); % 0 means areal and point emissions summed, 1 means areal and point emissions separated.
abs_perc             =flag_optim_oth(15); % Specify if the cost constraints expressed in absolute (0) or MRR% values (1). If cost-effectiveness, put 0.
tmp_thres_cost(1)    =flag_optim_oth(16); % Costs (Meuro over CLE) for Pareto curve (CLE and MRR automatically created). If cost-effectiveness, only this value is considered.
tmp_thres_cost(2)    =flag_optim_oth(17); % Costs (Meuro over CLE) for Pareto curve (CLE and MRR automatically created). If cost-effectiveness, this value is not considered.
tmp_thres_cost(3)    =flag_optim_oth(18); % Costs (Meuro over CLE) for Pareto curve (CLE and MRR automatically created). If cost-effectiveness, this value is not considered.
optAQINum            =flag_optim_oth(19); % Number of AQIs to consider during the optimization.

offset = 19; % The number of the last row read.

% Read other infos for optimization, as:
% - number o
% AQIs to be optimized;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FROM HERE ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create directory name
mkdir('output')

if (flag_mode_ce_mo==0 | flag_mode_ce_mo==1 | flag_mode_ce_mo==2)
    %load gains DB
    %load file "out_clemfr.txt", containing informations all quadrupes corinair
    %macrosector-sectorIiasa-activityIiasa-tech, and related informations
    % global_data=load('input/out_clemfr.txt');
    global_data=load(pathOCM);
    
    %extract part of gains DB related to emissions for which no reduction
    %technologies are available
    global_data_noc=global_data(find(global_data(:,4)==nocID),:);
    global_data(find(global_data(:,4)==nocID),:)=[];
    
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
domainData=importdata(pathFOD);
if size(domainData.data,2)==6
    flag_region_dom=domainData.data(:,3); %regional domain
    flag_optim_dom=domainData.data(:,4);%optimization domain
    flag_aqi_dom=domainData.data(:,5); %aqi computation domain
    pop=domainData.data(:,6); %population
else
    flag_region_dom=domainData.data(:,3);
    flag_optim_dom=domainData.data(:,3);%optimization domain
    flag_aqi_dom=domainData.data(:,4); %aqi computation domain
    pop=domainData.data(:,5); %population
end


% load domain coordinate
coordinate=domainData.data(:,1:2);
xutm=coordinate(:,1);
yutm=coordinate(:,2);

%compute nx and ny from coordinate
% stepsize=max(coordinate(2,1)-coordinate(1,1),coordinate(2,2)-coordinate(1,2));
stepsize=round(max(abs(coordinate(2,1)-coordinate(1,1)),abs(coordinate(2,2)-coordinate(1,2))));

nxny=round((max(coordinate)-min(coordinate))/stepsize+1);
nx=nxny(1,1);
ny=nxny(1,2);
ncel=nx*ny;

%load cells emissions
%emi: emissions of the starting case (2 rows in the case of cells
%outside domain, sector-activity-techs rows if inside domain)
%if the cell is inside the domain, the correspondent file contains both
%emissions (first 6 columns) and activity level (seventh column)
%emi contains two columnes: column 1 is the aggregated data, column 2 the
%detailed data

if (flag_mode_ce_mo==0 | flag_mode_ce_mo==1 | flag_mode_ce_mo==2)
    pathEMI=strcat(pathEMI,'TP1/');
end
%UNIBS(ET)20131001 - flag_optim_dom changed with flag_region_dom in input
[emi,flag_ADS_tp]=MAINload_emi(ncel,pathEMI,flag_region_dom,flag_mode_ce_mo);

%in case of RIAT, matrices and input are cut and flipped (due to an initial
%problem in the data). Otherwise no change to the data is performed
if coordinateFlip==1
    ny=ny-3; %problem in RIAT about number of cells
    ncel=nx*ny;
    [dataerase]=find(coordinate(:,1)>874.28 | coordinate(:,2)>5210.7);
    flag_optim_dom(dataerase,:)=[];
    coordinate(dataerase,:)=[];
    pop(dataerase,:)=[];
    flag_aqi_dom(dataerase,:)=[];
    %UNIBS(ET)20131001 - correct regional flag
    flag_region_dom(dataerase,:)=[];
    
    %redefine xutm and yutm
    xutm=coordinate(:,1);
    yutm=coordinate(:,2);
    %process files to get the following format:
    % xutm1 yutm1
    % xutm1 yutm2
    % ...
    %while now format is
    % xutm1 yutm1
    % xutm2 yutm1
    % ...
    %new optim flag
    flag_optim_dom1=reshape(flag_optim_dom,nx,ny)';
    flag_optim_dom2=reshape(flag_optim_dom1,nx*ny,1);
    flag_optim_dom=flag_optim_dom2;
    
    %UNIBS(ET)20131001 - new regional flag
    flag_region_dom1=reshape(flag_region_dom,nx,ny)';
    flag_region_dom2=reshape(flag_region_dom1,nx*ny,1);
    flag_region_dom=flag_region_dom2;
    
    %new pop file
    pop1=reshape(pop,nx,ny)';
    pop2=reshape(pop1,nx*ny,1);
    pop=pop2;
    
    %new aqi compuation domain file
    flag_aqi_dom1=reshape(flag_aqi_dom,nx,ny)';
    flag_aqi_dom2=reshape(flag_aqi_dom1,nx*ny,1);
    flag_aqi_dom=flag_aqi_dom2;
    
    %new xutm
    xutm1=reshape(xutm,nx,ny)';
    xutm2=reshape(xutm1,nx*ny,1);
    
    %new xutm
    yutm1=reshape(yutm,nx,ny)';
    yutm2=reshape(yutm1,nx*ny,1);
    
    %new coordinates
    coordinate=[xutm2 yutm2];
    xutm=xutm2;
    yutm=yutm2;
    %clipping data
    emi(dataerase,:)=[];
    %taking into account both aggregated and detailed emissions - this part
    %has been added to manage "boundary cells" (cells partly inside and
    %partly outside the optimization domain)
    emiAGG=emi(:,1);
    emiDET=emi(:,2);
    emiAGG1=reshape(emiAGG,nx,ny)';
    emiAGG2=reshape(emiAGG1,nx*ny,1);
    emiDET1=reshape(emiDET,nx,ny)';
    emiDET2=reshape(emiDET1,nx*ny,1);
    emi=[];
    emi=[emiAGG2 emiDET2];
end

if coordinateFlip==2 %ASPA CASE
    flag_optim_dom1=reshape(flag_optim_dom,ny,nx);
    flag_optim_dom2=reshape(flipud(flag_optim_dom1),nx*ny,1);
    flag_optim_dom=flag_optim_dom2;
    
    %UNIBS(ET)20131001 - new regional flag
    flag_region_dom1=reshape(flag_region_dom,nx,ny)';
    flag_region_dom2=reshape(flag_region_dom1,nx*ny,1);
    flag_region_dom=flag_region_dom2;
    
    %new pop file
    pop1=reshape(pop,ny,nx);
    pop2=reshape(flipud(pop1),nx*ny,1);
    pop=pop2;
    
    %new aqi compuation domain file
    flag_aqi_dom1=reshape(flag_aqi_dom,ny,nx);
    flag_aqi_dom2=reshape(flipud(flag_aqi_dom1),nx*ny,1);
    flag_aqi_dom=flag_aqi_dom2;
    
    %new xutm
    xutm1=reshape(xutm,ny,nx);
    xutm2=reshape(flipud(xutm1),nx*ny,1);
    
    %new xutm
    yutm1=reshape(yutm,ny,nx);
    yutm2=reshape(flipud(yutm1),nx*ny,1);
    
    %new coordinates
    coordinate=[xutm2 yutm2];
    xutm=xutm2;
    yutm=yutm2;
    
    emiAGG=emi(:,1);
    emiDET=emi(:,2);
    emiAGG1=reshape(emiAGG,ny,nx);
    emiAGG2=reshape(flipud(emiAGG1),nx*ny,1);
    emiDET1=reshape(emiDET,ny,nx);
    emiDET2=reshape(flipud(emiDET1),nx*ny,1);
    emi=[];
    emi=[emiAGG2 emiDET2];
end

%define number of optimization cells: to be used to apply D and d in the
%routine MAINcompute_aqi (1 is for cells inside optim domain, 2 for border
%cells)
ncelopt=length(find(flag_optim_dom==1 | flag_optim_dom==2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%part in common with INIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load PM10 to PM25 relationship
pm10aveToExceed=load(pathPM);

%output file name
nomeOUTPUT=pathOUF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE LOG FILE WITH STATUS OF THE PROGRAMME
mkdir(pathOUF)
fileStatus=strcat(pathOUF,'/statusLOG.csv');
fidStatus=fopen(fileStatus,'wt');

fileExit=strcat(pathOUF,'/exitLOG.csv');
fidExit=fopen(fileExit,'wt');

%                    DATA STRUCTURES DEFINITION
%      (TO MANAGE ANN, D AND d FOR ALL THE POSSIBLE CONFIGURATIONS)

% DSuperSet(1) and nnSuperSet(1): annual
% DSuperSet(2) and nnSuperSet(2): winter
% DSuperSet(3) and nnSuperSet(3): summer
DSuperSet = struct('DSet', {});
nnSuperSet = struct('nnSet', {});
bcSuperSet = struct('bcSet', {});
for h=1:3,
    DSetHorizon = struct('D', {}, 'd', {}, 'Dp', {}, 'dp', {});
    % ET 20140331 - New classes for nnSetHorizon
    nnSetHorizon = struct('ps_input', {}, 'ps_target', {}, 'net', {}, 'icells', {},'Class', {}, 'PRECs', {},'ps_pca',{},'ArPt',{});
    % ET 20140331 - set bcSetHorizon - structure for basecases
    bcSetHorizon= struct('emi_bc',{}, 'conc_bc', {});
    
    for i=1:AQINum,
        if (isequal(strtrim(pathANN(h).ANNs(i,:)),'-999')==0)
            if flag_mode_ce_mo==3 %in case aggregated scenario analysis, do not load Dd
                D=struct('D', {-999}, 'd', {-999}, 'Dp', {-999}, 'dp', {-999});
            else
                D=load(strtrim(pathDd(h).Dd(i,:)));
            end
            %nn=load(strtrim((pathANN(h).ANNs(i,:))));
            [nn]=net_read(strtrim((pathANN(h).ANNs(i,:)))); %20140617 ET - read model from txt file
            % ET 20140331 load basecases for all nets
            if strcmp(nn.Class,'Delta')==1
                nnp=strtrim((pathANN(h).ANNs(i,:))); %net path
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
for o=1:optAQINum,
    DOptSet = [DOptSet, DSuperSet(aqi_horizon(o)).DSet(aqi_obj(o)+1)];
    nnOptSet = [nnOptSet, nnSuperSet(aqi_horizon(o)).nnSet(aqi_obj(o)+1)];
    bcOptSet = [bcOptSet, bcSuperSet(aqi_horizon(o)).bcSet(aqi_obj(o)+1)];
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
%UNIBS(ET)20131001 - added flag_region_dom in input
if (flag_mode_ce_mo==0 | flag_mode_ce_mo==1 | flag_mode_ce_mo==2)
    [actlev_final,actlev_final_sum,base_emi_low,base_emi_high,base_emi_low_noc,...
        base_emi_high_noc]=MAINbase_emi(emi,global_data,flag_optim_dom,flag_region_dom,global_data_noc);
    
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
        MAINcreate_constraints_matrix(global_data,CLE,flag_constraints,pathLST);
    
    
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint1_fix','A','B');
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
    if flag_constraints==1
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
    costPerMacrosectorMask = zeros(num_of_technologies,MNum);
    for m=1:MNum,
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
    % macrosector when technologies are applied at a given application
    % rate (see buildSolution function)
    
    costPerMacrosector = zeros(num_of_technologies,MNum);
    for m=1:MNum,
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
    costconstrPerMacrosectorLHS = zeros(num_of_technologies,MNum);
    for m=1:MNum,
        costconstrPerMacrosectorLHS(:,m) = ...
            costPerMacrosector(:,m) - costconstr * MPBudget(m);
    end;
    
    % rhs term
    costconstrPerMacrosectorRHS = zeros(MNum,1);
    for m=1:MNum,
        costconstrPerMacrosectorRHS(m,1) = ...
            minimumCostPerMacrosector(MId(m)) - ...
            MPBudget(m) * (sum(minimumCostPerMacrosector));
    end
    
    if(MNum > 0)
        A = [A; costconstrPerMacrosectorLHS'];
        B = [B; costconstrPerMacrosectorRHS];
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       TRAFFIC CONSTRAINTS
    
    % read configuration file
    flag_optim_oth = load(pathTRD);
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
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint2_fix','A1');
    
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
    Update();
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint3_fix','A','A1','costconstr');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           OPTIMIZATION
    % 'active-set', 'trust-region-reflective', 'interior-point', 'levenberg-marquardt','sqp'.
    % opt=optimset('LargeScale','off',...'Algorithm','sqp',...
    % NOTE: SQP is the best!!!
    opt=optimset('LargeScale','off','Algorithm','sqp',...
        'Display','Iter',...
        'Diagnostic','off',...
        'FunValCheck','off',...
        'TolFun',conv_value);
    
    % compute CLE values
    x_CLE_free=CLE;
    
    %ENR20130415
    strStatus='PROGRESSION: computing CLE results...';
    disp(strStatus);
    fprintf(fidStatus, '%s\n',strStatus);
    
    sol_CLE = buildSolution(x_CLE_free, fixed_x, inhibited_x, ...
        costconstr'*(x_CLE_free) + deltaCost, ...
        costconstr'*(x_CLE_free) + deltaCost);
    
    solutionList = [solutionList, sol_CLE];
    
    solutionSet = [solutionSet, ...
        struct('BUDGET', sol_CLE.BUDGET, 'SOL_LIST', solutionList)];
    
end
%20130424 END THE PART SPECIFIC TO MO, CE, DETAILED SA

% if multi-objective
%save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint14_fix', 'solutionSet');
if flag_mode_ce_mo==0       %multi-objective
    
    % fix the starting point for the optimization
    x0=CLE;
    
    %ENR20130415
    strStatus='PROGRESSION: computing MFR optimal results...';
    disp(strStatus);
    fprintf(fidStatus, '%s\n',strStatus);
    
    % compute best solution
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint15a_fix', 'A','B','LB','UB','x0');
    x_BEST_free = MAINcompute_sol(A,B,LB,UB,x0);
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint24_fix', 'x_BEST_free');
    
    sol_BEST = buildSolution(x_BEST_free, fixed_x, inhibited_x, ...
        costconstr'*(x_BEST_free) + deltaCost, ...
        costconstr'*(x_BEST_free) + deltaCost);
    
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint25_fix', 'sol_BEST');
    %     startingcputime=cputime;
    
    % if (flag_paretopoints > 0) you manage the case with a constraint on
    % the total computation time, and not the generation of a fixed number
    % of pareto curve points
    if(flag_paretopoints == 0)
        % 3 points are added to the pareto curve
        
        % always use CLE as starting point, at the moment
        x0=CLE;
        
        % pareto curve points definition
        for i = 1:3
            
            % different threshold (budget) if costs are in absolute or
            % percentage values
            if abs_perc==0
                thres_cost = sol_CLE.COST + tmp_thres_cost(i);
            else
                thres_cost=(sol_CLE.COST + ...
                    (sol_BEST.COST - sol_CLE.COST) * tmp_thres_cost(i));
            end
            
            % insert the cost constraint
            B1 = [B;thres_cost - deltaCost];
            
            %ENR20130415
            strStatus = sprintf('PROGRESSION: computing point %d optimal results...',i+1);
            disp(strStatus)
            fprintf(fidStatus, '%s\n',strStatus);
            
            % compute best solution
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint15b_fix', 'A1','B1','LB','UB','x0');
            x_sol_free = MAINcompute_sol(A1,B1,LB,UB,x0);
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint26a_fix', 'x_sol_free');
            
            sol = buildSolution(x_sol_free, fixed_x, inhibited_x, ...
                costconstr'*(x_sol_free) + deltaCost, ...
                thres_cost);
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint26b_fix', 'sol');
            
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
    if(flag_paretopoints > 0 )
        
        time = toc;
        % computation time to be dedicated to pareto curve
        maxTime=flag_paretopoints;
        
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
                
                %%save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint15c_fix', 'A1','B1','LB','UB','x0');
                x_sol_free = MAINcompute_sol(A1,B1,LB,UB,x0);
                %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint27a_fix', 'x_sol_free');
                
                sol = buildSolution(x_sol_free, fixed_x, inhibited_x, ...
                    costconstr'*(x_sol_free) + deltaCost, ...
                    thres_cost);
                %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint27b_fix', 'sol');

                
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
    
elseif flag_mode_ce_mo==1   %cost-effectiveness
    
    % fix the starting point for the optimization
    x0=CLE;
    
    % insert the cost constraint
    thres_cost = sol_CLE.COST + tmp_thres_cost(1);
    B1 = [B; thres_cost - deltaCost];
    
    strStatus='PROGRESSION: computing the optimal cost-effectiveness results...';
    disp(strStatus);
    fprintf(fidStatus, '%s\n',strStatus);
    
    %compute solution
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint15d_fix', 'A1','B1','LB','UB','x0');
    x_sol_free = MAINcompute_sol(A1,B1,LB,UB,x0);
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint28a_fix', 'x_sol_free');
    
    sol = buildSolution(x_sol_free, fixed_x, inhibited_x, ...
        costconstr'*(x_sol_free) + deltaCost, ...
        thres_cost);
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint28b_fix', 'sol');
    
    solutionList = struct('X', {}, 'X_free', {}, 'AQIs', {}, ...
        'COST', {}, 'BUDGET', {}, ...
        'COSTPERMACROSECTOR', {}, ...
        'BUDGETPERMACROSECTOR', {});
    
    solutionList = [solutionList, sol];
    
    solutionSet = [solutionSet, ...
        struct('BUDGET', thres_cost, 'SOL_LIST', solutionList)];
    
elseif flag_mode_ce_mo==2   % scenario mode
    
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
    
    sol = buildSolution(x_sol, fixed_x_scenariomode, [], ...
        originalcostconstr'*(x_sol), ...
        originalcostconstr'*(x_sol));
    
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint29a_fix', 'sol');
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
elseif flag_mode_ce_mo==3   %aggregated scenario analysis
    
    strStatus='PROGRESSION: computing the aggregated scenario analysis results...';
    disp(strStatus);
    fprintf(fidStatus, '%s\n',strStatus);
    
    %save('C:\data\work\projects\riat\RiatPlus-v6.1b\datasave\sa_p1_fix', 'flag_ADS_tp', 'xutm', 'yutm');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LOAD EMISSION REDUCTION PERCENTAGES
    %     aggsa=importdata(strcat(pathEMI,'aggregated_ads_perc_reductions.txt'));
    aggsa=importdata(pathADS);
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\sa1_fix', 'aggsa');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %COMPUTE REDUCED EMISSIONS (BOTH AREAL AND POINT...WE EXPECT BOTH TYPE
    %OF FILES)
    %conisder the seasons - yea, win, sum
    for sea=1:3
        for i=1:ncel
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
    mkdir(strcat(pathOUF,'\maps_aqi'));
    mkdir(strcat(pathOUF,'\maps_emi'));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %SAVING EMISSIONS, OF CLE AND OF REDUCED EMISSIONS
    for indpar=1:2 %2 points to be produced
        %UNIBS(ET)20130927
        %loop over seasons introduced only if data for that season exists
        for sea2=1:3
            if flag_ADS_tp(sea2)==2
                emiADS=emiCleRed{sea,indpar};
                for emiind=1:6
                    emiADSFinal=emiADS(:,emiind)+emiADS(:,emiind+6); %sum areal and point
                    emiADSFinal(flag_optim_dom==0)=-999;
                    file=strcat(pathOUF,'/maps_emi/',nomeemi{emiind},nomesea{sea2},'_point',int2str(indpar),'.csv');
                    fid=fopen(file,'wt');
                    fprintf(fid,'%s %c %s %c %s \n','xutm[km]', delim, 'yutm[km]', delim, strcat(nomeemi{emiind},emiUM{emiind}));
                    dlmwrite(file,[xutm yutm emiADSFinal],'-append','roffset', 0,'precision','%6.2f');
                    fclose(fid);
                end
            end
        end
    end
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\sa3_fix', 'emiADSFinal');
    
    %LOOP OVER AQIS - CONSIDER ONLY YEARLY EMISSIONS AND AQIS
    for ii = 1:3, %season
        for jj = 1:AQINum, %aqi
            %UNIBS(ET)20130927
            %Compute AQI only if data for that season exists
            if flag_ADS_tp(ii)==2
                if (isequal(strtrim(pathANN(ii).ANNs(jj,:)),'-999')==0)
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %COMPUTE AQI PER CELL, RELATED TO INITIAL AND REDUCED CASES
                    %CONSTRAINTS ON ANNS IDENTIFICATION BOUNDS ARE CHECKED HERE
                    [aqiCLEini]=MAINaggregated_scenario_mode(emiCleRed{ii,1},ii,jj);
                    [aqiCLEred]=MAINaggregated_scenario_mode(emiCleRed{ii,2},ii,jj);
                    
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
                        
                        file=strcat(pathOUF,'/maps_aqi/',nomeaqi{jj},nomesea{ii},'_point',int2str(indpar),'.csv');
                        fid=fopen(file,'wt');
                        fprintf(fid,'%s %c %s %c %s \n','xutm[km]', delim, 'yutm[km]', delim, strcat(nomeaqi{jj},aqiUM{jj}));
                        dlmwrite(file,[xutm yutm aqiFinal],'-append','roffset', 0,'precision','%6.2f');
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

if (flag_mode_ce_mo==0 | flag_mode_ce_mo==1 | flag_mode_ce_mo==2)
    CLE = global_data(:,19)/100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save results in mat file
 save(nomeOUTPUT,'paretoSolutions','DSuperSet','nnSuperSet','solutionSet');
toc

strStatus='PROGRESSION: Computation part is finished.';
disp(strStatus);
fprintf(fidStatus, '%s\n',strStatus);

% fprintf(fidStatus, '%s\n','Computation part is finished.');

if flag_mode_ce_mo<3
    % %postprocessing of results
    %UNIBS(ET)20131001 - added flag_region_dom in input
    %UNIBS(ET)20140403 - added bcSuperSet
    POSTmain(f2,nnSuperSet,DSuperSet,bcSuperSet,pathANN, flag_region_dom, flag_optim_dom, flag_aqi_dom,...
        coordinate,emi,global_data,base_emi_low,base_emi_high,CLE, ...
        base_emi_low_noc,base_emi_high_noc,pathOCM,paretoSolutions,...
        areal_point,flag_reg_net,pathOUF,pop,nocID,pathAR,pm10aveToExceed,...
        cell_threshold_set,fidStatus,fidExit);
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
        global_data_tc=load(pathOCM);
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
            defaultTrafficConstraintsSubmatrixA(find(~(global_data_tc(:,4)==nocID)),:);
        
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
            defaultTrafficConstraintsSubmatrixB(find(~(global_data_tc(:,4)==nocID)),:);
        
        trafficConstraintsSubmatrixB = trafficConstraintsSubmatrixB';
        
        
        % define the final matrix
        trafficConstraintsMatrix = ...
            [trafficConstraintsSubmatrixA;...
            trafficConstraintsSubmatrixB];
        
        
        clear defaultTrafficConstraintsMatrix;
        clear global_data_tc;
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Update()
        
        %if techs change emissions less than threshold (sum for all
        %pollutants) the techs are disregarded
        threshold = 1e-6;
        
        final_active_technologies = zeros(num_of_technologies,1);
        
        for ii = 1:optAQINum,
            
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
            if (areal_point==1)
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
                
                if (((flag_constraints == 0) && (flag_usab_tech(tt,1) == 0)) || ...
                        ((flag_constraints == 0) && (flag_usab_tech(tt,1) == 1) && (totalReductionFactor(tt,1) <= threshold)) || ...
                        ((flag_constraints == 0) && (flag_usab_tech(tt,1) == 1) && (abs(UB(tt)-LB(tt))<=1e-4)) || ...
                        ((flag_constraints == 1) && (flag_usab_tech(tt,1) == 0) && (flag_repl_tech(tt,1) == 1)) || ...
                        ((flag_constraints == 1) && (flag_usab_tech(tt,1) == 1) && (flag_repl_tech(tt,1) == 1) && (totalReductionFactor(tt,1) <= threshold)) || ...
                        ((flag_constraints == 1) && (flag_usab_tech(tt,1) == 1) && (flag_repl_tech(tt,1) == 1) && (abs(UB(tt)-LB(tt))<=1e-4)))
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
        
        for ii = 1:optAQINum,
            %in case of areal+point summed, D and d contain infos both or areal
            %and point, together
            %in case of areal and point separated, D and d are only related to
            %areal emissions (the point infos Dp and dp will be used later on)
            DOptSet(ii).d = DOptSet(ii).d - (DOptSet(ii).D * inhibited_x);
            if (areal_point==1)
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
        for ii = 1:optAQINum,
            DOptSet(ii).D = DOptSet(ii).D(:,(final_active_technologies==1));
            if (areal_point==1)
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
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\clue_fix', 'DOptSet');

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function sol = buildSolution(x_free, fixed_x, inhibited_x, ...
            cost, budget)
        
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
        
        cost_Per_Macrosector = zeros(MNum,1);
        for mm=1:MNum,
            cost_Per_Macrosector(mm,1) = ...
                costPerMacrosector(:,mm)'*x;
        end;
        
        %computing the AQIS for the 3 temporal periods, 6 aqis, 3
        %aggregation types
        aqi = ones(3,AQINum,3) * -999;
        for ii = 1:3,
            for jj = 1:AQINum,
                if (isequal(strtrim(pathANN(ii).ANNs(jj,:)),'-999')==0)
                    %20140403 ET - Basecase inputs added
                    aqi(ii,jj,1) = MAINcompute_aqi(x,nnSuperSet(ii).nnSet(jj),DSuperSet(ii).DSet(jj),0,cell_threshold_set(jj),jj-1,1,bcSuperSet(ii).bcSet(jj));
                    aqi(ii,jj,2) = MAINcompute_aqi(x,nnSuperSet(ii).nnSet(jj),DSuperSet(ii).DSet(jj),1,cell_threshold_set(jj),jj-1,1,bcSuperSet(ii).bcSet(jj));
                    aqi(ii,jj,3) = MAINcompute_aqi(x,nnSuperSet(ii).nnSet(jj),DSuperSet(ii).DSet(jj),2,cell_threshold_set(jj),jj-1,1,bcSuperSet(ii).bcSet(jj));
                    %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint12_fix', 'aqi');
                end
                
            end
        end
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint13_fix', 'aqi');
        
        sol = struct('X', x, 'X_free', x_free, 'AQIs', aqi, ...
            'COST', cost, 'BUDGET', budget, ...
            'COSTPERMACROSECTOR', cost_Per_Macrosector, ...
            'BUDGETPERMACROSECTOR', budget * MPBudget);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%30140403 ET - Added BCset input
    function aqi_per_cell=MAINcompute_aqi_per_cell(x,NN,DD,BCset)
        
        
        % in case of areal+point summed, D and d contain infos both or areal
        % and point, together
        % in case of areal and point separated, D and d are only related to
        % areal emissions (the point infos Dp and dp will be used later on)
        %this will be used when AL will be a decision variable
        %       disp(strcat('PROGRESSION: Creating preprocessed emissions, for AQI -->',strtrim(pathANN(k).ANNs(indaqi,:))))
        %                 [D,d,Dp,dp,CellPerc]=INITquadrant(emi,global_data,...
        %                     base_emi_low, base_emi_high, base_emi_low_noc,...
        %                     base_emi_high_noc,flag_region_dom,nx,ny,NN,areal_point);
        
        
        D=DD.D;
        d=DD.d;
        
        % E is the reduced emissions, depending on application rates
        % e are the initial emissions
        
        % calculation starts from sparse matrices
        E = D * sparse(x);
        E = d - E;
        E_full = full(E);
        
        %%save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint5_fix','E_full');
        
        
        %from global matrix, extract quadrant informations
        s1_NOX = E_full((ncelopt*0)+1:ncelopt*1);
        s2_NOX = E_full((ncelopt*1)+1:ncelopt*2);
        s3_NOX = E_full((ncelopt*2)+1:ncelopt*3);
        s4_NOX = E_full((ncelopt*3)+1:ncelopt*4);
        s1_VOC = E_full((ncelopt*4)+1:ncelopt*5);
        s2_VOC = E_full((ncelopt*5)+1:ncelopt*6);
        s3_VOC = E_full((ncelopt*6)+1:ncelopt*7);
        s4_VOC = E_full((ncelopt*7)+1:ncelopt*8);
        s1_NH3 = E_full((ncelopt*8)+1:ncelopt*9);
        s2_NH3 = E_full((ncelopt*9)+1:ncelopt*10);
        s3_NH3 = E_full((ncelopt*10)+1:ncelopt*11);
        s4_NH3 = E_full((ncelopt*11)+1:ncelopt*12);
        s1_PM10 = E_full((ncelopt*12)+1:ncelopt*13);
        s2_PM10 = E_full((ncelopt*13)+1:ncelopt*14);
        s3_PM10 = E_full((ncelopt*14)+1:ncelopt*15);
        s4_PM10 = E_full((ncelopt*15)+1:ncelopt*16);
        s1_PM25 = E_full((ncelopt*16)+1:ncelopt*17);
        s2_PM25 = E_full((ncelopt*17)+1:ncelopt*18);
        s3_PM25 = E_full((ncelopt*18)+1:ncelopt*19);
        s4_PM25 = E_full((ncelopt*19)+1:ncelopt*20);
        s1_SO2 = E_full((ncelopt*20)+1:ncelopt*21);
        s2_SO2 = E_full((ncelopt*21)+1:ncelopt*22);
        s3_SO2 = E_full((ncelopt*22)+1:ncelopt*23);
        s4_SO2 = E_full((ncelopt*23)+1:ncelopt*24);
        
        %create input vector with rigth order of input
        NH3_all=[s1_NH3,s2_NH3,s3_NH3,s4_NH3];
        NOX_all=[s1_NOX,s2_NOX,s3_NOX,s4_NOX];
        PM10_all=[s1_PM10,s2_PM10,s3_PM10,s4_PM10];
        PM25_all=[s1_PM25,s2_PM25,s3_PM25,s4_PM25];
        SO2_all=[s1_SO2,s2_SO2,s3_SO2,s4_SO2];
        VOC_all=[s1_VOC,s2_VOC,s3_VOC,s4_VOC];
        
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\emis_fix','NH3_all', ...
                   % 'NOX_all', 'PM10_all', 'PM25_all', 'SO2_all', 'VOC_all');
        %CASE OF AREAL+POINT, SUMMED
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint6_fix','NH3_all', ...
        %            'NOX_all', 'PM10_all', 'PM25_all', 'SO2_all', 'VOC_all');
        %disp('******');
        %disp('NN.PRECs');
        %disp(NN.PRECs);
        %disp('******');
        if areal_point==0
            
            %ET20140313
            %input selection based on ANNs features
            
            %read net and select precursors
            emissioni=[];
            for i_prec=1:size(NN.PRECs,2) %areal
                if strcmp(NN.PRECs(1,i_prec),'NH3')==1
                    emissioni=[emissioni, NH3_all];
                elseif strcmp(NN.PRECs(1,i_prec),'NOX')==1
                    emissioni=[emissioni, NOX_all];
                elseif strcmp(NN.PRECs(1,i_prec),'PM10')==1
                    emissioni=[emissioni, PM10_all];
                elseif strcmp(NN.PRECs(1,i_prec),'PM25')==1
                    emissioni=[emissioni, PM25_all];
                elseif strcmp(NN.PRECs(1,i_prec),'SO2')==1
                    emissioni=[emissioni, SO2_all];
                elseif strcmp(NN.PRECs(1,i_prec),'VOC')==1
                    emissioni=[emissioni, VOC_all];
                end
            end
            
            %read net and create Delta emissions
            %emi delta
            if strcmp(NN.Class,'Delta')==1
                %20140403 ET - Filtering basecase cells
                BCemi=BCset.emi_bc(flag_optim_dom==1 | flag_optim_dom==2,:);
                BCconc=BCset.conc_bc(flag_optim_dom==1 | flag_optim_dom==2,:);
                
                emissioni2=(emissioni-BCemi)./BCemi;
                emissioni=emissioni2;
                test_nan=find(isnan(emissioni));
                emissioni(test_nan)=0.;
            end
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPointEmissioni_a_fix','emissioni');

            % OLD PROCEDURE
            %ENR 20130314 - managing 03 ANNs with 2 input
            %             %for ozone ANN only NOX and VOC as input
            %             %for ozone ANN only NOX and VOC as input
            %             %consider linear and non linear case
            %             %IF IT IS A ANN: size(NN.net,1)==1 (SO CHECK NN.net.inputs)
            %             %IF IT IS A LIN: check dimension of size(NN.net,1), that is the
            %             %matrix of weights
            %             if (size(NN.net,1)==24 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==24))% ANN/linear, 6 input
            %                 emissioni=[NH3_all,NOX_all,PM10_all,PM25_all,SO2_all,VOC_all];
            %
            %             elseif (size(NN.net,1)==8 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==8))%ANNlinear, 2 input
            %                 emissioni=[NOX_all,VOC_all];
            %
            %             end
            %             %             emissioni=[NH3_all,NOX_all,PM10_all,PM25_all,SO2_all,VOC_all];
            %
            %CASE FOR POINT EMISSIONS (AREAL INFOS HAVE BEEN PROCESSED IN
            %THE LINES ABOVE)
        elseif areal_point==1
            Dp=DD.Dp;
            dp=DD.dp;
            
            Ep = Dp * sparse(x);
            Ep = dp - Ep;
            Ep_full = full(Ep);
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\var_buildEmission_new','Ep_full');
            
            s1_NOXp = Ep_full((ncelopt*0)+1:ncelopt*1);
            s2_NOXp = Ep_full((ncelopt*1)+1:ncelopt*2);
            s3_NOXp = Ep_full((ncelopt*2)+1:ncelopt*3);
            s4_NOXp = Ep_full((ncelopt*3)+1:ncelopt*4);
            s1_VOCp = Ep_full((ncelopt*4)+1:ncelopt*5);
            s2_VOCp = Ep_full((ncelopt*5)+1:ncelopt*6);
            s3_VOCp = Ep_full((ncelopt*6)+1:ncelopt*7);
            s4_VOCp = Ep_full((ncelopt*7)+1:ncelopt*8);
            s1_NH3p = Ep_full((ncelopt*8)+1:ncelopt*9);
            s2_NH3p = Ep_full((ncelopt*9)+1:ncelopt*10);
            s3_NH3p = Ep_full((ncelopt*10)+1:ncelopt*11);
            s4_NH3p = Ep_full((ncelopt*11)+1:ncelopt*12);
            s1_PM10p = Ep_full((ncelopt*12)+1:ncelopt*13);
            s2_PM10p = Ep_full((ncelopt*13)+1:ncelopt*14);
            s3_PM10p = Ep_full((ncelopt*14)+1:ncelopt*15);
            s4_PM10p = Ep_full((ncelopt*15)+1:ncelopt*16);
            s1_PM25p = Ep_full((ncelopt*16)+1:ncelopt*17);
            s2_PM25p = Ep_full((ncelopt*17)+1:ncelopt*18);
            s3_PM25p = Ep_full((ncelopt*18)+1:ncelopt*19);
            s4_PM25p = Ep_full((ncelopt*19)+1:ncelopt*20);
            s1_SO2p = Ep_full((ncelopt*20)+1:ncelopt*21);
            s2_SO2p = Ep_full((ncelopt*21)+1:ncelopt*22);
            s3_SO2p = Ep_full((ncelopt*22)+1:ncelopt*23);
            s4_SO2p = Ep_full((ncelopt*23)+1:ncelopt*24);
            
            NH3_allp=[s1_NH3p,s2_NH3p,s3_NH3p,s4_NH3p];
            NOX_allp=[s1_NOXp,s2_NOXp,s3_NOXp,s4_NOXp];
            PM10_allp=[s1_PM10p,s2_PM10p,s3_PM10p,s4_PM10p];
            PM25_allp=[s1_PM25p,s2_PM25p,s3_PM25p,s4_PM25p];
            SO2_allp=[s1_SO2p,s2_SO2p,s3_SO2p,s4_SO2p];
            VOC_allp=[s1_VOCp,s2_VOCp,s3_VOCp,s4_VOCp];

            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\emis_p_fix','NH3_allp', ...
            %        'NOX_allp', 'PM10_allp', 'PM25_allp', 'SO2_allp', 'VOC_allp');
            %ENR20130313
            %for ozone ANN only NOX and VOC as input
            %consider linear and non linear case
            
            %ET20140313
            %input selection based on ANNs features
            
            %read net and select precursors
            emissioni=[];
            for i_prec=1:size(NN.PRECs,2) %areal
                if strcmp(NN.PRECs(1,i_prec),'NH3')==1
                    emissioni=[emissioni, NH3_all];
                elseif strcmp(NN.PRECs(1,i_prec),'NOX')==1
                    emissioni=[emissioni, NOX_all];
                elseif strcmp(NN.PRECs(1,i_prec),'PM10')==1
                    emissioni=[emissioni, PM10_all];
                elseif strcmp(NN.PRECs(1,i_prec),'PM25')==1
                    emissioni=[emissioni, PM25_all];
                elseif strcmp(NN.PRECs(1,i_prec),'SO2')==1
                    emissioni=[emissioni, SO2_all];
                elseif strcmp(NN.PRECs(1,i_prec),'VOC')==1
                    emissioni=[emissioni, VOC_all];
                end
            end
            for i_prec=1:size(NN.PRECs,2) %point
                
                if strcmp(NN.PRECs(1,i_prec),'NH3')==1
                    emissioni=[emissioni, NH3_allp];
                elseif strcmp(NN.PRECs(1,i_prec),'NOX')==1
                    emissioni=[emissioni, NOX_allp];
                elseif strcmp(NN.PRECs(1,i_prec),'PM10')==1
                    emissioni=[emissioni, PM10_allp];
                elseif strcmp(NN.PRECs(1,i_prec),'PM25')==1
                    emissioni=[emissioni, PM25_allp];
                elseif strcmp(NN.PRECs(1,i_prec),'SO2')==1
                    emissioni=[emissioni, SO2_allp];
                elseif strcmp(NN.PRECs(1,i_prec),'VOC')==1
                    emissioni=[emissioni, VOC_allp];
                end
            end
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPointEmissioni_step1_fix','emissioni');
            %read net and sum Areal and Point
            if strcmp(NN.ArPt,'Aggregated')==1
                emissioni=emissioni(:,1:(size(emissioni,2)/2))+emissioni(:,(size(emissioni,2)/2)+1:size(emissioni,2));
            end
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPointEmissioni_step2_fix','emissioni');
            
            %read net and create Delta emissions
            %emi delta
            if strcmp(NN.Class,'Delta')==1
                %20140403 ET - Filtering basecase cells
                BCemi=BCset.emi_bc(flag_optim_dom==1 | flag_optim_dom==2,:);
                BCconc=BCset.conc_bc(flag_optim_dom==1 | flag_optim_dom==2,:);
                
                emissioni2=(emissioni-BCemi)./BCemi;
                emissioni=emissioni2;
                test_nan=find(isnan(emissioni));
                emissioni(test_nan)=0.;
            end
            
            % OLD PROCEDURE
            %             if (size(NN.net,1)==48 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==48))% ANN/linear, 6 input
            %                 emissioni=[NH3_all,NOX_all,PM10_all,PM25_all,SO2_all,VOC_all,...
            %                     NH3_allp,NOX_allp,PM10_allp,PM25_allp,SO2_allp,VOC_allp];
            %
            %             elseif (size(NN.net,1)==16 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==16))%ANNlinear, 2 input
            %                 emissioni=[NOX_all,VOC_all,NOX_allp,VOC_allp];
            %
            %             end
            %
            %                         elseif size(NN.net,1)==1 %ANNs
            %                             if NN.net.inputs{1}.size==48 %ANNs, 6 input
            %                                 emissioni=[NH3_all,NOX_all,PM10_all,PM25_all,SO2_all,VOC_all,...
            %                                     NH3_allp,NOX_allp,PM10_allp,PM25_allp,SO2_allp,VOC_allp];
            %
            %                             elseif NN.net.inputs{1}.size==16 %ANNs, 2 input
            %                                 emissioni=[NOX_all,VOC_all,NOX_allp,VOC_allp];
            %                             end
            %                         end
            %
            %                         emissioni=[NH3_all,NOX_all,PM10_all,PM25_all,SO2_all,VOC_all,...
            %                             NH3_allp,NOX_allp,PM10_allp,PM25_allp,SO2_allp,VOC_allp];
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPointEmissioni_step3_fix','emissioni');

        end
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\emissioni_fix','emissioni');
        %load network
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\emissioni_fix','emissioni');
        input_rete2=emissioni';
        
        %in case it is necessary to process quadrant emissions (if too close
        %to domain boundary, it is necessary to increment emissions with
        %the assumptions that part of the quadrant in which emissions are
        %not available, still contain same emission average)
        if strcmp(pathAR,'-1')==0
            %load path for Area_ratio variable...file contain ratios  with
            %the order as NORTH SOUTH NWEST EAST
            
            %remove cells not in flag_optim_dom
            
            %repmat to get same "emissioni'"
            %dimensions
            AR=load(pathAR);                             %load not in optim
            AR=AR.Ratio;                                 %rename variable
            AR(flag_optim_dom==0,:)=[];                  %remove cells outside flag_optim_dom
            
            %ENR20130313
            %manage case with 2 input for ozone, or 6 input for ozone and
            %PM ANNs
            %             if NN.net.inputs{1}.size==48
            %                 ARreordrepm=repmat(AR,1,12);            %repmat
            %             elseif NN.net.inputs{1}.size==16
            %                 ARreordrepm=repmat(AR,1,4);            %repmat
            %             end
            
            
            %             if (size(NN.net.IW{1,1},2)==48 || (size(NN.net.IW{1,1},2)==1 && NN.net.inputs{1}.size==48))% ANN/linear, 6 input
            %                 ARreordrepm=repmat(AR,1,12);            %repmat
            %
            %             elseif (size(NN.net.IW{1,1},2)==16 || (size(NN.net.IW{1,1},2)==1 && NN.net.inputs{1}.size==16))%ANNlinear, 2 input
            %                 ARreordrepm=repmat(AR,1,4);            %repmat
            %
            %             end
            
            if (size(NN.net.IW{1,1},2)==48)% ANN/linear, 6 input
                ARreordrepm=repmat(AR,1,12);            %repmat
                
            elseif (size(NN.net.IW{1,1},2)==16)%ANNlinear, 2 input
                ARreordrepm=repmat(AR,1,4);            %repmat
                
            end
            
            
            %             if size(NN.net,1)==48 % linear, 6 input
            %                 ARreordrepm=repmat(AR,1,12);            %repmat
            %
            %             elseif size(NN.net,1)==16 %linear, 2 input
            %                 ARreordrepm=repmat(AR,1,4);            %repmat
            %
            %             elseif size(NN.net,1)==1 %ANNs
            % %                 if NN.net.inputs{1}.size==48 %ANNs, 6 input
            %                     ARreordrepm=repmat(AR,1,12);            %repmat
            %
            %                 elseif NN.net.inputs{1}.size==16 %ANNs, 2 input
            %                     ARreordrepm=repmat(AR,1,4);            %repmat
            %                 end
            %             end
            
            
            %             ARreord=[AR(:,2) AR(:,3) AR(:,1) AR(:,4)];   %reorder
            %             ARreordrepm=repmat(AR,1,12);            %repmat
            input_rete2=emissioni'./ARreordrepm';        %rewrite input_rete2 dividing emissions by area_ratio
        end
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\MAINcompute_aqi_per_cell_2_fix','input_rete2');
        
        %define if lin o net
        if size(NN.net,1)==1
            flag_reg_net=1;
        else
            flag_reg_net=0;
        end
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint8_fix','input_rete2');
        
        switch flag_reg_net
            case 0 %linear case
                %20140404 and 20140612(ET) - Adapted to the new linear models
                icel=0;
                input_rete=input_rete2';
                for dimen1=1:size(input_rete,2)
                    for dimen2=dimen1:size(input_rete,2)
                        icel=icel+1;
                        TMPnet(:,icel)=input_rete(:,dimen1).*input_rete(:,dimen2);
                    end
                end
                aqi_per_cell=NN.net'*[input_rete TMPnet]';
                %if delta net transform net output in abs values
                if strcmp(NN.Class,'Delta')==1
                    aqi_per_cell=BCconc'+(aqi_per_cell.*BCconc');
                end
                myNet=NN.net;
                %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\MAINcompute_aqi_per_cell_3_fix','myNet','aqi_per_cell');
                %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\MAINcompute_aqi_per_cell_4_fix','aqi_per_cell');
                
            case 1 %nonlinear case
                %20111215 - run the network only on the optimization cells
                
                %20141205-ET - mapminmax replaced with a handwritten normalization
                %to deal with the problems of standalone windows executables
                ymin=NN.ps_input.ymin;
                ymax=NN.ps_input.ymax;
                xmin=NN.ps_input.xmin;
                xmax=NN.ps_input.xmax;
                G1=(ymax-ymin)./(xmax-xmin);
                
                input_rete_norm2=((input_rete2-repmat(xmin,1,size(input_rete2,2))).*repmat(G1,1,size(input_rete2,2))+repmat(ymin,size(input_rete2,1),size(input_rete2,2)));
                myNetFix=NN;
                %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\MAINcompute_aqipercell_3_fix','myNetFix', 'input_rete_norm2');
                %NN.ps_input.no_change=0; %to be compatible among different matlab versions
                %[input_rete_norm2]=mapminmax('apply',input_rete2,NN.ps_input);
                
                %%%
                %check for results outside bounds - force to be inside bounds
                input_rete_norm2(find(input_rete_norm2<-1))=-1;
                input_rete_norm2(find(input_rete_norm2>1))=1;
                %%%
                
                %20141121-ET - Use new sim_exe function instead of matlab
                %sim function
                output_rete_norm2=sim_exe(NN,input_rete_norm2);
                %%%
                
                %20141205-ET - mapminmax replaced with a handwritten normalization
                %to deal with the problems of standalone windows
                %executables, matlab transformation included
                G=NN.ps_target.gain;
                offst=NN.ps_target.offset;
                ymin=NN.ps_target.ymin;
                ymax=NN.ps_target.ymax;
                xmin=NN.ps_target.xmin;
                xmax=NN.ps_target.xmax;
                
                aqi_per_cell=((output_rete_norm2-ymin)*G*(xmax-xmin))/(ymax-ymin)+xmin+((xmax-xmin)*(-ymin+offst)/(ymax-ymin));
                %                 NN.ps_target.no_change=0; %to be compatible among different matlab versions
                %                 aqi_per_cell=mapminmax('reverse',output_rete_norm2,NN.ps_target);
                
                %check for results outside bounds - force to be inside bounds
                aqi_per_cell(find(aqi_per_cell<xmin))=xmin;
                aqi_per_cell(find(aqi_per_cell>xmax))=xmax;
                %%%
                
                %30130402  ET -  if delta net transform net output in abs values
                if strcmp(NN.Class,'Delta')==1
                    aqi_per_cell=BCconc'+(aqi_per_cell.*BCconc');
                end
                %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\MAINcompute_aqi_per_cell_nl_fix','aqi_per_cell','output_rete_norm2');
        end
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint9_fix','aqi_per_cell');
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%30140403 ET - Added BCset input
    function aqi_val = ...
            MAINcompute_aqi(x,nn,D,obj_function_type,cell_threshold,AQIId,REAL,BCset)
        
        % compute aqi per cell
        %30140403 ET - Added BCset input
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint4_fix',...
        %'x', 'D', 'obj_function_type', 'cell_threshold', 'AQIId', 'REAL','BCset');
        
        aqi_per_cell = MAINcompute_aqi_per_cell(x,nn,D,BCset);
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint10_fix', 'aqi_per_cell');
        
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
                indpop=pop;
                indpop(flag_aqi_dom==0)=[];
                aqi_val=sum(aqi_per_cell'.*indpop)/sum(indpop);
        end
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint11_fix', 'aqi_val');
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function aqi_val=MAINcompute_aqi_mp(x,optAQIBestValues)
        
        optAQIValues = zeros(optAQINum,1);
        
        for ii = 1:optAQINum,
            % 20140403 ET - Added bcSuperSet(ii) input
            myD=DOptSet(ii);
            mybc=bcOptSet(ii);
            mynn=nnOptSet(ii);
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint21_fix', 'myD','mybc', 'mynn');

            optAQIValues(ii,1) = MAINcompute_aqi(x,nnOptSet(ii),DOptSet(ii),...
                aqi_obj_function_type(ii),...
                cell_threshold_set(aqi_obj(ii)+1),...
                aqi_obj(ii),...
                0,bcOptSet(ii));
            
            optAQIcle(ii,1)=sol_CLE.AQIs(aqi_horizon(ii),aqi_obj(ii)+1,aqi_obj_function_type(ii)+1);
            
        end
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint22_fix', 'optAQIcle', 'optAQIBestValues');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %NEW VERSION
        %20130430 - testing the claudio proposal to improve objective
        %function computation
        delta_init=optAQIcle-optAQIBestValues;
        delta_prod=ones(optAQINum,1);
        for ii=1:optAQINum
            for kk=1:optAQINum
                if kk~=ii
                    delta_prod(ii,1)=delta_prod(ii,1).*delta_init(kk,1);
                end
            end
        end

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
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint_23_fix', 'aqi_weights', 'delta_prod', 'aqi_weights_init', 'aqi_val', 'delta_prod');
        
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

    function x_sol = MAINcompute_sol(A,B,LB,UB,x0)
        % 20140403 ET - Added bcSuperSet(ii) input to MAINcompute_aqi calls
        if(optAQINum == 1)
            myD=DOptSet(1);
            mybcOptSet=bcOptSet(1);
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint18a_fix', 'myD','mybcOptSet');
            
            x_sol=fmincon(@(x) MAINcompute_aqi(x,nnOptSet(1),DOptSet(1),...
                aqi_obj_function_type(1),...
                cell_threshold_set(aqi_obj(1)+1),...
                aqi_obj(1),...
                0,bcOptSet(1)),...
                x0, A, B, [], [], LB, UB, [], opt);
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\var_x_sol_fix', 'x_sol');
        else
            
            optAQIBestValues = zeros(optAQINum,1);
            
            for ii = 1:optAQINum,
                
                myD=DOptSet(ii);
                mybcOptSet=bcOptSet(ii);
                %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint18b_fix', 'myD','mybcOptSet');
                x_aqi=fmincon(@(x) MAINcompute_aqi(x,nnOptSet(ii),DOptSet(ii),...
                    aqi_obj_function_type(ii),...
                    cell_threshold_set(aqi_obj(ii)+1),...
                    aqi_obj(ii),...
                    0,bcOptSet(ii)),...
                    x0, A, B, [], [], LB, UB, [], opt);
                
                %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\var_x_aqi_fix', 'x_aqi');
                optAQIBestValues(ii,1) = MAINcompute_aqi(x_aqi,nnOptSet(ii),DOptSet(ii),...
                    aqi_obj_function_type(ii),...
                    cell_threshold_set(aqi_obj(ii)+1),...
                    aqi_obj(ii),....
                    0,bcOptSet(ii));
            end
            
            opt=optimset(opt,'TolFun',1e-6,'MaxIter',70);
            
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\var_optAQIBestValues_fix','optAQIBestValues');
            %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint19_fix','optAQIBestValues',...
            %    'x0', 'A', 'B', 'LB', 'UB', 'opt');
            x_sol=fmincon(@(x) MAINcompute_aqi_mp(x,optAQIBestValues),...
                x0, A, B, [], [], LB, UB, [], opt);
            
            opt=optimset(opt,'TolFun',conv_value,'MaxIter',400);
            
        end
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varpoint23_fix', 'x_sol','opt');
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [spi1,spi2,spi3,spi4]=MAINquadrant_emi(valori_emi,r,nx,ny)
        
        %compute quadrant emissions
        
        %domain dimensions and other variables
        dimx=nx;
        dimy=ny;
        prova_emi=valori_emi(:,1);
        emi=reshape(prova_emi,dimy,dimx);
        
        %enlarged emission, with increased zeros all around (used to be able to
        %perform simple products among emi and factor matrix
        emi_enlarged=zeros(size(emi,1)+2*(r-2),size(emi,2)+2*(r-2));
        emi_enlarged(r-2+1:r-2+1+size(emi,1)-1,r-2+1:r-2+1+size(emi,2)-1)=emi;
        [sa sb]=size(emi_enlarged);
        
        %
        %left quadrant
        %factors to perform products of cell emissions
        fattore=tril(ones(r+1,r+1))-diag(repmat(0.5,r+1,1));
        sotto=flipud(fattore);
        fattore=[fattore;sotto(2:end,:)];
        fattore(r+1,r+1)=0.25;
        
        %initialize variable, and set dimensions
        spicchio3=zeros(dimy,dimx);
        
        %down quadrant
        %factors to perform products of cell emissions
        fattoreD=rot90(fattore,3);
        
        %initialize variable, and set dimensions
        spicchio1=zeros(dimy,dimx);
        
        %right quadrant
        %factors to perform products of cell emissions
        fattoreR=rot90(fattore,2);
        
        %initialize variable, and set dimensions
        spicchio4=zeros(dimy,dimx);
        
        %up quadrant
        %factors to perform products of cell emissions
        fattoreU=rot90(fattore,1);
        
        %initialize variable, and set dimensions
        spicchio2=zeros(dimy,dimx);
        
        %loop to aggregate emissions
        for i=r+1:sa-r
            for j=r+1:sb-r
                %left quadrant
                spicchio3(i-r+2,j-r+2)=sum(sum(emi_enlarged(i-r:i+r,j-r:j).*fattore));
                %up quadrant (down as matrix, up geographically speaking)
                spicchio1(i-r+2,j-r+2)=sum(sum(emi_enlarged(i-r:i,j-r:j+r).*fattoreD));
                %right quadrant
                spicchio4(i-r+2,j-r+2)=sum(sum(emi_enlarged(i-r:i+r,j:j+r).*fattoreR));
                %down quadrant (up as matrix, down geographically speaking)
                spicchio2(i-r+2,j-r+2)=sum(sum(emi_enlarged(i:i+r,j-r:j+r).*fattoreU));
            end
        end
        
        spi1=reshape(spicchio1,dimx*dimy,1);
        spi2=reshape(spicchio2,dimx*dimy,1);
        spi3=reshape(spicchio3,dimx*dimy,1);
        spi4=reshape(spicchio4,dimx*dimy,1);
        
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTING AGGREGATED SCENARIO MODE SOLUTIONS
    function [aqi_per_cell]=MAINaggregated_scenario_mode(emiTMP,ii,jj)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %ANNS CONSIDERED
        NN=nnSuperSet(ii).nnSet(jj);
        %save('C:\data\work\projects\riat\RiatPlus-v6.1b\datasave\sa_p2_fix', 'emiTMP','NN');
        icells=NN.icells;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %AREAL CASE
        [s1_NOX,s2_NOX,s3_NOX,s4_NOX]=MAINquadrant_emi(emiTMP(:,1),icells,nx,ny);
        [s1_VOC,s2_VOC,s3_VOC,s4_VOC]=MAINquadrant_emi(emiTMP(:,2),icells,nx,ny);
        [s1_NH3,s2_NH3,s3_NH3,s4_NH3]=MAINquadrant_emi(emiTMP(:,3),icells,nx,ny);
        [s1_PM10,s2_PM10,s3_PM10,s4_PM10]=MAINquadrant_emi(emiTMP(:,4),icells,nx,ny);
        [s1_PM25,s2_PM25,s3_PM25,s4_PM25]=MAINquadrant_emi(emiTMP(:,5),icells,nx,ny);
        [s1_SO2,s2_SO2,s3_SO2,s4_SO2]=MAINquadrant_emi(emiTMP(:,6),icells,nx,ny);
        NH3_all=[s1_NH3,s2_NH3,s3_NH3,s4_NH3];
        %save('C:\data\work\projects\riat\RiatPlus-v6.1b\datasave\sa_NH3_all_fix', 'NH3_all');
        NOX_all=[s1_NOX,s2_NOX,s3_NOX,s4_NOX];
        PM10_all=[s1_PM10,s2_PM10,s3_PM10,s4_PM10];
        PM25_all=[s1_PM25,s2_PM25,s3_PM25,s4_PM25];
        SO2_all=[s1_SO2,s2_SO2,s3_SO2,s4_SO2];
        VOC_all=[s1_VOC,s2_VOC,s3_VOC,s4_VOC];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CONSIDER CASE WITH 6 OR 2 PERCURSOR EMISSIONS AS INPUT,
        %AREAL CASE. ALWAYS THERE ARE 4 QUADRANTS TO CONSIDER WIND
        %DIRETCIONS
        if areal_point==0
            %read net and select precursors
            emissioni=[];
            for i_prec=1:size(NN.PRECs,2) %areal
                if strcmp(NN.PRECs(1,i_prec),'NH3')==1
                    emissioni=[emissioni, NH3_all];
                elseif strcmp(NN.PRECs(1,i_prec),'NOX')==1
                    emissioni=[emissioni, NOX_all];
                elseif strcmp(NN.PRECs(1,i_prec),'PM10')==1
                    emissioni=[emissioni, PM10_all];
                elseif strcmp(NN.PRECs(1,i_prec),'PM25')==1
                    emissioni=[emissioni, PM25_all];
                elseif strcmp(NN.PRECs(1,i_prec),'SO2')==1
                    emissioni=[emissioni, SO2_all];
                elseif strcmp(NN.PRECs(1,i_prec),'VOC')==1
                    emissioni=[emissioni, VOC_all];
                end
            end
            
            
            %             if (size(NN.net,1)==24 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==24))% ANNlinear, 6 input
            %                 emissioni=[NH3_all,NOX_all,PM10_all,PM25_all,SO2_all,VOC_all];
            %             elseif (size(NN.net,1)==8 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==8))%ANNlinear, 2 input
            %                 emissioni=[NOX_all,VOC_all];
            %             end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %IF ALSO POINT SOURCES
        elseif areal_point==1
            %point
            [s1_NOXp,s2_NOXp,s3_NOXp,s4_NOXp]=MAINquadrant_emi(emiTMP(:,7),icells,nx,ny);
            [s1_VOCp,s2_VOCp,s3_VOCp,s4_VOCp]=MAINquadrant_emi(emiTMP(:,8),icells,nx,ny);
            [s1_NH3p,s2_NH3p,s3_NH3p,s4_NH3p]=MAINquadrant_emi(emiTMP(:,9),icells,nx,ny);
            [s1_PM10p,s2_PM10p,s3_PM10p,s4_PM10p]=MAINquadrant_emi(emiTMP(:,10),icells,nx,ny);
            [s1_PM25p,s2_PM25p,s3_PM25p,s4_PM25p]=MAINquadrant_emi(emiTMP(:,11),icells,nx,ny);
            [s1_SO2p,s2_SO2p,s3_SO2p,s4_SO2p]=MAINquadrant_emi(emiTMP(:,12),icells,nx,ny);
            NH3_allp=[s1_NH3p,s2_NH3p,s3_NH3p,s4_NH3p];
            %save('C:\data\work\projects\riat\RiatPlus-v6.1b\datasave\sa_NH3_allp_fix', 'NH3_allp');
            NOX_allp=[s1_NOXp,s2_NOXp,s3_NOXp,s4_NOXp];
            PM10_allp=[s1_PM10p,s2_PM10p,s3_PM10p,s4_PM10p];
            PM25_allp=[s1_PM25p,s2_PM25p,s3_PM25p,s4_PM25p];
            SO2_allp=[s1_SO2p,s2_SO2p,s3_SO2p,s4_SO2p];
            VOC_allp=[s1_VOCp,s2_VOCp,s3_VOCp,s4_VOCp];
            
            
            
            %create input structure, to be used in the ANNs
            %read net and select precursors
            emissioni=[];
            for i_prec=1:size(NN.PRECs,2) %areal
                if strcmp(NN.PRECs(1,i_prec),'NH3')==1
                    emissioni=[emissioni, NH3_all];
                elseif strcmp(NN.PRECs(1,i_prec),'NOX')==1
                    emissioni=[emissioni, NOX_all];
                elseif strcmp(NN.PRECs(1,i_prec),'PM10')==1
                    emissioni=[emissioni, PM10_all];
                elseif strcmp(NN.PRECs(1,i_prec),'PM25')==1
                    emissioni=[emissioni, PM25_all];
                elseif strcmp(NN.PRECs(1,i_prec),'SO2')==1
                    emissioni=[emissioni, SO2_all];
                elseif strcmp(NN.PRECs(1,i_prec),'VOC')==1
                    emissioni=[emissioni, VOC_all];
                end
            end
            for i_prec=1:size(NN.PRECs,2) %point
                
                if strcmp(NN.PRECs(1,i_prec),'NH3')==1
                    emissioni=[emissioni, NH3_allp];
                elseif strcmp(NN.PRECs(1,i_prec),'NOX')==1
                    emissioni=[emissioni, NOX_allp];
                elseif strcmp(NN.PRECs(1,i_prec),'PM10')==1
                    emissioni=[emissioni, PM10_allp];
                elseif strcmp(NN.PRECs(1,i_prec),'PM25')==1
                    emissioni=[emissioni, PM25_allp];
                elseif strcmp(NN.PRECs(1,i_prec),'SO2')==1
                    emissioni=[emissioni, SO2_allp];
                elseif strcmp(NN.PRECs(1,i_prec),'VOC')==1
                    emissioni=[emissioni, VOC_allp];
                end
            end
            
            %             if (size(NN.net,1)==48 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==48))% ANN/linear, 6 input
            %                 emissioni=[NH3_all,NOX_all,PM10_all,PM25_all,SO2_all,VOC_all,...
            %                     NH3_allp,NOX_allp,PM10_allp,PM25_allp,SO2_allp,VOC_allp];
            %             elseif (size(NN.net,1)==16 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==16))%ANNlinear, 2 input
            %                 emissioni=[NOX_all,VOC_all,NOX_allp,VOC_allp];
            %             end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %ANN INPUT DATA
        %keep only optim domain
        emissioni(find(flag_optim_dom==0),:)=[];
        %20130820 - consider only if cell completerly in PAD
        %         emissioni(find(flag_optim_dom==0 | flag_optim_dom==2),:)=[];
        input_rete2=emissioni';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %IN CASE QUADRANT EMISSIONS TOO CLOSE TO CTM BOUNDARY
        if strcmp(pathAR,'-1')==0
            AR=load(pathAR);                             %load not in optim
            AR=AR.Ratio;                                 %rename variable
            AR(flag_optim_dom==0,:)=[];                  %remove cells outside flag_optim_dom
            
            if (size(NN.net.IW{1,1},2)==48)% ANN/linear, 6 input
                ARreordrepm=repmat(AR,1,12);            %repmat
                
            elseif (size(NN.net.IW{1,1},2)==16)%ANNlinear, 2 input
                ARreordrepm=repmat(AR,1,4);            %repmat
            end
            
            %             if (size(NN.net,1)==48 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==48))% ANN/linear, 6 input
            %                 ARreordrepm=repmat(AR,1,12);            %repmat
            %             elseif (size(NN.net,1)==16 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==16))%ANNlinear, 2 input
            %                 ARreordrepm=repmat(AR,1,4);            %repmat
            %             end
            
            input_rete2=emissioni'./ARreordrepm';        %rewrite input_rete2 dividing emissions by area_ratio
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CHECK IF EMISSIONS ARE INSIDE ANNS IDENTIFICATION
        %BOUNDS...IF NOT, STOP RIAT+
        %compare ANNs input and ps_input, to see if we are in the
        %ANNs bounds
        %check in the PAD
        %         for indcel=1:ncel
        %                 if sum(input_rete2(:,indcel)<NN.ps_input.xmin)>0 | sum(input_rete2(:,indcel)>NN.ps_input.xmax)>0
        %                     find(input_rete2(:,indcel)<NN.ps_input.xmin)
        %                     find(input_rete2(:,indcel)>NN.ps_input.xmax)
        %                     figure;plot([input_rete2(:,indcel)-NN.ps_input.xmin  NN.ps_input.xmax-input_rete2(:,indcel)])
        %                     error('SR model identification bounds not respected...RIAT+ is terminated')
        %                 end
        %         end
        %matrix of min and max ANNS values
        minVal=repmat(NN.ps_input.xmin,1,size(emissioni,1));
        maxVal=repmat(NN.ps_input.xmax,1,size(emissioni,1));
        %save('C:\data\work\projects\riat\RiatPlus-v6.1b\datasave\sa_minmax_fix', 'minVal', 'maxVal');
        %check if my scenario is out of bounds
        %         indmin=find(input_rete2<minVal);  OK
        %         indmax=find(input_rete2>maxVal);  OK
        %20130828 - input_rete2>0 added to manage issue of the Alsace cells
        %adjacent to the CTM grid
        indmin=find(input_rete2<minVal & input_rete2>0);
        indmax=find(input_rete2>maxVal & input_rete2>0);
        %check bounds problem.....if limited problems, simulate the
        %scenario, otherwise exit
        %         ii
        %         jj
        %         max((input_rete2(indmax)-maxVal(indmax))./(maxVal(indmax))*100)
        %         min((input_rete2(indmin)-minVal(indmin))./minVal(indmin)*100)
        
        if length(indmax)>0
            %             if max(input_rete2(indmax)-maxVal(indmax))>1
            if max((input_rete2(indmax)-maxVal(indmax))./(maxVal(indmax))*100)>30
                
                Error=3;
                
                fprintf(fidExit, int2str(Error));
                strStatus='PROGRESSION: Aggregated scenario analysis finished.';
                disp(strStatus);
                fprintf(fidStatus, '%s\n',strStatus);
                error('SR model identification bounds not respected...RIAT+ is terminated')
            end
        end
        
        if length(indmin)>0
            %             if max(input_rete2(indmin)-minVal(indmin))<-1
            if min((input_rete2(indmin)-minVal(indmin))./minVal(indmin)*100)<-30
                
                Error=3;
                fprintf(fidExit, int2str(Error));
                strStatus='PROGRESSION: Aggregated scenario analysis finished.';
                disp(strStatus);
                fprintf(fidStatus, '%s\n',strStatus);
                error('SR model identification bounds not respected...RIAT+ is terminated')
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %IF BOUNDS RESPECTED, RUN SR MODEL
        %if constraints respected for this AQI, create results
        %define if lin o net
        if size(NN.net,1)==1
            flag_reg_net=1;
        elseif size(NN.net,1)>1
            flag_reg_net=0;
        end
        
        switch flag_reg_net
            case 0 %linear case
                aqi_per_cell=(input_rete2'*NN.net)';
            case 1 %nonlinear case
                %20111215 - run the network only on the optimization cells
                
                %20141205-ET - mapminmax replaced with a handwritten normalization
                %to deal with the problems of standalone windows executables
                ymin=NN.ps_input.ymin;
                ymax=NN.ps_input.ymax;
                xmin=NN.ps_input.xmin;
                xmax=NN.ps_input.xmax;
                G1=(ymax-ymin)./(xmax-xmin);
                
                input_rete_norm2=((input_rete2-repmat(xmin,1,size(input_rete2,2))).*repmat(G1,1,size(input_rete2,2))+repmat(ymin,size(input_rete2,1),size(input_rete2,2)));
                %                 NN.ps_input.no_change=0; %to be compatible among different matlab versions
                %                 [input_rete_norm2]=mapminmax('apply',input_rete2,NN.ps_input);
                
                %check for results outside bounds - force to be inside bounds
                input_rete_norm2(find(input_rete_norm2<-1))=-1;
                input_rete_norm2(find(input_rete_norm2>1))=1;
                
                %20141121-ET - Use new sim_exe function instead of matlab
                %sim function
                output_rete_norm2=sim_exe(NN,input_rete_norm2);
                
                %20141205-ET - mapminmax replaced with a handwritten normalization
                %to deal with the problems of standalone windows
                %executables, matlab transformation included
                G=NN.ps_target.gain;
                offst=NN.ps_target.offset;
                ymin=NN.ps_target.ymin;
                ymax=NN.ps_target.ymax;
                xmin=NN.ps_target.xmin;
                xmax=NN.ps_target.xmax;
                
                aqi_per_cell=((output_rete_norm2-ymin)*G*(xmax-xmin))/(ymax-ymin)+xmin+((xmax-xmin)*(-ymin+offst)/(ymax-ymin));
                %                 NN.ps_target.no_change=0; %to be compatible among different matlab versions
                %                 aqi_per_cell=mapminmax('reverse',output_rete_norm2,NN.ps_target);
                
                %%%
                %check for results outside bounds - force to be inside bounds
                aqi_per_cell(find(aqi_per_cell<xmin))=xmin;
                aqi_per_cell(find(aqi_per_cell>xmax))=xmax;
                %%%
                
                
        end
    end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [output_net_norm]=sim_exe(NN,input_rete_norm2)
        %20141121-ET new function to substitute matlab function sim
        sz=size(input_rete_norm2,2);
        %neural network simulation
        s_a=NN.net.IW{1,1}*input_rete_norm2+repmat(NN.net.b{1},1,sz);
        s_a1=eval(strcat(NN.net.layers{1}.transferFcn,'(s_a)'));
        s_b=NN.net.LW{2}*s_a1+repmat(NN.net.b{2},1,sz);
        output_net_norm=eval(strcat(NN.net.layers{2}.transferFcn,'(s_b)'));
        
    end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%keep track of computing time
toc

end
