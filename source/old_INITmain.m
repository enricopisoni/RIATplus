%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%20130830-EP
% - added the words "PROGRESSION: " in the beginning of command window messages

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
%domain is defined by flag_region_dom==1 and flag_region_dom==2

%for border cells, in "MAINbase_emi" the contribution of the part of the
%cell that is outside the optimization domain, is added to the base_emi_low
%and base_emi_low

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NB: the INIT is able to create new D and d matrices, processing the emission
%DB (first line read in file flag_optim_path, using icells that is read from ANNs files
% (line 10-15 file flag_optim_path), and creating the
% D,d file using names from first 6 lines of Dd matrices (line 28-33
% flag_optim_path)
% Even if you process a different database (not yearly, but seasonal i.e.,
% you need to use these lines to define ANNs paths and D,d filenames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This version of the code is running both on linux and matlab
function []=old_INITmain(f1)
%#function network

%to run the routine, write
% INITmain('inputOPERAaspa/flag_optim_path.txt');

%20111206 generalization of quadrants: loading the ANNs, the tool loads the
%dimension of quadrants (icells) contained in the ANNs
%So D and d are computed, depending on the dimension of quadrants (icells)
%, that can be different for each considered pollutant

%load configuration files
% f1='inputOPERAaspa/flag_optim_path.txt';
% f2='input/POSTflag_optim_path.txt';

%to keep track of computing time
tic

%remove warning
%needed otherwise the compile toolbox is not working
warning off MCR:ClassInUseAtExit
warning off MATLAB:ClassInstanceExists
warning('off','all')
%read file input/flag_optim_path
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
flag_reg_net         =flag_optim_oth(2); % DUMMY VARIABLE - NOT READ HERE. IN RIAT, THIS WAS 0 means you consider linear model, 1 neural network model, for all AQIs.
flag_constraints     =flag_optim_oth(3); % flag_constraints=0 means all techs are not replaceable (keep LB as they are).
flag_mode_ce_mo      =flag_optim_oth(4); % 0 means multi-objective, 1 means cost-effectiveness.
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
if aqi_weights(1)==1
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

%output file name
% nomeOUTPUT=pathOUF;

%create directory name
mkdir('output')

%load NETWORKS
nnSuperSet = struct('nnSet', {});
for h=1:3,
    nnSetHorizon = struct('ps_input', {}, 'ps_target', {}, 'net', {}, 'icells', {},'Class', {}, 'PRECs', {},'ps_pca',{},'ArPt',{});
    for i=1:AQINum,
        if (isequal(strtrim(pathANN(h).ANNs(i,:)),'-999')==0)
            %nn=load(strtrim((pathANN(h).ANNs(i,:))));
            [nn]=net_read(strtrim((pathANN(h).ANNs(i,:)))); %20140617 ET - read model from txt file

        else
            nn = struct('ps_input', {-999}, 'ps_target', {-999}, 'net', {-999}, 'icells', {-999},'Class', {-999}, 'PRECs', {-999},'ps_pca',{-999},'ArPt',{-999});
        end
        nnSetHorizon= [nnSetHorizon, nn];
    end
    nnSuperSet = [nnSuperSet, struct('nnSet', nnSetHorizon)];
end


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

%COLUMNS OF OUT_CLEMFR FILE
%macrosector sector activity tecnology lowHighEmissionFlag
%removalEfficiencies (NOx,VOC,NH3,PM10,PM25,SO2)
%unabatedEmissionFactores (NOx,VOC,NH3,PM10,PM25,SO2)
%unitcost arCLE_year arPOT_year flagReplace flagOptim flagNTM_TM arADSdet
%atREFyear


%COLUMNS OF OUT_CLEMFR FILE
%macrosector sector activity tecnology lowHighEmissionFlag
%removalEfficiencies (NOx,VOC,NH3,PM10,PM25,SO2)
%unabatedEmissionFactores (NOx,VOC,NH3,PM10,PM25,SO2)
%unitcost arCLE_year arPOT_year flagReplace flagOptim ar2005

%hl_flag=1 (low), =2(high);
%replace=1 (not replaceable)
%optim=1 (to be used in the optimization)

%NOTE: AR and RE are in %

%rows are quadruples * 2 (low and high). Quadruples are related to
%GAINS activities, without NOC

%select optimization domain: %1 means cell in optimization domain, 0 means
%outside (2 is border cell, partly inside and partly outside)
%UNIBS(ET)20131001 - added new column for regional domain
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
%EPcorrection
stepsize=round(max(abs(coordinate(2,1)-coordinate(1,1)),abs(coordinate(2,2)-coordinate(1,2))));

nxny=round((max(coordinate)-min(coordinate))/stepsize+1);
nx=nxny(1,1);
ny=nxny(1,2);
ncel=nx*ny;

%FROM HERE INSERT NEW LOOP
%new pathEMI updated
pathEMIyea=strcat(pathEMI,'TP1/');
pathEMIwin=strcat(pathEMI,'TP2/');
pathEMIsum=strcat(pathEMI,'TP3/');

pathVEC={pathEMIyea,pathEMIwin,pathEMIsum};

%if summer and winter directory do not exist, only create yearly matrices
%of precomputed emissions
if (exist(pathVEC{2},'dir')==0 & (exist(pathVEC{3},'dir')==0))
    indini=1;
    indfini=1;
    %if both summer and winter directory exist, create summer and winter
    %precomputed emissions preprocessing data, and yearl preprocessed
    %emissions as the sum of winter and summer
elseif (exist(pathVEC{2},'dir')==7 & (exist(pathVEC{3},'dir')==7))
    indini=1;
    indfini=3;
else
    error('No data available for emission preprocessing')
end

for k=indini:indfini %year, winter, summmer
    
    emi={};
    pathEMI=pathVEC{k};
    coordinate=domainData.data(:,1:2);
    
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
    
    
    
    
    
    %     stepsize=max(coordinate(2,1)-coordinate(1,1),coordinate(2,2)-coordinate(1,2));
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
    [emi,~]=MAINload_emi(ncel,pathEMI,flag_region_dom,flag_mode_ce_mo);
    
    %load PM10 to PM25 relationship
    %     pm10aveToExceed=load(pathPM);
    
    %in case of RIAT, matrices and input are cut and flipped (due to an initial
    %problem in the data). Otherwise no change to the data is performed
    if coordinateFlip==1
        ny=ny-3; %problem in RIAT about number of cells
        ncel=nx*ny;
        [dataerase]=find(coordinate(:,1)>874.28 | coordinate(:,2)>5210.7);
        flag_region_dom(dataerase,:)=[];
        coordinate(dataerase,:)=[];
        pop(dataerase,:)=[];
        flag_aqi_dom(dataerase,:)=[];
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
        %new regional flag
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
    
    %imagesc(flipud(reshape(flag_region_dom,ny,nx))) % spatial map
    if coordinateFlip==2 %ASPA CASE
        flag_region_dom1=reshape(flag_region_dom,ny,nx);
        flag_region_dom2=reshape(flipud(flag_region_dom1),nx*ny,1);
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
    ncelopt=length(find(flag_region_dom==1 | flag_region_dom==2));
    
    %load or compute base emissions for each cell sum of emissions, for each cell,
    %separated for low and high rows extracte from emi (one row per each sec-act) and then summed
    %emissions of the part of the border cell outside optimziation domain, is
    %summed up to base_emi_low and base_emi_high
    %actlev_final=activity level of all sector-activity-techs
    %actlev_final_sum=optimization domain wide sum of activity levels
    %UNIBS(ET)20131002 - added flag_region_dom in input
    [actlev_final,actlev_final_sum,base_emi_low,base_emi_high,base_emi_low_noc,...
        base_emi_high_noc]=MAINbase_emi(emi,global_data,flag_optim_dom,flag_region_dom,global_data_noc);
    
    %create matrices D and d, for quadrant computation
    %in case of areal+point emissions, only D and d are created
    %in case ot areal and point emissions separated, D, d, and also Dp and dp are created
    %in this second cass, the matrices contain only the relative areal or point
    %contribution (D and d areal, Dp and dp point)
    %CREATE D and d only if ANN is available for the considered AQI
    for indaqi = 1:size(pathANN(1).ANNs,1), %aqis
        if (isequal(strtrim(pathANN(k).ANNs(indaqi,:)),'-999')==0) %check if AQI is defined
            if exist(strcat((strtrim(pathDd(k).Dd(indaqi,:))),'.mat'))==0 %check if Dd is already existing
                
                disp(strcat('PROGRESSION: Creating preprocessed emissions, for AQI -->',strtrim(pathANN(k).ANNs(indaqi,:))))
                [D,d,Dp,dp]=INITquadrant(emi,global_data,...
                    base_emi_low, base_emi_high, base_emi_low_noc,...
                    base_emi_high_noc,flag_region_dom,nx,ny,strtrim(pathANN(k).ANNs(indaqi,:)),areal_point);
                save(strtrim(pathDd(k).Dd(indaqi,:)),'D','d','Dp','dp');
            else
                disp(strcat('PROGRESSION: Preprocessed emissions already created, for AQI -->',strtrim(pathANN(k).ANNs(indaqi,:))))
            end
        end
    end
end

toc
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DDsparse,ddsparse,DDsparsep,ddsparsep]=INITquadrant(emi,global_data,...
    base_emi_low, base_emi_high, base_emi_low_noc,...
    base_emi_high_noc,flag_region_dom,nx,ny,pathANNpm10,areal_point)

%in case of areal+point emissions, only DDsparse and ddsparse contain
%useful informations (the DDsparsep and ddsparsep are empty matrices)

%routine to create matrices D and d, to be used to compute
%quadrant emissions, applying:
% d - D * x
DDsparsep=[];
ddsparsep=[];
DDsparse=[];
ddsparse=[];
%to define icells, load pathANNpm10
%#function network
[net]=net_read(pathANNpm10); %20140617 ET - read model from txt file
%net=load(pathANNpm10);
icells=net.icells;

%defininf variables
np=size(base_emi_low,2);    %number of pollutants
nc=length(flag_region_dom);  %number of cells
ncv=size(global_data,1);    %number of control variables
re=global_data(:,6:11)/100; %removal efficiencies
quad=4;

%info to cut ddsparse and ddsparsep, considering only optimization cells
%this clipping allows to speed up the procedure
indopt=flag_region_dom; %cells in optimization domain
indoptrep=repmat(indopt,quad*np,1); %index repeated for pollutants and quadrants
indiciMAT=reshape(indopt,ny,nx);

%initialize matrices
emi_red(1:nc,1:np,1:ncv)=0;
emi_rem(1:nc,1:np)=0;


%loop over cells to create emission matrices
for i=1:length(flag_region_dom)
    
    %case of cells without optimization domain
    %emission order: NOX, COV, NH3, PM10, PM25, SO2.
    if flag_region_dom(i)==0
        
        %INITIAL EMISSIONS
        %case areal and point summed:
        if areal_point==0
            emi_rem(i,:)=base_emi_low(i,:)+base_emi_high(i,:);
            %case areal and point separated
        elseif areal_point==1
            emi_rem_low(i,:)=base_emi_low(i,:);
            emi_rem_high(i,:)=base_emi_high(i,:);
        end
        
    else     %for cells inside optimization domain
        
        %EMISSIONS REDUCED (BOTH "AREAL+POINT" AND "AREAL AND POINT")
        %define emission matrix
        emi_red_tmp{i}=emi{i,2}(1:ncv,1:6);
        %compute emissions*RE
        emi_red(i,1:np,1:ncv)=(emi_red_tmp{i}.*re)';
        
        %INITIAL EMISSIONS
        %compute final component of D matrix
        if areal_point==0
            emi_rem(i,:)=base_emi_low(i,:)+base_emi_high(i,:)+...
                base_emi_low_noc(i,:)+base_emi_high_noc(i,:);
            %case areal and point separated
        elseif areal_point==1
            emi_rem_low(i,:)=base_emi_low(i,:)+...
                base_emi_low_noc(i,:);
            emi_rem_high(i,:)=base_emi_high(i,:)+...
                base_emi_high_noc(i,:);
        end
    end
end

%COMPUTE d: INITIAL EMISSIONS
%grid emissions are as NOX, COV, NH3, PM10, PM25, SO2
%areal+point emissions
if areal_point==0
    %compute quadrants
    [s1_NOX,s2_NOX,s3_NOX,s4_NOX]=INITquadrant_emi(emi_rem(:,1),icells,nx,ny,'NOX',indiciMAT);
    [s1_VOC,s2_VOC,s3_VOC,s4_VOC]=INITquadrant_emi(emi_rem(:,2),icells,nx,ny,'VOC',indiciMAT);
    [s1_NH3,s2_NH3,s3_NH3,s4_NH3]=INITquadrant_emi(emi_rem(:,3),icells,nx,ny,'NH3',indiciMAT);
    [s1_PM10,s2_PM10,s3_PM10,s4_PM10]=INITquadrant_emi(emi_rem(:,4),icells,nx,ny,'PM10',indiciMAT);
    [s1_PM25,s2_PM25,s3_PM25,s4_PM25]=INITquadrant_emi(emi_rem(:,5),icells,nx,ny,'PM25',indiciMAT);
    [s1_SO2,s2_SO2,s3_SO2,s4_SO2]=INITquadrant_emi(emi_rem(:,6),icells,nx,ny,'SO2',indiciMAT);
    %create matrix d
    dd=[s1_NOX;s2_NOX;s3_NOX;s4_NOX;s1_VOC;s2_VOC;s3_VOC;s4_VOC;...
        s1_NH3;s2_NH3;s3_NH3;s4_NH3;s1_PM10;s2_PM10;s3_PM10;s4_PM10;...
        s1_PM25;s2_PM25;s3_PM25;s4_PM25;s1_SO2;s2_SO2;s3_SO2;s4_SO2];
    %create sparse matrix
    ddsparse=sparse(dd);
    %keep optimization cells (inside optimization==1 and border cells==2)
    ddsparse=ddsparse((indoptrep==1 | indoptrep==2),:);
    clearvars dd s1_NOX s2_NOX s3_NOX s4_NOX s1_VOC s2_VOC s3_VOC s4_VOC ...
        s1_NH3 s2_NH3 s3_NH3 s4_NH3 s1_PM10 s2_PM10 s3_PM10 s4_PM10 ...
        s1_PM25 s2_PM25 s3_PM25 s4_PM25 s1_SO2 s2_SO2 s3_SO2 s4_SO2
    
    %areal and point emissions separated
elseif areal_point==1
    %areal
    [s1_NOX,s2_NOX,s3_NOX,s4_NOX]=INITquadrant_emi(emi_rem_low(:,1),icells,nx,ny,'NOx_d',indiciMAT);
    [s1_VOC,s2_VOC,s3_VOC,s4_VOC]=INITquadrant_emi(emi_rem_low(:,2),icells,nx,ny,'VOC_d',indiciMAT);
    [s1_NH3,s2_NH3,s3_NH3,s4_NH3]=INITquadrant_emi(emi_rem_low(:,3),icells,nx,ny,'NH3_d',indiciMAT);
    [s1_PM10,s2_PM10,s3_PM10,s4_PM10]=INITquadrant_emi(emi_rem_low(:,4),icells,nx,ny,'PM10_d',indiciMAT);
    [s1_PM25,s2_PM25,s3_PM25,s4_PM25]=INITquadrant_emi(emi_rem_low(:,5),icells,nx,ny,'PM25_d',indiciMAT);
    [s1_SO2,s2_SO2,s3_SO2,s4_SO2]=INITquadrant_emi(emi_rem_low(:,6),icells,nx,ny,'SO2_d',indiciMAT);
    %par1_Right=[s1_NOX;s2_NOX;s3_NOX;s4_NOX];
    %%save('C:\data\work\projects\riat\s1_NOX', 's1_NOX');
    %point
    [s1_NOXp,s2_NOXp,s3_NOXp,s4_NOXp]=INITquadrant_emi(emi_rem_high(:,1),icells,nx,ny,'NOx_dp',indiciMAT);
    [s1_VOCp,s2_VOCp,s3_VOCp,s4_VOCp]=INITquadrant_emi(emi_rem_high(:,2),icells,nx,ny,'VOC_dp',indiciMAT);
    [s1_NH3p,s2_NH3p,s3_NH3p,s4_NH3p]=INITquadrant_emi(emi_rem_high(:,3),icells,nx,ny,'NH3p',indiciMAT);
    [s1_PM10p,s2_PM10p,s3_PM10p,s4_PM10p]=INITquadrant_emi(emi_rem_high(:,4),icells,nx,ny,'PM10_dp',indiciMAT);
    [s1_PM25p,s2_PM25p,s3_PM25p,s4_PM25p]=INITquadrant_emi(emi_rem_high(:,5),icells,nx,ny,'PM25_dp',indiciMAT);
    [s1_SO2p,s2_SO2p,s3_SO2p,s4_SO2p]=INITquadrant_emi(emi_rem_high(:,6),icells,nx,ny,'SO2_dp',indiciMAT);
    %create matrix d
    dd=[s1_NOX;s2_NOX;s3_NOX;s4_NOX;s1_VOC;s2_VOC;s3_VOC;s4_VOC;...
        s1_NH3;s2_NH3;s3_NH3;s4_NH3;s1_PM10;s2_PM10;s3_PM10;s4_PM10;...
        s1_PM25;s2_PM25;s3_PM25;s4_PM25;s1_SO2;s2_SO2;s3_SO2;s4_SO2];
    ddp=[s1_NOXp;s2_NOXp;s3_NOXp;s4_NOXp;s1_VOCp;s2_VOCp;s3_VOCp;s4_VOCp;...
        s1_NH3p;s2_NH3p;s3_NH3p;s4_NH3p;s1_PM10p;s2_PM10p;s3_PM10p;s4_PM10p;...
        s1_PM25p;s2_PM25p;s3_PM25p;s4_PM25p;s1_SO2p;s2_SO2p;s3_SO2p;s4_SO2p];
    %create sparse matrix
	%save('C:\data\work\projects\riat\all_bef_sparse_right', 'dd');
    ddsparse=sparse(dd);
    ddsparsep=sparse(ddp);
    %keep optimization cells (inside optimization==1 and border cells==2)
	%save('C:\data\work\projects\riat\all_bef_filter_right', 'ddsparse');
    ddsparse=ddsparse((indoptrep==1 | indoptrep==2),:);
    ddsparsep=ddsparsep((indoptrep==1 | indoptrep==2),:);
	%save('C:\data\work\projects\riat\all_final_right', 'ddsparse');
    
    clearvars dd ddp s1_NOX s2_NOX s3_NOX s4_NOX s1_VOC s2_VOC s3_VOC s4_VOC ...
        s1_NH3 s2_NH3 s3_NH3 s4_NH3 s1_PM10 s2_PM10 s3_PM10 s4_PM10 ...
        s1_PM25 s2_PM25 s3_PM25 s4_PM25 s1_SO2 s2_SO2 s3_SO2 s4_SO2 ...
        s1_NOXp s2_NOXp s3_NOXp s4_NOXp s1_VOCp s2_VOCp s3_VOCp s4_VOCp ...
        s1_NH3p s2_NH3p s3_NH3p s4_NH3p s1_PM10p s2_PM10p s3_PM10p s4_PM10p ...
        s1_PM25p s2_PM25p s3_PM25p s4_PM25p s1_SO2p s2_SO2p s3_SO2p s4_SO2p
end


%NOW START TO CONSIDER EMISSIONS THAT ARE REMOVED DUE TO APPLICATION OF TECHNOLOGIES

%define D matrix as empty matrix (only considering PAD cells)
DD=zeros(size(find(indoptrep==1 | indoptrep==2),1),size(global_data,1));
if areal_point==1
    DDp=zeros(size(find(indoptrep==1 | indoptrep==2),1),size(global_data,1));
end

%loop to create matrix D
%areal+point emissions
if areal_point==0
    for i=1:size(global_data,1)
        %         i
        [s1_NOX_red,s2_NOX_red,s3_NOX_red,s4_NOX_red]=INITquadrant_emi(emi_red(:,1,i),icells,nx,ny,'NOx_D',indiciMAT);
        [s1_VOC_red,s2_VOC_red,s3_VOC_red,s4_VOC_red]=INITquadrant_emi(emi_red(:,2,i),icells,nx,ny,'VOC_D',indiciMAT);
        [s1_NH3_red,s2_NH3_red,s3_NH3_red,s4_NH3_red]=INITquadrant_emi(emi_red(:,3,i),icells,nx,ny,'NH3_D',indiciMAT);
        [s1_PM10_red,s2_PM10_red,s3_PM10_red,s4_PM10_red]=INITquadrant_emi(emi_red(:,4,i),icells,nx,ny,'PM10_D',indiciMAT);
        [s1_PM25_red,s2_PM25_red,s3_PM25_red,s4_PM25_red]=INITquadrant_emi(emi_red(:,5,i),icells,nx,ny,'PM25_D',indiciMAT);
        [s1_SO2_red,s2_SO2_red,s3_SO2_red,s4_SO2_red]=INITquadrant_emi(emi_red(:,6,i),icells,nx,ny,'SO2_D',indiciMAT);
        %full list of results
        tmp=[s1_NOX_red;s2_NOX_red;s3_NOX_red;s4_NOX_red;s1_VOC_red;s2_VOC_red;s3_VOC_red;s4_VOC_red;...
            s1_NH3_red;s2_NH3_red;s3_NH3_red;s4_NH3_red;s1_PM10_red;s2_PM10_red;s3_PM10_red;s4_PM10_red;...
            s1_PM25_red;s2_PM25_red;s3_PM25_red;s4_PM25_red;s1_SO2_red;s2_SO2_red;s3_SO2_red;s4_SO2_red];
        %only extract PAD cells
        DD(:,i)=tmp((indoptrep==1 | indoptrep==2));
    end
    %create sparse matrix
    DDsparse=sparse(DD);
    %keep optimization cells (inside optimization==1 and border cells==2)
    %     DDsparse=DDsparse((indoptrep==1 | indoptrep==2),:);
    
    clearvars DD s1_NOX s2_NOX s3_NOX s4_NOX s1_VOC s2_VOC s3_VOC s4_VOC ...
       s1_NH3 s2_NH3 s3_NH3 s4_NH3 s1_PM10 s2_PM10 s3_PM10 s4_PM10 ...
       s1_PM25 s2_PM25 s3_PM25 s4_PM25 s1_SO2 s2_SO2 s3_SO2 s4_SO2 tmp
    
    %areal and point emissions separated
elseif areal_point==1
    
    %     DDp=zeros(size(ddp,1),size(global_data,1));
    for i=1:size(global_data,1)
        [i size(global_data,1)];  %display status
        [s1_NOX_red,s2_NOX_red,s3_NOX_red,s4_NOX_red]=INITquadrant_emi(emi_red(:,1,i),icells,nx,ny,'NOx_D',indiciMAT);
        [s1_VOC_red,s2_VOC_red,s3_VOC_red,s4_VOC_red]=INITquadrant_emi(emi_red(:,2,i),icells,nx,ny,'VOC_D',indiciMAT);
        [s1_NH3_red,s2_NH3_red,s3_NH3_red,s4_NH3_red]=INITquadrant_emi(emi_red(:,3,i),icells,nx,ny,'NH3_D',indiciMAT);
        [s1_PM10_red,s2_PM10_red,s3_PM10_red,s4_PM10_red]=INITquadrant_emi(emi_red(:,4,i),icells,nx,ny,'PM10_D',indiciMAT);
        [s1_PM25_red,s2_PM25_red,s3_PM25_red,s4_PM25_red]=INITquadrant_emi(emi_red(:,5,i),icells,nx,ny,'PM25_D',indiciMAT);
        [s1_SO2_red,s2_SO2_red,s3_SO2_red,s4_SO2_red]=INITquadrant_emi(emi_red(:,6,i),icells,nx,ny,'SO2_D',indiciMAT);
        
        if global_data(i,5)==1 %areal emissions
            tmp1=[s1_NOX_red;s2_NOX_red;s3_NOX_red;s4_NOX_red;s1_VOC_red;s2_VOC_red;s3_VOC_red;s4_VOC_red;...
                s1_NH3_red;s2_NH3_red;s3_NH3_red;s4_NH3_red;s1_PM10_red;s2_PM10_red;s3_PM10_red;s4_PM10_red;...
                s1_PM25_red;s2_PM25_red;s3_PM25_red;s4_PM25_red;s1_SO2_red;s2_SO2_red;s3_SO2_red;s4_SO2_red];
            DD(:,i)=tmp1((indoptrep==1 | indoptrep==2));
            
        elseif global_data(i,5)==2 %point emissions
            tmp2=[s1_NOX_red;s2_NOX_red;s3_NOX_red;s4_NOX_red;s1_VOC_red;s2_VOC_red;s3_VOC_red;s4_VOC_red;...
                s1_NH3_red;s2_NH3_red;s3_NH3_red;s4_NH3_red;s1_PM10_red;s2_PM10_red;s3_PM10_red;s4_PM10_red;...
                s1_PM25_red;s2_PM25_red;s3_PM25_red;s4_PM25_red;s1_SO2_red;s2_SO2_red;s3_SO2_red;s4_SO2_red];
            DDp(:,i)=tmp2((indoptrep==1 | indoptrep==2));
        end
        
        clearvars tmp1 tmp2 s1_NOX s2_NOX s3_NOX s4_NOX s1_VOC s2_VOC s3_VOC s4_VOC ...
        s1_NH3 s2_NH3 s3_NH3 s4_NH3 s1_PM10 s2_PM10 s3_PM10 s4_PM10 ...
        s1_PM25 s2_PM25 s3_PM25 s4_PM25 s1_SO2 s2_SO2 s3_SO2 s4_SO2
        
    end
    %create sparse matrix
    DDsparse=sparse(DD);
    DDsparsep=sparse(DDp);
%     %keep optimization cells (inside optimization==1 and border cells==2)
%         DDsparse=DDsparse((indoptrep==1 | indoptrep==2),:);
%         DDsparsep=DDsparsep((indoptrep==1 | indoptrep==2),:);
end

%cut ddsparse and ddsparsep


disp('PROGRESSION: matrices d and D created');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spi1,spi2,spi3,spi4]=INITquadrant_emi(valori_emi,r,nx,ny,inq_deb,indiciMAT)

%compute quadrant emissions

%domain dimensions and other variables

dimx=nx;
dimy=ny;
prova_emi=valori_emi(:,1);
emi=reshape(prova_emi,dimy,dimx);

%enlarged emission, with increased zeros all around (used to be able to
%perform simple products among emi and factor matrix
emi_enlarged=zeros(size(emi,1)+2*(r-2),size(emi,2)+2*(r-2));
ZeroMap=emi_enlarged;
emi_enlarged(r-2+1:r-2+1+size(emi,1)-1,r-2+1:r-2+1+size(emi,2)-1)=emi;
[sa sb]=size(emi_enlarged);
indiciMAT_enlarged=zeros(size(emi,1)+2*(r-2),size(emi,2)+2*(r-2));
indiciMAT_enlarged(r-2+1:r-2+1+size(indiciMAT,1)-1,r-2+1:r-2+1+size(indiciMAT,2)-1)=indiciMAT;
%
%left quadrant
% %factors to perform products of cell emissions
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
% %factors to perform products of cell emissions
 fattoreR=rot90(fattore,2);

%initialize variable, and set dimensions
spicchio4=zeros(dimy,dimx);

%up quadrant
% %factors to perform products of cell emissions
 fattoreU=rot90(fattore,1);

%initialize variable, and set dimensions
spicchio2=zeros(dimy,dimx);

% % % initialize fattore
%  fattoreS1=zeros(((sa-r)*((sa-r)-r+2-1))+((sb-r)-r+2),dimy*dimx);
%  fattoreS2=fattoreS1;
%  fattoreS3=fattoreS1;
%  fattoreS4=fattoreS1;

if sum(valori_emi)~=0
    
    %loop to aggregate emissions
    %create index matrix
    %INDx=[];
    INDx=zeros(size(r+1:sb-r,2)*((sa-r)-(r)),2);
    cic=1;
    for iv=r+1:(sa-r)
        
        i1=(size(r+1:sb-r,2)*(cic-1))+1;
        i2=((size(r+1:sb-r,2)*cic));
        INDx(i1:i2,:)=[ones(size(r+1:sb-r,2),1)*iv,(r+1:sb-r)'];
        cic=cic+1;
    end
    
    for i=1:size(INDx,1)
        if indiciMAT_enlarged(INDx(i,1),INDx(i,2))>0
            %left quadrant
            %spicchio3(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1)-r:INDx(i,1)+r,INDx(i,2)-r:INDx(i,2)))); %rectangle
            spicchio3(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1)-r:INDx(i,1)+r,INDx(i,2)-r:INDx(i,2)).*fattore)); %quadrant

            %up quadrant (down as matrix, up geographically speaking)
            %spicchio1(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1)-r:INDx(i,1),INDx(i,2)-r:INDx(i,2)+r))); %rectangle
            spicchio1(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1)-r:INDx(i,1),INDx(i,2)-r:INDx(i,2)+r).*fattoreD));

            %right quadrant
            %spicchio4(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1)-r:INDx(i,1)+r,INDx(i,2):INDx(i,2)+r)));%rectangle
            spicchio4(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1)-r:INDx(i,1)+r,INDx(i,2):INDx(i,2)+r).*fattoreR));

            %down quadrant (up as matrix, down geographically speaking)
            %spicchio2(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1):INDx(i,1)+r,INDx(i,2)-r:INDx(i,2)+r))); %rectangle
            spicchio2(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1):INDx(i,1)+r,INDx(i,2)-r:INDx(i,2)+r).*fattoreU));

            %areal sum
            %spicchio(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1)-r:INDx(i,1)+r,INDx(i,2)-r:INDx(i,2)+r)));
        end
    end
    
    
end

%ORIGINAL
spi1=reshape(spicchio1,dimx*dimy,1);
spi2=reshape(spicchio2,dimx*dimy,1);
spi3=reshape(spicchio3,dimx*dimy,1);
spi4=reshape(spicchio4,dimx*dimy,1);

par1_Right=[spi1;spi2;spi3;spi4];
%%save('C:\data\work\projects\riat\par1_Right', 'par1_Right');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%