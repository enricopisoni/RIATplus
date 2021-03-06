%emisType 0 areal, 1 point, 2 both
function [ finalEmis ] = firstguess_buildEmission(emis,emisP,NN,bcSet, refInfo,refData,emissionSource,periodIndex,aqiIndex)

finalEmis=[];
%precNames=refData.nnSuperSet(periodIndex).nnSet(aqiIndex).PRECs;
% ToDO: modify FG_read
% CHECK size function on REAL precNames var 
precNames={'NOx';'NMVOC';'NH3';'PM10';'PM25';'SOx'};
% 20160419 Add two flags (Delta/Aggregated) inside SR/First Guess
isDelta=0==1 ;
isAggregated=0==1;
%isDelta=strcmp(NN.Class,'Delta')==1 ;
%isAggregated=strcmp(NN.ArPt,'Aggregated')==1;
%save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\var_buildEmission_new','emis','emisP');

for i_prec=1:size(precNames,1) %areal do always
    if strcmp(precNames(i_prec),'NH3')==1
        finalEmis=[finalEmis, emis.NH3_all];
    elseif strcmp(precNames(i_prec),'NOx')==1
        finalEmis=[finalEmis, emis.NOX_all];
    elseif strcmp(precNames(i_prec),'PM10')==1
        finalEmis=[finalEmis, emis.PM10_all];
    elseif strcmp(precNames(i_prec),'PM25')==1
        finalEmis=[finalEmis, emis.PM25_all];
    elseif strcmp(precNames(i_prec),'SOx')==1
        finalEmis=[finalEmis, emis.SO2_all];
    elseif strcmp(precNames(i_prec),'NMVOC')==1
        finalEmis=[finalEmis, emis.VOC_all];
    end
end

if (emissionSource == 'P')
    for i_prec=1:size(precNames,1) %areal
        if strcmp(precNames(i_prec),'NH3')==1
            finalEmis=[finalEmis, emisP.NH3_all];
        elseif strcmp(precNames(i_prec),'NOx')==1
            finalEmis=[finalEmis, emisP.NOX_all];
        elseif strcmp(precNames(i_prec),'PM10')==1
            finalEmis=[finalEmis, emisP.PM10_all];
        elseif strcmp(precNames(i_prec),'PM25')==1
            finalEmis=[finalEmis, emisP.PM25_all];
        elseif strcmp(precNames(i_prec),'SOx')==1
            finalEmis=[finalEmis, emisP.SO2_all];
        elseif strcmp(precNames(i_prec),'NMVOC')==1
            finalEmis=[finalEmis, emisP.VOC_all];
        end
    end
end

%save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPointEmissioni_step1_new','finalEmis');

if isAggregated
    finalEmis=finalEmis(:,1:(size(finalEmis,2)/2))+finalEmis(:,(size(finalEmis,2)/2)+1:size(finalEmis,2));
end
%save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPointEmissioni_step2_new','finalEmis');

% build right ordered emissions...
%refInfo.isAggregated;
%refInfo.isDelta;

%read net and sum Areal and Point
%if strcmp(NN.ArPt,'Aggregated')==1

%read net and create Delta emissions
%emi delta
%if strcmp(NN.Class,'Delta')==1

%emi delta
if isDelta
    flag_optim_dom=commonDataInfo.domainInfo.flag_optim_dom;
    flag_opt_filt=flag_optim_dom(flag_region_dom==1 | flag_region_dom==2,1);
    %20140403 ET - Filtering basecase cells
    BCemi=BCset.emi_bc(flag_optim_dom==1 | flag_optim_dom==2,:);
    BCconc=BCset.conc_bc(flag_optim_dom==1 | flag_optim_dom==2,:);
    emissioni=emissioni(flag_opt_filt==1 | flag_opt_filt==2,:);
    
    emissioni2=(emissioni-BCemi)./BCemi;
    finalEmis=emissioni2;
    test_nan=find(isnan(finalEmis));
    finalEmis(test_nan)=0.;
end



%save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\emissioni_new','finalEmis');

end