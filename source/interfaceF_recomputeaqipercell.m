function [aqi_per_cell]=interfaceF_recomputeaqipercell(type, aqi_per_cell, periodIndex, aqiIndex, pathBc, flag_optim_dom, nx, ny)

if isequal(type, 'FIRSTGUESS')
    %                 sum(emissioni)
    %             aqi_per_cell=firstguess_do_aqi_per_cell(emissioni, NN, aggregationInfo, commonDataInfo, periodIndex, aqiIndex);
    fName=strtrim(pathBc(periodIndex).Bc(aqiIndex,:));
    if (strcmp(fName,'-999') == 0)
        % check compatibility of elements
        concentration=firstguess_read_Bc(fName);
        conc2=reshape(concentration',nx*ny,1);
        conc2(flag_optim_dom==0)=[];
        aqi_per_cell=conc2'-aqi_per_cell;
    end
%     else
%         aqi_per_cell=aqi_per_cell;
%     end
end

end