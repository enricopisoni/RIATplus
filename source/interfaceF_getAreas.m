function [areaIn, areaOut]=interfaceF_getAreas(type, pathArea, flag_region_dom)

if isequal(type, 'FIRSTGUESS')
    % load
    % areas=importdata(commonDataInfo.dirs.pathArea);
    areas=importdata(pathArea);
    areaIn=areas.data(:,5);
    areaOut=areas.data(:,6);
    %put 1 in the two variables, to avoid division per zero
    %     areaIn(areaIn==0)=1;
    %     areaOut(areaOut==0)=1;
else
    areaIn=repmat(1,1,length(flag_region_dom))';%repmat(1,length(flag_region_dom));
    areaOut=repmat(1,1,length(flag_region_dom))';%repmat(1,length(flag_region_dom));
end

end

