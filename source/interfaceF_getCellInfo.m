function [ncel, nx, ny]=interfaceF_getCellInfo(type, commonDataInfo)

%if isequal(aggregationInfo.type, 'FIRSTGUESS')
if isequal(type, 'FIRSTGUESS')
        [ncel, nx, ny]=calcCellNo('latlon', commonDataInfo);
else
        [ncel, nx, ny]=calcCellNo('utm', commonDataInfo);
end

end

