function [emis]=interfaceF_getEmisCLE(type, savedCLE, CLE, D, d, x, geometryIntermediateData, commonDataInfo)

if isequal(type, 'FIRSTGUESS')
    if size(x,1) == size(savedCLE,1)
        thisCLE=savedCLE;
    else
        thisCLE=CLE;
    end
    %             E = D * sparse(x);
    %             E = d - E;
    %             E_full = full(E);
    Escenario = D * sparse(x);
    Escenario = d - Escenario;
    %sum(D)
    
    Ecle = D * sparse(thisCLE);
    Ecle = d - Ecle;
    E_full = full(Ecle) - full(Escenario);
    %             plot(E_full)
    emis=firstguess_rebuildOrderedEmis(E_full, geometryIntermediateData, commonDataInfo);
else
    % 20160420 nn/quadrant version
    E = D * sparse(x);
    E = d - E;
    E_full = full(E);
    %for quadrant case
    %from global matrix, extract quadrant informations
    emis=interface_rebuildOrderedEmis(E_full, geometryIntermediateData, commonDataInfo);
end

end