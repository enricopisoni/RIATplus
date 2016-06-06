%function [IndicEq]=EquaIndic(ic,ir,rf,nx,ny,nSc,Indic)
function [IndicEq]=EquaIndic(ic,ir,rf,nx,ny,Indic)
%IndicEq: Indicator values for constant F

%data (Indic) written as:
%x1,yn ... xn,yn
%x1,y1 ... xn,y1

% IndicEq=zeros(length(-rf:rf)*length(-rf:rf)*nSc,1);
iEq=1;

%for iSc=1:nSc %loop on scenarios
    %look for data around your reference cell.
    for jr=-rf:rf %move vertically on the matrix
        for jc=-rf:rf %move horizontally on the matrix
            if (ic+jc >= 1 && ic+jc <= nx && ir+jr >= 1 && ir+jr <= ny)%save infos only if you are in the matrix bounds
                %IndicEq(iEq,1)=Indic(ir+jr,ic+jc,iSc);
                IndicEq(iEq,1)=Indic(ir+jr,ic+jc);
                iEq=iEq+1;
            end
        end
    end
%end
end


