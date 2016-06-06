function [ D ] = interface_get_build_solution_D_indexed(refInfo,index1,index2)

stringToRun='';
stringToRun=strcat(stringToRun, refInfo.get_build_solution_D_indexedFunction);
stringToRun=strcat(stringToRun, '(');
stringToRun=strcat(stringToRun,'refInfo,');
stringToRun=strcat(stringToRun,'index1,');
stringToRun=strcat(stringToRun,'index2');
stringToRun=strcat(stringToRun,')');
[D]=eval(stringToRun);

end