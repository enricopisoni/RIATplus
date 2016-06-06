function [ setInfo  ] = interface_getOrderedPoll( commonData, refInfo)

stringToRun='';
stringToRun=strcat(stringToRun, refInfo.getOrderedPollFunction);
stringToRun=strcat(stringToRun, '(');
stringToRun=strcat(stringToRun,'commonData');
stringToRun=strcat(stringToRun,'refInfo');
stringToRun=strcat(stringToRun,')');
[setInfo]=eval(stringToRun);