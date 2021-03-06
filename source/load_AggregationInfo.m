function [aggregationInfo]=load_AggregationInfo(confFile)

 [type,configurationFiles]=init_ConfigurationVar(confFile);
 
 commonDataInfo=load_AggregationDataSetting(configurationFiles.commonInfoFile);
 geometryDataInfo=load_AggregationDataSetting(configurationFiles.geometryInfoFile);
 mathDataInfo=load_AggregationDataSetting(configurationFiles.mathInfoFile);
 domainApplyTypeDataInfo=load_AggregationDataSetting(configurationFiles.domainApplyTypeInfoFile);
 fifthConfDataInfo=load_AggregationDataSetting(configurationFiles.fifthConfInfoFile);
 sixthConfDataInfo=load_AggregationDataSetting(configurationFiles.sixthConfInfoFile);

 aggregationInfo.type=type;
 aggregationInfo.commonDataInfo=commonDataInfo;
 aggregationInfo.geometryDataInfo=geometryDataInfo;
 aggregationInfo.mathDataInfo=mathDataInfo;
 aggregationInfo.domainApplyTypeDataInfo=domainApplyTypeDataInfo;
 aggregationInfo.fifthConfDataInfo=fifthConfDataInfo;
 aggregationInfo.sixthConfDataInfo=sixthConfDataInfo;

end
