This is a README file for MMI data compiling and processing 
which corresponds to the manuscript 
"Scalable and robust regression models for continuous proportional data"
by anonymous authors

MMI data is compiled from three different sources (1,2,3)
which are altogether joined with COMID, the unique lake identifier 

1. National Lakes Assessment 2017

Download data from https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys
A. Download MMI data
Lakes 2017 - Benthic Macroinvertebrates - nla2017_BentMMI.csv
B. Download site info data (contains lake identifier COMID)
Lakes 2017	Site Information - nla_2017_site_information-data.csv
A and B are joined with SITE_ID column


2. LakeCat database

Lake watershed covariate can be directly downloaded from
https://www.epa.gov/national-aquatic-resource-surveys/lakecat-web-tool-table-view
Instructions: 
A. Select request type - conterminous US
B. Select metric names - 
  * agkffact
  * bfi
  * cbnf
  * conif
  * pctcrop2016, pcthay2016: sum becomes crophay
  * fert
  * manure
  * pestic1997
  * pcturbmd2016, pcturbhi2016: sum becomes urbmdhi
C. Select area of interest - Watershed
D. Click "generate csv". 


3. NHDPlusV2: lake coordinates and area information

Download gdb data (NHDPlusV21_National_Seamless_Flattened_Lower48.gdb) from 
https://www.epa.gov/waterdata/nhdplus-national-data
NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07.7z
A. lake coordinate is extracted from centroid of lake polygon (sfc_MULTIPOLYGON) object
B. lake area is available in LakeArea column.




