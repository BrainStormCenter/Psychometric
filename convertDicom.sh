#!/bin/sh


###############################
#
#		CREATED ON:		2017-10-10
#
#		CREATED BY:		JASON CRAGGS
#
#		USAGE:			CONVERT RAW DCM FILES TO NII FORMAT
#
#		MODIFIED ON:		2017_10_19
#  	MODIFIED ON:	  2017_12_06
#
###############################
#
#
#dicomData=/Volumes/Data/Imaging/Pilot_Saline/Craggs_pilot_saline/1.3.12.2.1107.5.2.32.35338.30000017100413592759300000001/
#niixOutput=/Volumes/Data/Imaging/ConnTest2/Pilot_Saline/

#       AS OF December 6, 2017 - THIS SCRIPT HAS TO BE UPDATED FOR EACH INDIVIDUAL
#       DATASET THAT IS CONVERTED. THIS WILL HAVE TO BE UPDATED FOR USE ON THE
#       CLUSTER 


sub="Sub031"
visit="_v1"
cond="_pos"

dicomData=/Volumes/Data/Imaging/R01/rawdata/${sub}${visit}"/"
niixOutput=/Volumes/Data/Imaging/R01/preprocessed/_Jason_0/${sub}${visit}${cond}"/"

dcm2niix -b y -z y -o ${niixOutput} -f ${sub}_%p_s%2s ${dicomData}


echo ${dicomData}
echo ${niixOutput}
cd ${niixOutput}

#pwd
#			RENAME OUTPUT FILES
# mv ${sub}"_t1_"*".json"  "T1"${visit}".json"
# mv ${sub}"_resting_bold1_"*".json"  "RSrun1"${visit}".json"
# mv ${sub}"_resting_bold2_"*".json"  "RSrun2"${visit}".json"
# mv ${sub}"_resting_bold3_"*".json"  "RSrun3"${visit}".json"
# mv ${sub}"_resting_bold4_"*".json"  "RSrun4"${visit}".json"


mv ${sub}"_t1_"*".gz"  "T1"${visit}".gz"
mv ${sub}"_resting_bold1_"*".gz"  "RSrun1"${visit}".gz"
mv ${sub}"_resting_bold2_"*".gz"  "RSrun2"${visit}".gz"
mv ${sub}"_resting_bold3_"*".gz"  "RSrun3"${visit}".gz"
mv ${sub}"_resting_bold4_"*".gz"  "RSrun4"${visit}".gz"


#		CREATE NEW DIRECTORY AND MOVE EXTRA FILES
mkdir AdditionalScans
mv ${sub}* ./AdditionalScans



#mv ${niixOutput}"/"${sub}"_resting_bold2_"*".json"  "RSrun1"${visit}".json"
#RSrun1_v1


#			EXAMPLE FROM MRIroGL
#dcm2niix" -b y -z y -f "%p"  "MyDicomFolder"
#			EXAMPLE FROM PAUL WRIGHT
#dcm2niix -m y -z y -o ${niftidirectory} -f SubID001_%t_s%s_%d ${dicomdirectory}
