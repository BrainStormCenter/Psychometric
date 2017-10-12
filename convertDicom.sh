#!/bin/sh


###############################
#
#		CREATED ON:		2017-10-10
#
#		CREATED BY:		JASON CRAGGS
#
#		USAGE:			CONVERT RAW DCM FILES TO NII FORMAT
#
#
###############################
#
#
dicomData=/Volumes/Data/Imaging/Pilot_Saline/Craggs_pilot_saline/1.3.12.2.1107.5.2.32.35338.30000017100413592759300000001/
niixOutput=/Volumes/Data/Imaging/ConnTest2/Pilot_Saline/

dcm2niix -b y -z y -o ${niixOutput} -f Sub001_%p_s%2s ${dicomData}




#			EXAMPLE FROM MRIroGL
#dcm2niix" -b y -z y -f "%p"  "MyDicomFolder"
#			EXAMPLE FROM PAUL WRIGHT
#dcm2niix -m y -z y -o ${niftidirectory} -f SubID001_%t_s%s_%d ${dicomdirectory}
