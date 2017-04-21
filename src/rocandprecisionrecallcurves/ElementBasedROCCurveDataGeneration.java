/**
 * 
 */
package rocandprecisionrecallcurves;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import common.Commons;
import auxiliary.FileOperations;
import enumtypes.AssociationMeasureType;
import enumtypes.DataDrivenExperimentCellLineType;
import enumtypes.DataDrivenExperimentDnaseOverlapExclusionType;
import enumtypes.DataDrivenExperimentElementNameType;
import enumtypes.DataDrivenExperimentGeneType;
import enumtypes.DataDrivenExperimentTPMType;
import enumtypes.GenerateRandomDataMode;
import enumtypes.IsochoreFamilyMode;
import enumtypes.ToolType;
import gnu.trove.map.TObjectFloatMap;
import gnu.trove.map.hash.TObjectFloatHashMap;

/**
 * @author Burçak Otlu
 * @date Mar 24, 2017
 * @project GLANETBioinformatics 
 *
 */
public class ElementBasedROCCurveDataGeneration {
	
	public static void fillWithDataDrivenExperimentElementNameTypesIncludingAmbigiousElements(
			DataDrivenExperimentCellLineType cellLineType,
			List<DataDrivenExperimentElementNameType> activatorElementList,
			List<DataDrivenExperimentElementNameType> repressorElementList){
		
		
		//Activators
		activatorElementList.add(DataDrivenExperimentElementNameType.POL2);
		activatorElementList.add(DataDrivenExperimentElementNameType.H2AZ);
		activatorElementList.add(DataDrivenExperimentElementNameType.H3K27AC);
		activatorElementList.add(DataDrivenExperimentElementNameType.H3K4ME2);
		activatorElementList.add(DataDrivenExperimentElementNameType.H3K4ME3);		
		activatorElementList.add(DataDrivenExperimentElementNameType.H3K79ME2);		
		activatorElementList.add(DataDrivenExperimentElementNameType.H3K9AC);		
		activatorElementList.add(DataDrivenExperimentElementNameType.H3K4ME1); 		// <--- Ambigious	
		
		//Let's remove this ambigious
		activatorElementList.add(DataDrivenExperimentElementNameType.H4K20ME1); 	// <--- Ambigious
		
		//We have to consider H3K36ME3 otherwise AUC increases enormously
	    activatorElementList.add(DataDrivenExperimentElementNameType.H3K36ME3); 	// <--- Ambigious		
						
		//Repressors 
		repressorElementList.add(DataDrivenExperimentElementNameType.H3K27ME3);		
		
		//We have to consider H3K9me3 otherwise woGCM has the highest AUC in Expressed scenario
		repressorElementList.add(DataDrivenExperimentElementNameType.H3K9ME3); 		// <--- Ambigious		
		
		
		switch (cellLineType){
		
			case K562:
				//Activators
				activatorElementList.add(DataDrivenExperimentElementNameType.H3K9ACB);		
				break;
				
			default:
				break;
		
		}//End  of switch
							
	}
	
	public static String createTag(
			ToolType toolType,
			DataDrivenExperimentCellLineType cellLine,
			DataDrivenExperimentGeneType geneType,
			DataDrivenExperimentTPMType tpmType,
			DataDrivenExperimentDnaseOverlapExclusionType exclusionType,
			AssociationMeasureType measureType,
			IsochoreFamilyMode isochoreFamilyMode,
			GenerateRandomDataMode generateRandomDataMode){
		
		String tag = "";
		
		
		switch(toolType){
		
			case GLANET:
				tag = 	cellLine.convertEnumtoString() + "_" +
						measureType.convertEnumtoShortString() + "_" + 
						isochoreFamilyMode.convertEnumtoShortString() + "_" + 
						generateRandomDataMode.convertEnumtoShortString();
				break;
				
			case GAT:
				
				if (generateRandomDataMode.isGenerateRandomDataModeWithGC()){
					
					tag = 	cellLine.convertEnumtoString() + "_" +
							measureType.convertEnumtoShortString() + "_" + 
							IsochoreFamilyMode.DO_USE_ISOCHORE_FAMILY.convertEnumtoShortString();
				
				}else if (generateRandomDataMode.isGenerateRandomDataModeWithoutMapabilityandGc()){
					
					tag = 	cellLine.convertEnumtoString() + "_" +
							measureType.convertEnumtoShortString() + "_" + 
							IsochoreFamilyMode.DO_NOT_USE_ISOCHORE_FAMILY.convertEnumtoShortString();
				}
				
				tag = tag + "_" + toolType.convertEnumtoString();
				break;
				
			default:
				break;
				
		}//End of switch
				
		return tag;
		
	}
	
	
	
	public static Float getEmpiricalPValue(String strLine){
		
		int indexofFirstTab;
		int indexofSecondTab;
		int indexofThirdTab;
		int indexofFourthTab;
		int indexofFifthTab;
		int indexofSixthTab;
		int indexofSeventhTab;
		int indexofEigthTab;
		int indexofNinethTab;
		int indexofTenthTab;
		
		indexofFirstTab 	= strLine.indexOf('\t');
		indexofSecondTab 	= (indexofFirstTab > 0)?strLine.indexOf('\t', indexofFirstTab + 1):-1;
		indexofThirdTab 	= (indexofSecondTab > 0)?strLine.indexOf('\t', indexofSecondTab + 1):-1;
		indexofFourthTab 	= (indexofThirdTab > 0)?strLine.indexOf('\t', indexofThirdTab + 1):-1;
		indexofFifthTab 	= (indexofFourthTab > 0)?strLine.indexOf('\t', indexofFourthTab + 1):-1;
		indexofSixthTab 	= (indexofFifthTab > 0)?strLine.indexOf('\t', indexofFifthTab + 1):-1;
		indexofSeventhTab 	= (indexofSixthTab > 0)?strLine.indexOf('\t', indexofSixthTab + 1):-1;
		indexofEigthTab 	= (indexofSeventhTab > 0)?strLine.indexOf('\t', indexofSeventhTab + 1):-1;
		indexofNinethTab 	= (indexofEigthTab > 0)?strLine.indexOf('\t', indexofEigthTab + 1):-1;
		indexofTenthTab 	= (indexofNinethTab > 0)?strLine.indexOf('\t', indexofNinethTab + 1):-1;
		
		//Between indexofNinethTab and indexofTenthTab
		Float empiricalPValue = Float.parseFloat(strLine.substring(indexofNinethTab+1, indexofTenthTab));
		
		return empiricalPValue;
	}
	
	
	//For each output file 
	public static void readResultFile(
			String resultFileName, 
			List<DataDrivenExperimentElementNameType> activatorElementList,
			List<DataDrivenExperimentElementNameType> repressorElementList,
			TObjectFloatMap<DataDrivenExperimentElementNameType> activator2EmpPValueMap,
			TObjectFloatMap<DataDrivenExperimentElementNameType> repressor2EmpPValueMap){
		
		
		FileReader fileReader = null;
		BufferedReader bufferedReader = null;
		String strLine = null;
		
		int indexofFirstUnderscore;
		String elementName = null;
		Float empiricalPValue = null;
		
		DataDrivenExperimentElementNameType element = null;
		
		
		try {
			fileReader = FileOperations.createFileReader(resultFileName);
			bufferedReader = new BufferedReader(fileReader);
					
			
			while((strLine = bufferedReader.readLine())!=null){
				
				if (!strLine.startsWith("#")){
					
					indexofFirstUnderscore = strLine.indexOf('_');
					elementName= strLine.substring(0,indexofFirstUnderscore);
					
					element = DataDrivenExperimentElementNameType.convertStringtoEnum(elementName);
					
					if (activatorElementList.contains(element)){
						
						//gets its empirical p value
						empiricalPValue = getEmpiricalPValue(strLine);
						activator2EmpPValueMap.put(element,empiricalPValue);
						
												
					}else if (repressorElementList.contains(element)){
						
						empiricalPValue = getEmpiricalPValue(strLine);
						repressor2EmpPValueMap.put(element,empiricalPValue);
						
					}				
					
				}//End of IF not an comment line
				
			}//End of while
			
			//Close
			bufferedReader.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}	
		
	}
	
	
	//We use empirical pValue
	public static void prepareTabDelimitedDataFileForGLANET(
			String mainDirectory, 
			String runDirectory, 
			int numberofRuns, 
			DataDrivenExperimentGeneType geneType, 
			String tag, 
			List<DataDrivenExperimentElementNameType> activatorElementList, 
			List<DataDrivenExperimentElementNameType> repressorElementList,
			BufferedWriter elementSpecificBufferedWriterROCCurve){
		
		String eachRunDirectory = null;
		File histoneDirectory = null;
		File tfDirectory = null;
		
		String resultFileName = null;	
		
		boolean thereIsSuchAFile = false;
		
		int numberofFilesReadforHistone = 0;
		int numberofFilesReadforTF = 0;
			
		DataDrivenExperimentElementNameType element = null;
		float pValue = 0f;
		
		//For each run
		for(int i=0; i< numberofRuns; i++){
			
			//Create Instances
			TObjectFloatMap<DataDrivenExperimentElementNameType> activator2EmpPValueMap  =  new TObjectFloatHashMap<DataDrivenExperimentElementNameType>();
			TObjectFloatMap<DataDrivenExperimentElementNameType> repressor2EmpPValueMap  =  new TObjectFloatHashMap<DataDrivenExperimentElementNameType>();
					
			eachRunDirectory = mainDirectory + runDirectory + i + System.getProperty("file.separator") + Commons.ENRICHMENT + System.getProperty("file.separator");
			
			histoneDirectory  = new File(eachRunDirectory + Commons.HISTONE + System.getProperty("file.separator"));
			tfDirectory  = new File(eachRunDirectory + Commons.TF + System.getProperty("file.separator"));
			
			//Histone
			if(histoneDirectory.exists() && histoneDirectory.isDirectory()){
				
				//Initialize
				thereIsSuchAFile = false;
				
				for(File eachHistoneFile : histoneDirectory.listFiles()){
					
					//Read cellLineFiltered file
					if(!eachHistoneFile.isDirectory() && 
							(eachHistoneFile.getAbsolutePath().contains("_2017_"))){
						
						//formerly it was  contains("Run" + i + ".txt")
						
						thereIsSuchAFile = true;
						numberofFilesReadforHistone++;

						resultFileName = eachHistoneFile.getAbsolutePath();								
						readResultFile(resultFileName,activatorElementList,repressorElementList,activator2EmpPValueMap,repressor2EmpPValueMap);						
						break;

					}// End of IF EnrichmentFile under EnrichmentDirectory
					
				}//End of for each file in histone Directory
				
				//Check
				if(!thereIsSuchAFile){
					System.out.println("There is no result file for directory " + histoneDirectory);
				}
				
				
			}//End of IF histone directory exists
			
			
			//TF
			if(tfDirectory.exists() && tfDirectory.isDirectory()){
				
				//initialize
				thereIsSuchAFile = false;
				
				for(File eachTFFile : tfDirectory.listFiles()){
					
					//read cellLineFiltered file
					if(!eachTFFile.isDirectory() && 
							(eachTFFile.getAbsolutePath().contains("_2017_"))){
						
						//formerly it was  contains("Run" + i + ".txt")
						
						thereIsSuchAFile = true;
						numberofFilesReadforTF++;

						resultFileName = eachTFFile.getAbsolutePath();								
						readResultFile(resultFileName,activatorElementList,repressorElementList,activator2EmpPValueMap,repressor2EmpPValueMap);						
						break;

					}// End of IF EnrichmentFile under EnrichmentDirectory
					
				}//End of for each file in histone Directory
				
				//Check
				if(!thereIsSuchAFile){
					System.out.println("There is no result file for directory " + tfDirectory);
				}
				
			}//End of IF TF directory exists
			
			//For all missing activator elements put a pValue of 1.
			for(Iterator<DataDrivenExperimentElementNameType> itr =activatorElementList.iterator();itr.hasNext();){
				
				element = itr.next();
				
				if (!activator2EmpPValueMap.containsKey(element)){
					activator2EmpPValueMap.put(element,1.0f);					
				}
				
			}//End of FOR
			
			
			//For all missing repressor elements put a pValue of 1.
			for(Iterator<DataDrivenExperimentElementNameType> itr =repressorElementList.iterator();itr.hasNext();){
				
				element = itr.next();
				
				if (!repressor2EmpPValueMap.containsKey(element)){
					repressor2EmpPValueMap.put(element,1.0f);					
				}
				
			}//End of FOR
			
			
			try {

				//write to a file
				//elementName tab empiricalNValue
				//with elementName in ascending order
				//try in laptop then move the jar to levreks
				
				//Write in a sorted way
				for(int j=0; j<activatorElementList.size();j++){
					
					element = activatorElementList.get(j);
					
					pValue = activator2EmpPValueMap.get(element);
					
					if (geneType.isExpressingProteinCodingGenes()){
						elementSpecificBufferedWriterROCCurve.write("0" + "\t" + "Enriched" + "\t" + element.convertEnumtoString() + "\t" + pValue + "\t" + tag + System.getProperty("line.separator"));						
												
					}else {
						elementSpecificBufferedWriterROCCurve.write("1" + "\t" + "NotEnriched" + "\t" + element.convertEnumtoString() + "\t" + pValue + "\t" + tag + System.getProperty("line.separator"));	
					}
										
				}
				
				for(int j=0; j<repressorElementList.size();j++){
					
					element = repressorElementList.get(j);
					
					pValue = repressor2EmpPValueMap.get(element);
					
					if (geneType.isExpressingProteinCodingGenes()){
						elementSpecificBufferedWriterROCCurve.write("1" + "\t" + "NotEnriched" + "\t" + element.convertEnumtoString() + "\t" +  pValue + "\t" + tag +  System.getProperty("line.separator"));		
					}else {
						elementSpecificBufferedWriterROCCurve.write("0" + "\t" + "Enriched" + "\t" + element.convertEnumtoString() + "\t" + pValue + "\t" + tag + System.getProperty("line.separator"));	
					}
				}
				
				
			
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			//Free space
			activator2EmpPValueMap = null;
			repressor2EmpPValueMap = null;

		}//End of For each run number
		
		//Check
		if (numberofFilesReadforHistone!=1000){
			System.out.println("numberofFilesReadforHistone: " + numberofFilesReadforHistone);
		}
		
		//Check
		if (numberofFilesReadforTF!=1000){
			System.out.println("numberofFilesReadforTF: " + numberofFilesReadforTF);
		}

		
	}
	
	
	//For collecting GAT ROC Curve Data
	public static void processLine(
			String strLine,
			String tag,
			DataDrivenExperimentGeneType geneType,
			List<DataDrivenExperimentElementNameType> activatorElementList, 
			List<DataDrivenExperimentElementNameType> repressorElementList,
			BufferedWriter bufferedWriter) throws IOException{
		
		String cellLineNameElementName = null;
		String elementName = null;
		DataDrivenExperimentElementNameType elementNameType = null;
	
		float empiricalPValue;
		float ln2Fold;
		
		int indexofUnderscore;

		int indexofFirstTab;
		int indexofSecondTab;
		int indexofThirdTab;
		int indexofFourthTab;
		int indexofFifthTab;
		int indexofSixthTab;
		int indexofSeventhTab;
		int indexofEigthTab;
		int indexofNinethTab;
		int indexofTenthTab;
		
		indexofFirstTab = strLine.indexOf( '\t');
		indexofSecondTab = ( indexofFirstTab > 0)?strLine.indexOf( '\t', indexofFirstTab + 1):-1;
		indexofThirdTab = ( indexofSecondTab > 0)?strLine.indexOf( '\t', indexofSecondTab + 1):-1;
		indexofFourthTab = ( indexofThirdTab > 0)?strLine.indexOf( '\t', indexofThirdTab + 1):-1;
		indexofFifthTab = ( indexofFourthTab > 0)?strLine.indexOf( '\t', indexofFourthTab + 1):-1;
		indexofSixthTab = ( indexofFifthTab > 0)?strLine.indexOf( '\t', indexofFifthTab + 1):-1;
		indexofSeventhTab = ( indexofSixthTab > 0)?strLine.indexOf( '\t', indexofSixthTab + 1):-1;
		indexofEigthTab = ( indexofSeventhTab > 0)?strLine.indexOf( '\t', indexofSeventhTab + 1):-1;
		indexofNinethTab = ( indexofEigthTab > 0)?strLine.indexOf( '\t', indexofEigthTab + 1):-1;
		indexofTenthTab = ( indexofNinethTab > 0)?strLine.indexOf( '\t', indexofNinethTab + 1):-1;
		
		if(indexofFirstTab!=-1 && indexofSecondTab!=-1 && indexofNinethTab!=-1 && indexofTenthTab!=-1){
			
			cellLineNameElementName = strLine.substring(indexofFirstTab + 1, indexofSecondTab);
			indexofUnderscore = cellLineNameElementName.indexOf("_");
			elementName = cellLineNameElementName.substring(indexofUnderscore+1);
			
			elementNameType = DataDrivenExperimentElementNameType.convertStringtoEnum(elementName);			
				
			//In case of enrichment ln2Fold must be positive and pValue must be less than significance level
			//In case of depletion ln2Fold must be negative and pValue must be less than significance level
					
			//ln2fold is between 8thTab and 9thTab
			ln2Fold = Float.parseFloat(strLine.substring(indexofEigthTab + 1, indexofNinethTab));
			
			//pValue is between 9thTab and 10thTab
			//pValue is both used for enrichment and depletion
			empiricalPValue = Float.parseFloat(strLine.substring(indexofNinethTab + 1, indexofTenthTab));
			
			
			//TODO Shall I consider ln2Fold?
			//Decision: Yes
			//If ln2fold is positive take empiricalPValue as it is
			//If ln2fold is negative take 1-empiricalPValue
			if(ln2Fold<0){
				empiricalPValue = 1-empiricalPValue;
			}
			
			//Update and Accumulate
			if (activatorElementList.contains(elementNameType)){
				
				switch (geneType) {
				
					case  EXPRESSING_PROTEINCODING_GENES:
						bufferedWriter.write("0" + "\t" + "Enriched" + "\t" +  elementNameType.convertEnumtoString() + "\t"  + empiricalPValue + "\t" + tag + System.getProperty("line.separator"));																				
						break;
					
						
					case NONEXPRESSING_PROTEINCODING_GENES:
						bufferedWriter.write("1" + "\t" + "NotEnriched" + "\t" + elementNameType.convertEnumtoString() + "\t"  + empiricalPValue + "\t" + tag + System.getProperty("line.separator"));													
						break;
	
					default:
						break;
				}//End of switch
				
			}else if (repressorElementList.contains(elementNameType)){
				
				switch (geneType) {
				
					case  EXPRESSING_PROTEINCODING_GENES:					
						bufferedWriter.write("1" + "\t" + "NotEnriched" + "\t" + elementNameType.convertEnumtoString() + "\t"  + empiricalPValue + "\t" + tag + System.getProperty("line.separator"));						
						break;
						
					case NONEXPRESSING_PROTEINCODING_GENES:
						bufferedWriter.write("0" + "\t" + "Enriched" + "\t" + elementNameType.convertEnumtoString() + "\t"  + empiricalPValue + "\t" + tag + System.getProperty("line.separator"));						
						break;
	
					default:
						break;
				}//End of switch
				
			}
				
			
		}//End of IF valid strLine control
		
	}

	
	public static void prepareTabDelimitedDataFileForGAT(
			String mainDirectory, 
			String runDirectory, 
			int numberofRuns, 
			DataDrivenExperimentGeneType geneType, 
			String tag, 
			List<DataDrivenExperimentElementNameType> activatorElementList, 
			List<DataDrivenExperimentElementNameType> repressorElementList,
			BufferedWriter elementSpecificBufferedWriterROCCurve){
		
		String strLine = null;
		File gatTSVFile = null;
		
		FileReader gatTSVFileReader = null;
		BufferedReader gatTSVBufferedReader = null;
			
		try{

			// For each run
			for(int i = 0; i <numberofRuns; i++){
								
				gatTSVFile = new File(mainDirectory + System.getProperty("file.separator")  + runDirectory + i + Commons.TSV);
									
				if (gatTSVFile.exists() && gatTSVFile.isFile()){
					
					//check it
					gatTSVFileReader = FileOperations.createFileReader(gatTSVFile.getAbsolutePath());
					gatTSVBufferedReader = new BufferedReader(gatTSVFileReader);
					
					//example lines
//					track	annotation	observed	expected	CI95low	CI95high	stddev	fold	l2fold	pvalue	qvalue	track_nsegments	track_size	track_density	annotation_nsegments	annotation_size	annotation_density	overlap_nsegments	overlap_size	overlap_density	percent_overlap_nsegments_track	percent_overlap_size_track	percent_overlap_nsegments_annotation	percent_overlap_size_annotation
//					merged	GM12878_H3K4ME3	5042	10527.8922	6860	14481	2316.527	0.479	-1.062	6.10E-03	1.22E-02	500	276698	8.94E-03	27920	101790620	3.29E+00	24	5042	1.63E-04	4.8	1.8222	0.086	0.005
//					merged	GM12878_H3K4ME1	8039	15521.7919	11019	20190	2779.2138	0.5179	-0.9491	1.40E-03	5.60E-03	500	276698	8.94E-03	41605	156124633	5.04E+00	24	8039	2.60E-04	4.8	2.9053	0.0577	0.0051
//					merged	GM12878_H3K4ME2	9057	13720.137	9604	18156	2602.1074	0.6601	-0.5991	2.96E-02	3.95E-02	500	276698	8.94E-03	45162	134962712	4.36E+00	35	9057	2.93E-04	7	3.2732	0.0775	0.0067
//					merged	GM12878_H2AZ	6764	8783.6121	5543	12379	2094.3924	0.7701	-0.3769	1.71E-01	1.71E-01	500	276698	8.94E-03	40693	89430210	2.89E+00	27	6764	2.19E-04	5.4	2.4445	0.0664	0.0076

					// Skip HeaderLine
					strLine = gatTSVBufferedReader.readLine();
										
					// Read gatTSVFile
					while( ( strLine = gatTSVBufferedReader.readLine()) != null){
						
						processLine(strLine,tag,geneType,activatorElementList, repressorElementList, elementSpecificBufferedWriterROCCurve);
						
					}// End of WHILE
					
					
					// Close
					gatTSVBufferedReader.close();
					
				}//End of IF gatTSVFile exists
				
				else{
					//Write down this file does not exists
					//Exit from loop
					System.out.println("There is no file for this run: " + i + "\t" +  gatTSVFile.getAbsolutePath());
				
				}
						
			}// End of FOR each run
			
			
		}catch( IOException e){
			e.printStackTrace();
		}
		
	}

	
	public static void prepareElementSpecificROCCurveData(
			String mainDirectory,
			int numberofRuns,
			ToolType toolType,
			DataDrivenExperimentCellLineType cellLine,
			DataDrivenExperimentDnaseOverlapExclusionType exclusionType,
			DataDrivenExperimentTPMType tpmType){
		
		FileWriter elementSpecificFileWriterROCCurve = null;
		BufferedWriter elementSpecificBufferedWriterROCCurve = null;
		
		List<DataDrivenExperimentElementNameType> activatorElementList = new ArrayList<DataDrivenExperimentElementNameType>();
		List<DataDrivenExperimentElementNameType> repressorElementList = new ArrayList<DataDrivenExperimentElementNameType>();
		
		List<DataDrivenExperimentGeneType> geneTypeList= new ArrayList<DataDrivenExperimentGeneType>();		
		DataDrivenExperimentGeneType geneType = null;
		
		List<DataDrivenExperimentDnaseOverlapExclusionType> NonExp_DnaseExclusionTypeList= new ArrayList<DataDrivenExperimentDnaseOverlapExclusionType>();		
		List<DataDrivenExperimentTPMType> Exp_tpmTypeList= new ArrayList<DataDrivenExperimentTPMType>();	
		
		List<AssociationMeasureType> associationMeasureTypeList= new ArrayList<AssociationMeasureType>();				
		List<IsochoreFamilyMode> isochoreFamilyModeList= new ArrayList<IsochoreFamilyMode>();
		List<GenerateRandomDataMode> generateRandomDataModeList= new ArrayList<GenerateRandomDataMode>();
		
		AssociationMeasureType measureType = null;
		IsochoreFamilyMode isochoreFamilyMode = null;
		GenerateRandomDataMode generateRandomDataMode = null;
		
		String runDirectory = null;
		
		String tag = "";
		String distinguisingFileName = 	cellLine.convertEnumtoString() + "_" + 
										exclusionType.convertEnumtoString() + "_" + 
										tpmType.convertEnumtoString();
				
		
		
		//Fill starts
		geneTypeList.add(DataDrivenExperimentGeneType.NONEXPRESSING_PROTEINCODING_GENES);
		geneTypeList.add(DataDrivenExperimentGeneType.EXPRESSING_PROTEINCODING_GENES);
		
		NonExp_DnaseExclusionTypeList.add(exclusionType);
		Exp_tpmTypeList.add(tpmType);
		
		associationMeasureTypeList.add(AssociationMeasureType.EXISTENCE_OF_OVERLAP);
		associationMeasureTypeList.add(AssociationMeasureType.NUMBER_OF_OVERLAPPING_BASES);
		
		isochoreFamilyModeList.add(IsochoreFamilyMode.DO_USE_ISOCHORE_FAMILY);
		isochoreFamilyModeList.add(IsochoreFamilyMode.DO_NOT_USE_ISOCHORE_FAMILY);
		
		generateRandomDataModeList.add(GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITH_GC_CONTENT);
		generateRandomDataModeList.add(GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITH_MAPPABILITY);
		generateRandomDataModeList.add(GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITH_MAPPABILITY_AND_GC_CONTENT);
		generateRandomDataModeList.add(GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITHOUT_MAPPABILITY_AND_GC_CONTENT);
		//Fill ends
				
		try {
			
			elementSpecificFileWriterROCCurve = FileOperations.createFileWriter(mainDirectory + System.getProperty("file.separator") + toolType.convertEnumtoString() + "_ROC_Curve_Data" + System.getProperty("file.separator") + toolType.convertEnumtoString() + "_ROC_Curve_Data_" +  distinguisingFileName + ".txt");
			elementSpecificBufferedWriterROCCurve = new BufferedWriter(elementSpecificFileWriterROCCurve);
			
			//Header line ROC Curve
			elementSpecificBufferedWriterROCCurve.write("D" + "\t" + "D.str" + "\t" + "Element" + "\t" + "M"  + "\t" + "Parameter" + System.getProperty("line.separator")); 	

			//Initialize
			fillWithDataDrivenExperimentElementNameTypesIncludingAmbigiousElements(cellLine,activatorElementList,repressorElementList);

			for(Iterator<DataDrivenExperimentGeneType> geneTypeItr = geneTypeList.iterator(); geneTypeItr.hasNext();){
				
				geneType = geneTypeItr.next();
				
				if (geneType.isNonExpressingProteinCodingGenes()){
					exclusionType = NonExp_DnaseExclusionTypeList.get(0);							
					tpmType = DataDrivenExperimentTPMType.TOPUNKNOWN;
				}else if (geneType.isExpressingProteinCodingGenes()){						
					exclusionType = DataDrivenExperimentDnaseOverlapExclusionType.NO_DISCARD;							
					tpmType = Exp_tpmTypeList.get(0);
				}
				
				
				for(Iterator<AssociationMeasureType> measureTypeItr = associationMeasureTypeList.iterator();measureTypeItr.hasNext();){
					
					measureType = measureTypeItr.next();
									
					switch(toolType){
					
						case GLANET:
							/*******************************************************************************/
							/*******************************GLANET starts***********************************/
							/*******************************************************************************/
							
							//If you give one specific isochoreFamily there will be only one isochoreFamily in the isochoreFamilyModeList
							for(Iterator<IsochoreFamilyMode> isochoreFamilyModeItr =  isochoreFamilyModeList.iterator();isochoreFamilyModeItr.hasNext();) {
								
								isochoreFamilyMode = isochoreFamilyModeItr.next();
								
								for(Iterator<GenerateRandomDataMode> generateRandomDataModeItr = generateRandomDataModeList.iterator(); generateRandomDataModeItr.hasNext();){
									
									generateRandomDataMode =generateRandomDataModeItr.next();
									
									tag = createTag(toolType,
													cellLine,
													geneType,
													tpmType,
													exclusionType,
													measureType,
													isochoreFamilyMode,
													generateRandomDataMode);
									
									//Here we read isochore specific runDirectory
									runDirectory= "Output" + System.getProperty("file.separator") + cellLine + "_" + geneType.convertEnumtoString() + "_" + tpmType.convertEnumtoString() + "_" + exclusionType.convertEnumtoString() + "_" + generateRandomDataMode.convertEnumtoShortString() + "_" + isochoreFamilyMode.convertEnumtoShortString() + "_" + measureType.convertEnumtoShortString() + "Run";
									
									//We use empirical p-Value
									prepareTabDelimitedDataFileForGLANET(
											mainDirectory, 
											runDirectory, 
											numberofRuns, 
											geneType, 
											tag, 
											activatorElementList, 
											repressorElementList,
											elementSpecificBufferedWriterROCCurve);								
								
								}//End of for each generateRandomDataMode														
								
							}//End of each isochoreFamilyMode
							/*******************************************************************************/
							/*******************************GLANET ends*************************************/
							/*******************************************************************************/									
							break;
						
						case GAT:
							//TODO
							/*******************************************************************************/
							/*********************************GAT starts************************************/
							/*******************************************************************************/									
							//Initialization
							//For GAT I wrote wIF and woIF by myself in createTag method
							//Because GAT output filenames does not has wIF woIF in their names.
							
							//If isochoreFamilyMode is wIF collect only wGC(EOO) and wGCM(NOOB) results
							//If isochoreFamilyMode is woIF collect woGCM results (same for EOO and NOOB)
							//If thereis no isochoreFamilyMode set then collect all the results, simply pool them
							generateRandomDataModeList.clear();
							generateRandomDataModeList.add(GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITH_GC_CONTENT);
							generateRandomDataModeList.add(GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITHOUT_MAPPABILITY_AND_GC_CONTENT);											
							
														
							//GAT in fact only achieves wIF and woIF
							//For wIF I have named the files with wGC (EOO) and wGCM (NOOB)									
							//For woIF I have named the files with woGCM
							
							for(Iterator<GenerateRandomDataMode> generateRandomDataModeItr = generateRandomDataModeList.iterator(); generateRandomDataModeItr.hasNext();){
								
								generateRandomDataMode =generateRandomDataModeItr.next();
								
								tag = createTag(toolType,
										cellLine,
										geneType,
										tpmType,
										exclusionType,
										measureType,
										isochoreFamilyMode,
										generateRandomDataMode);
								
								//GAT DDE Output EOO uses wGC and woGCM in filenames for wGC and woGCM, respectively.
								//GAT DDE Output NOOB uses wGCM and woGCM in filenames for wGC and woGCM, respectively.																				
								if (measureType.isAssociationMeasureNumberOfOverlappingBases() && generateRandomDataMode.isGenerateRandomDataModeWithGC()){
									runDirectory= "Output" + System.getProperty("file.separator") + toolType.convertEnumtoString() + "_" +  cellLine + "_" + geneType.convertEnumtoString() + "_" + tpmType.convertEnumtoString() + "_" + exclusionType.convertEnumtoString() + "_" + GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITH_MAPPABILITY_AND_GC_CONTENT.convertEnumtoShortString() + "_" + measureType.convertEnumtoShortString() + "_" + Commons.DDE_RUN;
									
								}else{
									runDirectory= "Output" + System.getProperty("file.separator") + toolType.convertEnumtoString() + "_" +  cellLine + "_" + geneType.convertEnumtoString() + "_" + tpmType.convertEnumtoString() + "_" + exclusionType.convertEnumtoString() + "_" + generateRandomDataMode.convertEnumtoShortString() + "_" + measureType.convertEnumtoShortString() + "_" + Commons.DDE_RUN;
				
								}
								
								
								//We use empirical p-value
								prepareTabDelimitedDataFileForGAT(
										mainDirectory, 
										runDirectory, 
										numberofRuns, 
										geneType, 
										tag, 
										activatorElementList, 
										repressorElementList,
										elementSpecificBufferedWriterROCCurve);								
							
							}//End of for each generateRandomDataMode											
							/*******************************************************************************/
							/*********************************GAT ends**************************************/
							/*******************************************************************************/

							break;
						
						default:
							break;
					
					}//End of switch
					
				}//End of for each association measure type (EOO, NOOB)
								
			}//End of for each geneType (NonExp, Exp)
								
			//Close
			elementSpecificBufferedWriterROCCurve.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		

		
	}

	
	public static void main(String[] args) {

		//We need 6 arguments
		//main directory (parent of Output Directory)
		//numberofRuns		
		//tooltype (GLANET, GAT)
		//cell (GM12878, K562)
		//exclusionType (CompletelyDiscard,TakeTheLongest)
		//tpmType (Top5,Top20)
		
		String mainDirectory=args[0];				
		int numberofRuns = Integer.parseInt(args[1]);
		ToolType toolType = ToolType.convertStringtoEnum(args[2]);
		DataDrivenExperimentCellLineType cellLine = DataDrivenExperimentCellLineType.convertStringtoEnum(args[3]);		
		DataDrivenExperimentDnaseOverlapExclusionType exclusionType = DataDrivenExperimentDnaseOverlapExclusionType.convertStringtoEnum(args[4]);		
		DataDrivenExperimentTPMType tpmType = DataDrivenExperimentTPMType.convertStringtoEnum(args[5]);
		

		//We output element specific pValues
		prepareElementSpecificROCCurveData(
				mainDirectory,
				numberofRuns,
				toolType,
				cellLine,
				exclusionType,
				tpmType);
		
	}

}
