/**
 * 
 */
package typeIerrorandpower;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import auxiliary.FileOperations;
import common.Commons;
import enumtypes.AssociationMeasureType;
import enumtypes.DataDrivenExperimentCellLineType;
import enumtypes.DataDrivenExperimentDnaseOverlapExclusionType;
import enumtypes.DataDrivenExperimentElementNameType;
import enumtypes.DataDrivenExperimentGeneType;
import enumtypes.DataDrivenExperimentTPMType;
import enumtypes.GenerateRandomDataMode;
import enumtypes.IsochoreFamilyMode;
import enumtypes.ToolType;
import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.list.TFloatList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.TFloatObjectMap;
import gnu.trove.map.TObjectFloatMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TFloatObjectHashMap;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;

/**
 * @author Burçak Otlu
 * @date Mar 14, 2017
 * @project GLANETBioinformatics 
 * 
 * This class prepares TypeI errors and power for various significance levels.
 * 
 * This class prepares these data for GLANET and GAT.
 * 
 * This class will read the runs and collect number of enrichments per element for each significance level of alpha.  
 * For non-expressed genes scenario, it will provide typeI errors for activators and power for repressors. 
 * For expressed genes scenario, it will provide typeI errors for repressors and power for activators. 
 * 
 * It writes tags w.r.t (EOO,NOOB) (woIF,wIF) (wGC,wM,wGCM,woGCM)
 *
 */
public class DataPreparationForTypeIErrorAndPowerFigures {
	
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
		activatorElementList.add(DataDrivenExperimentElementNameType.H4K20ME1); 	// <--- Ambigious
		activatorElementList.add(DataDrivenExperimentElementNameType.H3K36ME3); 	// <--- Ambigious		
						
		//Repressors 
		repressorElementList.add(DataDrivenExperimentElementNameType.H3K27ME3);				
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

	public static void fillMaps(
			TFloatList significanceLevelList,
			TFloatObjectMap<TObjectIntMap<String>> significanceLevel2ElementNameTag2NumberofEnrichmentMapMap){
		
		Float significanceLevel= null;
		
		for(int i=0; i<=25; i++){
			significanceLevelList.add(0.01f*i);
		}
		
//		Float significanceLevel = 0.0f;
//		significanceLevelList.add(0.001f);
//		significanceLevelList.add(0.005f);
//		significanceLevelList.add(0.01f);
//		significanceLevelList.add(0.05f);
		
		
		//Initialize
		for(int i=0; i< significanceLevelList.size(); i++){
			
			significanceLevel = significanceLevelList.get(i); 			
			significanceLevel2ElementNameTag2NumberofEnrichmentMapMap.put(significanceLevel,new TObjectIntHashMap<String>());
			
		}//End of for each significance level
				
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
		indexofSecondTab 	= (indexofFirstTab>0)?strLine.indexOf('\t', indexofFirstTab + 1):-1;
		indexofThirdTab 	= (indexofSecondTab>0)?strLine.indexOf('\t', indexofSecondTab + 1):-1;
		indexofFourthTab 	= (indexofThirdTab>0)?strLine.indexOf('\t', indexofThirdTab + 1):-1;
		indexofFifthTab 	= (indexofFourthTab>0)?strLine.indexOf('\t', indexofFourthTab + 1):-1;
		indexofSixthTab 	= (indexofFifthTab>0)?strLine.indexOf('\t', indexofFifthTab + 1):-1;
		indexofSeventhTab 	= (indexofSixthTab>0)?strLine.indexOf('\t', indexofSixthTab + 1):-1;
		indexofEigthTab 	= (indexofSeventhTab>0)?strLine.indexOf('\t', indexofSeventhTab + 1):-1;
		indexofNinethTab 	= (indexofEigthTab>0)?strLine.indexOf('\t', indexofEigthTab + 1):-1;
		indexofTenthTab 	= (indexofNinethTab>0)?strLine.indexOf('\t', indexofNinethTab + 1):-1;
		
		//Between indexofNinethTab and indexofTenthTab
		Float empiricalPValue = Float.parseFloat(strLine.substring(indexofNinethTab+1, indexofTenthTab));
		
		return empiricalPValue;
	}
	
	
	public static void readResultFileAndFillMap(
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
	
	public static void prepareTabDelimitedDataFileForGLANETTypeIErrorAndPowerFigures(
			String mainDirectory, 
			String runDirectory, 
			int numberofRuns, 
			DataDrivenExperimentDnaseOverlapExclusionType exclusionType,
			DataDrivenExperimentTPMType tpmType,
			AssociationMeasureType measureType,
			IsochoreFamilyMode isochoreFamilyMode,
			GenerateRandomDataMode generateRandomDataMode,
			TFloatList significanceLevelList,
			List<DataDrivenExperimentElementNameType> activatorElementList, 
			List<DataDrivenExperimentElementNameType> repressorElementList,
			TObjectIntMap<String> elementNameTag2NumberofTotalRunsMap,
			TFloatObjectMap<TObjectIntMap<String>> significanceLevel2ElementNameTag2NumberofEnrichmentMapMap){
		
		String eachRunDirectory = null;
		File histoneDirectory = null;
		File tfDirectory = null;
		
		String resultFileName = null;	
		
		boolean thereIsSuchAFile = false;
		
		int numberofFilesReadforHistone = 0;
		int numberofFilesReadforTF = 0;
			
		DataDrivenExperimentElementNameType elementNameType = null;
		float pValue = 0f;
		
		Float significanceLevel = null;
		TObjectIntMap<String> elementNameTag2NumberofEnrichmentMap;
		
		String tag = exclusionType.convertEnumtoString() +  "_"  + tpmType.convertEnumtoString() + "_" + measureType.convertEnumtoShortString() + "_" + isochoreFamilyMode.convertEnumtoShortString() + "_" + generateRandomDataMode.convertEnumtoShortString();
		String elementNameTag = null;
		
		
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
						
						readResultFileAndFillMap(resultFileName,activatorElementList,repressorElementList,activator2EmpPValueMap,repressor2EmpPValueMap);
						
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
						
						readResultFileAndFillMap(resultFileName,activatorElementList,repressorElementList,activator2EmpPValueMap,repressor2EmpPValueMap);
						
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
				
				elementNameType = itr.next();
				
				if (!activator2EmpPValueMap.containsKey(elementNameType)){
					activator2EmpPValueMap.put(elementNameType,1.0f);					
				}
				
			}//End of FOR
			
			
			//For all missing repressor elements put a pValue of 1.
			for(Iterator<DataDrivenExperimentElementNameType> itr =repressorElementList.iterator();itr.hasNext();){
				
				elementNameType = itr.next();
				
				if (!repressor2EmpPValueMap.containsKey(elementNameType)){
					repressor2EmpPValueMap.put(elementNameType,1.0f);
				}
				
			}//End of FOR
			
			//Now accumulate the results of one unique run
			//Fill significanceLevel2ElementName2NumberofEnrichmentMapMap
			for(int j=0; j<activatorElementList.size();j++){
				
				elementNameType = activatorElementList.get(j);	
				elementNameTag = elementNameType.convertEnumtoString() + "_" + tag;
				pValue = activator2EmpPValueMap.get(elementNameType);
				
				//For each value of significanceLevel
				for(int k=0; k<significanceLevelList.size();k++){
					
					significanceLevel = significanceLevelList.get(k);
					
					elementNameTag2NumberofEnrichmentMap = significanceLevel2ElementNameTag2NumberofEnrichmentMapMap.get(significanceLevel);
					
					//You have to initialize
					if (!elementNameTag2NumberofEnrichmentMap.containsKey(elementNameTag)){
						elementNameTag2NumberofEnrichmentMap.put(elementNameTag, 0);
					}
					
					if (pValue <= significanceLevel){							
						if (elementNameTag2NumberofEnrichmentMap.get(elementNameTag)==0){
							elementNameTag2NumberofEnrichmentMap.put(elementNameTag, 1);
						}else{
							elementNameTag2NumberofEnrichmentMap.put(elementNameTag, elementNameTag2NumberofEnrichmentMap.get(elementNameTag)+1);
						}
					}//End of if pValue is less than significanceLevel
					
					//There is no enrichment
					//Nothing to do
					
				}//End of for each significanceLevel
				
				//update elementName2NumberofTotalRunsMap
				elementNameTag2NumberofTotalRunsMap.put(elementNameTag,elementNameTag2NumberofTotalRunsMap.get(elementNameTag)+1);
				
									
			}//End of for each activator element
			
			for(int j=0; j<repressorElementList.size();j++){
				
				elementNameType = repressorElementList.get(j);					
				elementNameTag = elementNameType.convertEnumtoString() + "_" + tag;
				pValue = repressor2EmpPValueMap.get(elementNameType);
				
				//For each value of significanceLevel
				for(int k=0; k<significanceLevelList.size();k++){
					
					significanceLevel = significanceLevelList.get(k);
					
					elementNameTag2NumberofEnrichmentMap = significanceLevel2ElementNameTag2NumberofEnrichmentMapMap.get(significanceLevel);
					
					//You have to initialize
					if (!elementNameTag2NumberofEnrichmentMap.containsKey(elementNameTag)){
						elementNameTag2NumberofEnrichmentMap.put(elementNameTag, 0);
					}
					
					if (pValue <= significanceLevel){							
						if (elementNameTag2NumberofEnrichmentMap.get(elementNameTag)==0){
							elementNameTag2NumberofEnrichmentMap.put(elementNameTag, 1);
						}else{
							elementNameTag2NumberofEnrichmentMap.put(elementNameTag, elementNameTag2NumberofEnrichmentMap.get(elementNameTag)+1);
						}
					}//End of if pValue is less than significanceLevel
					
					//There is no enrichment
					//Nothing to do
				
					
				}//End of for each significanceLevel
				
				//update elementName2NumberofTotalRunsMap
				elementNameTag2NumberofTotalRunsMap.put(elementNameTag,elementNameTag2NumberofTotalRunsMap.get(elementNameTag)+1);

			}//End of for each repressor element
				
							
			//Free space
			activator2EmpPValueMap = null;
			repressor2EmpPValueMap = null;

		}//End of For each run number
		
		//Check
		if (numberofFilesReadforHistone!=numberofRuns){
			System.out.println("numberofFilesReadforHistone: " + numberofFilesReadforHistone);
		}
		
		//Check
		if (numberofFilesReadforTF!=numberofRuns){
			System.out.println("numberofFilesReadforTF: " + numberofFilesReadforTF);
		}
		
		
//		//Debug starts
//		for (int i=0; i<significanceLevelList.size();i++){
//			
//			significanceLevel = significanceLevelList.get(i);
//			
//			System.out.println("significanceLevel" + "\t" + significanceLevel + "\t" + "numberof ElementNameTag"  + "\t" + significanceLevel2ElementNameTag2NumberofEnrichmentMapMap.get(significanceLevel).size());
//			
//		}
//		//Debug ends
		

	}


	
	public static void calculateTypeIErrorAndPower(
			ToolType toolType,
			DataDrivenExperimentCellLineType cellLine,
			DataDrivenExperimentGeneType geneType,			
			TFloatList significanceLevelList,
			List<DataDrivenExperimentElementNameType> activatorElementList,
			List<DataDrivenExperimentElementNameType> repressorElementList,
			TObjectIntMap<String> elementNameTag2NumberofTotalRunsMap,
			TFloatObjectMap<TObjectIntMap<String>> significanceLevel2ElementNameTag2NumberofEnrichmentMapMap,
			BufferedWriter bufferedWriterTypeIError,
			BufferedWriter bufferedWriterPower) throws IOException{
		
		
		TObjectIntMap<String> elementNameTag2NumberofEnrichmentMap = null;
		
		int numberofEnrichment = 0;
		int numberofTotalRuns =0;
		
		Float typeIError = 0.0f;
		Float power = 0.0f;
		
		Float significanceLevel = 0.0f;
		DataDrivenExperimentElementNameType elementNameType = null;
		
		String tag = null;
		String elementName = null;
		String elementNameTag = null;
		
		int indexofFirstUnderscore = 0;
		
		for(int i = 0; i<significanceLevelList.size(); i++){
			
			significanceLevel = significanceLevelList.get(i);
			
			elementNameTag2NumberofEnrichmentMap = significanceLevel2ElementNameTag2NumberofEnrichmentMapMap.get(significanceLevel);
			
			//debug starts
			//System.out.println("significanceLevel:"+ "\t" + significanceLevel + "\t" + "elementNameTag2NumberofEnrichmentMap.size:" +  "\t" + elementNameTag2NumberofEnrichmentMap.size());
			//debug ends
			
			for(TObjectIntIterator<String> itr=elementNameTag2NumberofEnrichmentMap.iterator();itr.hasNext();){
				
				itr.advance();
				elementNameTag =itr.key();
				
				numberofEnrichment = elementNameTag2NumberofEnrichmentMap.get(elementNameTag);
				numberofTotalRuns = elementNameTag2NumberofTotalRunsMap.get(elementNameTag);
			
				
				indexofFirstUnderscore = elementNameTag.indexOf('_');				
				elementName = elementNameTag.substring(0, indexofFirstUnderscore);
				elementNameType = DataDrivenExperimentElementNameType.convertStringtoEnum(elementName);
				tag = elementNameTag.substring(indexofFirstUnderscore+1);
				
				
				if (activatorElementList.contains(elementNameType)){
											
					//Scenario: Expressing Genes
					if (geneType.isExpressingProteinCodingGenes()){		
						
						power = (numberofEnrichment*1.0f)/numberofTotalRuns;
						
						bufferedWriterPower.write(
								toolType.convertEnumtoString()+ "\t" +
								cellLine.convertEnumtoString() + "\t" + 
								geneType.convertEnumtoShortString() + "\t" + 
								significanceLevel + "\t" + 
								elementName + "\t" +  
								tag + "\t" + 
								numberofEnrichment + "\t" + 
								numberofTotalRuns + "\t" + 
								power  + System.getProperty("line.separator")); 
					} 
					
					//Scenario: NonExpressing Genes
					else if (geneType.isNonExpressingProteinCodingGenes()){					
						
						typeIError = (numberofEnrichment*1.0f)/numberofTotalRuns;
							
						bufferedWriterTypeIError.write(
								toolType.convertEnumtoString()+ "\t" +
								cellLine.convertEnumtoString() + "\t" + 
								geneType.convertEnumtoShortString() + "\t" + 
								significanceLevel + "\t" + 
								elementName + "\t" + 
								tag + "\t" + 
								numberofEnrichment + "\t" + 
								numberofTotalRuns + "\t" + 
								typeIError  + System.getProperty("line.separator")); 
						
					}
					
				}//End of if element is activator
				
				else if (repressorElementList.contains(elementNameType)) {
					
					//Scenario: Expressing Genes
					if (geneType.isExpressingProteinCodingGenes()){		
						
						typeIError = (numberofEnrichment*1.0f)/numberofTotalRuns;
						
						bufferedWriterTypeIError.write(
								toolType.convertEnumtoString()+ "\t" +
								cellLine.convertEnumtoString() + "\t" + 
								geneType.convertEnumtoShortString() + "\t" + 
								significanceLevel + "\t" + 
								elementName + "\t" + 
								tag + "\t" + 
								numberofEnrichment + "\t" + 
								numberofTotalRuns + "\t" + 
								typeIError  + System.getProperty("line.separator")); 
					} 
					
					//Scenario: NonExpressing Genes
					else if (geneType.isNonExpressingProteinCodingGenes()){					
						
						power = (numberofEnrichment*1.0f)/numberofTotalRuns;
						
						bufferedWriterPower.write(
								toolType.convertEnumtoString()+ "\t" +
								cellLine.convertEnumtoString() + "\t" + 
								geneType.convertEnumtoShortString() + "\t" + 
								significanceLevel + "\t" + 
								elementName + "\t" + 
								tag + "\t" + 
								numberofEnrichment + "\t" + 
								numberofTotalRuns + "\t" + 
								power  + System.getProperty("line.separator")); 
						
					}
					
				}//End of if element is repressor
				
			}//End of for each elementNameTag
		
			
		}//End of for each significanceLevel
				
	}

	//Not called anymore
	public static int calculateNumberofTotalRuns(
			int numberofRuns,
			List<DataDrivenExperimentCellLineType> cellLineList,
			List<DataDrivenExperimentGeneType> geneTypeList,
			List<DataDrivenExperimentDnaseOverlapExclusionType> dnaseExclusionTypeList,
			List<DataDrivenExperimentTPMType> tpmTypeList,
			List<AssociationMeasureType> associationMeasureTypeList,
			List<IsochoreFamilyMode> isochoreFamilyModeList,
			List<GenerateRandomDataMode> generateRandomDataModeList){
		
		//Initialzie
		int numberofTotalRuns = numberofRuns;
		
		numberofTotalRuns = numberofTotalRuns * cellLineList.size();
		//We must consider only one geneType
		numberofTotalRuns = numberofTotalRuns * geneTypeList.size();		
		numberofTotalRuns = numberofTotalRuns * dnaseExclusionTypeList.size();
		numberofTotalRuns = numberofTotalRuns * tpmTypeList.size();
		numberofTotalRuns = numberofTotalRuns * associationMeasureTypeList.size();
		numberofTotalRuns = numberofTotalRuns * isochoreFamilyModeList.size();
		numberofTotalRuns = numberofTotalRuns * generateRandomDataModeList.size();
		
		return numberofTotalRuns;
		
	}
	
//	public static void fillPossibleTagsList(ToolType toolType,List<String> possibleTags){
//		
//		possibleTags.clear();
//		
//		switch(toolType){
//		
//			case GLANET:
//				possibleTags.add("EOO_woIF_wGC");
//				possibleTags.add("EOO_woIF_wM");
//				possibleTags.add("EOO_woIF_wGCM");
//				possibleTags.add("EOO_woIF_woGCM");
//
//				possibleTags.add("EOO_wIF_wGC");
//				possibleTags.add("EOO_wIF_wM");
//				possibleTags.add("EOO_wIF_wGCM");
//				possibleTags.add("EOO_wIF_woGCM");
//				
//				possibleTags.add("NOOB_woIF_wGC");
//				possibleTags.add("NOOB_woIF_wM");
//				possibleTags.add("NOOB_woIF_wGCM");
//				possibleTags.add("NOOB_woIF_woGCM");
//
//				possibleTags.add("NOOB_wIF_wGC");
//				possibleTags.add("NOOB_wIF_wM");
//				possibleTags.add("NOOB_wIF_wGCM");
//				possibleTags.add("NOOB_wIF_woGCM");
//				break;
//				
//			case GAT:
//				possibleTags.add("EOO_wIF_woGCM_GAT");
//				possibleTags.add("EOO_woIF_woGCM_GAT");
//				possibleTags.add("NOOB_wIF_woGCM_GAT");
//				possibleTags.add("NOOB_woIF_woGCM_GAT");
//				break;
//			default:
//				break;
//		
//		}//End of swicth
//		
//		
//	}
	
	
	public static void prepareDataFileForTypeIErrorFigures(
			ToolType toolType,
			String mainDirectory,
			int numberofRuns,
			DataDrivenExperimentCellLineType cellLine,
			DataDrivenExperimentGeneType geneType,
			DataDrivenExperimentDnaseOverlapExclusionType exclusionType,
			DataDrivenExperimentTPMType tpmType,
			AssociationMeasureType measureType,
			IsochoreFamilyMode isochoreFamilyMode,
			GenerateRandomDataMode generateRandomDataMode){
		
		
		FileWriter fileWriterTypeIError = null;
		BufferedWriter bufferedWriterTypeIError = null;
		
		FileWriter fileWriterPower = null;
		BufferedWriter bufferedWriterPower = null;

		String distinguisingFileName = "";
		String runDirectory = null;
		
		//Create instances
		List<DataDrivenExperimentElementNameType> activatorElementList = new ArrayList<DataDrivenExperimentElementNameType>();
		List<DataDrivenExperimentElementNameType> repressorElementList = new ArrayList<DataDrivenExperimentElementNameType>();
		TFloatList significanceLevelList = new TFloatArrayList();		
		
		TObjectIntMap<String> elementNameTag2NumberofTotalRunsMap = new TObjectIntHashMap<String>();
		TFloatObjectMap<TObjectIntMap<String>> significanceLevel2ElementNameTag2NumberofEnrichmentMapMap = new TFloatObjectHashMap<TObjectIntMap<String>>();

				

		/**************************************************************************************/
		//Set 1 considered cellLines
		List<DataDrivenExperimentCellLineType> cellLineList= new ArrayList<DataDrivenExperimentCellLineType>();
		if (cellLine!=null){
			cellLineList.add(cellLine);		
			distinguisingFileName = distinguisingFileName + cellLine.convertEnumtoString();
		}else{
			//pool for them
			cellLineList.add(DataDrivenExperimentCellLineType.GM12878);
			cellLineList.add(DataDrivenExperimentCellLineType.K562);						
			distinguisingFileName = distinguisingFileName + "Cell";			
		}		
		/**************************************************************************************/
		
		/**************************************************************************************/
		//Set 2 considered geneTypes
		List<DataDrivenExperimentGeneType> geneTypeList= new ArrayList<DataDrivenExperimentGeneType>();
		if (geneType!=null){
			geneTypeList.add(geneType);	
			distinguisingFileName = distinguisingFileName +  "_" + geneType.convertEnumtoShortString() ;
		}else{
			//pool for them
			geneTypeList.add(DataDrivenExperimentGeneType.NONEXPRESSING_PROTEINCODING_GENES);
			geneTypeList.add(DataDrivenExperimentGeneType.EXPRESSING_PROTEINCODING_GENES);			
			distinguisingFileName = distinguisingFileName + "_" + "Gene";			
		}
		/**************************************************************************************/
		

		/**************************************************************************************/
		//Set 3 considered dnase exclusion Types
		List<DataDrivenExperimentDnaseOverlapExclusionType> NE_DnaseExclusionTypeList= new ArrayList<DataDrivenExperimentDnaseOverlapExclusionType>();
		List<DataDrivenExperimentDnaseOverlapExclusionType> E_DnaseExclusionTypeList= new ArrayList<DataDrivenExperimentDnaseOverlapExclusionType>();
		List<DataDrivenExperimentDnaseOverlapExclusionType> dnaseExclusionTypeList= new ArrayList<DataDrivenExperimentDnaseOverlapExclusionType>();
		
		if (exclusionType!=null){
			
			if (geneTypeList.size()==1 && geneTypeList.contains(DataDrivenExperimentGeneType.NONEXPRESSING_PROTEINCODING_GENES)){
				NE_DnaseExclusionTypeList.add(exclusionType);				
			}else if (geneTypeList.size()==1 && geneTypeList.contains(DataDrivenExperimentGeneType.EXPRESSING_PROTEINCODING_GENES)){
				E_DnaseExclusionTypeList.add(exclusionType);				
			}
			
			distinguisingFileName = distinguisingFileName + "_" + exclusionType.convertEnumtoString() ;

						
		}else{
			//pool for them
			NE_DnaseExclusionTypeList.add(DataDrivenExperimentDnaseOverlapExclusionType.COMPLETELY_DISCARD_INTERVAL);
			NE_DnaseExclusionTypeList.add(DataDrivenExperimentDnaseOverlapExclusionType.PARTIALLY_DISCARD_INTERVAL_TAKE_THE_LONGEST_REMAINING_INTERVAL);
			
			E_DnaseExclusionTypeList.add(DataDrivenExperimentDnaseOverlapExclusionType.NO_DISCARD);
			
			
			distinguisingFileName = distinguisingFileName +  "_" + "Discard";
			
		}
		/**************************************************************************************/

	
		/**************************************************************************************/
		//Set 4 considered tpmType
		List<DataDrivenExperimentTPMType> NE_tpmTypeList= new ArrayList<DataDrivenExperimentTPMType>();
		List<DataDrivenExperimentTPMType> E_tpmTypeList= new ArrayList<DataDrivenExperimentTPMType>();		
		List<DataDrivenExperimentTPMType> tpmTypeList= new ArrayList<DataDrivenExperimentTPMType>();
		
		if (tpmType!=null){
			
			if (geneTypeList.size()==1 && geneTypeList.contains(DataDrivenExperimentGeneType.NONEXPRESSING_PROTEINCODING_GENES)){
				NE_tpmTypeList.add(tpmType);
			}else if (geneTypeList.size()==1 && geneTypeList.contains(DataDrivenExperimentGeneType.EXPRESSING_PROTEINCODING_GENES)){
				E_tpmTypeList.add(tpmType);
			}			
			distinguisingFileName = distinguisingFileName  + "_" + tpmType.convertEnumtoString();
		
		}else{
			//pool for them
			NE_tpmTypeList.add(DataDrivenExperimentTPMType.TOPUNKNOWN);
			
			E_tpmTypeList.add(DataDrivenExperimentTPMType.TOP5);
			E_tpmTypeList.add(DataDrivenExperimentTPMType.TOP20);
			
			distinguisingFileName = distinguisingFileName  + "_" + "Top";			
			
		}
		/**************************************************************************************/

		/**************************************************************************************/
		//Set 5 considered associationMeasureType
		List<AssociationMeasureType> associationMeasureTypeList= new ArrayList<AssociationMeasureType>();
		
		if (measureType!=null){			
			associationMeasureTypeList.add(measureType);	
			distinguisingFileName = distinguisingFileName + "_" + measureType.convertEnumtoShortString() ;

		}else{
			//pool for them
			associationMeasureTypeList.add(AssociationMeasureType.EXISTENCE_OF_OVERLAP);
			associationMeasureTypeList.add(AssociationMeasureType.NUMBER_OF_OVERLAPPING_BASES);
			
			distinguisingFileName = distinguisingFileName + "_" + "Measure";
			
		}
		/**************************************************************************************/

		/**************************************************************************************/
		//Set 6 considered isochoreFamilyMode
		List<IsochoreFamilyMode> isochoreFamilyModeList= new ArrayList<IsochoreFamilyMode>();
		
		if (isochoreFamilyMode!=null){
			isochoreFamilyModeList.add(isochoreFamilyMode);
			distinguisingFileName = distinguisingFileName + "_" + isochoreFamilyMode.convertEnumtoShortString();

		}else{
			//pool for them
			isochoreFamilyModeList.add(IsochoreFamilyMode.DO_USE_ISOCHORE_FAMILY);
			isochoreFamilyModeList.add(IsochoreFamilyMode.DO_NOT_USE_ISOCHORE_FAMILY);
			
			distinguisingFileName = distinguisingFileName + "_" + "IsochoreFamilyMode";
	
		}
		/**************************************************************************************/
		
		/**************************************************************************************/
		//Set 7 considered generateRandomDataMode
		List<GenerateRandomDataMode> generateRandomDataModeList= new ArrayList<GenerateRandomDataMode>();
		
		if (generateRandomDataMode!=null){
			generateRandomDataModeList.add(generateRandomDataMode);
			distinguisingFileName = distinguisingFileName + "_" + generateRandomDataMode.convertEnumtoShortString();

		}else{
			//pool for them
			generateRandomDataModeList.add(GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITH_GC_CONTENT);
			generateRandomDataModeList.add(GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITH_MAPPABILITY);
			generateRandomDataModeList.add(GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITH_MAPPABILITY_AND_GC_CONTENT);
			generateRandomDataModeList.add(GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITHOUT_MAPPABILITY_AND_GC_CONTENT);
			
			distinguisingFileName = distinguisingFileName + "_" + "GenerateRandomDataMode";
			
		}
		/**************************************************************************************/
		
		

		try {
			fileWriterTypeIError = FileOperations.createFileWriter(mainDirectory + System.getProperty("file.separator") + toolType.convertEnumtoString() + "_TypeIError_Power_Data" + System.getProperty("file.separator") + toolType.convertEnumtoString() +  "_TypeI_Error_Data_" +  distinguisingFileName + ".txt");
			bufferedWriterTypeIError = new BufferedWriter(fileWriterTypeIError);

			fileWriterPower = FileOperations.createFileWriter(mainDirectory + System.getProperty("file.separator") + toolType.convertEnumtoString() + "_TypeIError_Power_Data" + System.getProperty("file.separator") + toolType.convertEnumtoString() +"_Power_Data_" +  distinguisingFileName + ".txt");
			bufferedWriterPower = new BufferedWriter(fileWriterPower);

			//Header line
			bufferedWriterTypeIError.write(
					"Tool" + "\t" +
					"CellLine" + "\t" + 
					"Scenario" + "\t" + 
					"SignificanceLevel" + "\t" + 
					"ElementName" + "\t" + 
					"Method" + "\t" + 
					"NumberofEnrichment" + "\t" + 
					"NumberofTotalRuns" + "\t" + 
					"TypeIError"  + System.getProperty("line.separator")); 	
	
			//Header line
			bufferedWriterPower.write(
					"Tool" + "\t" +
					"CellLine" + "\t" + 
					"Scenario" + "\t" + 
					"SignificanceLevel" + "\t" + 
					"ElementName" + "\t" + 
					"Method" + "\t" + 
					"NumberofEnrichment" + "\t" + 
					"NumberofTotalRuns" + "\t" + 
					"Power"  + System.getProperty("line.separator")); 	

			
			for(Iterator<DataDrivenExperimentCellLineType> cellLineItr = cellLineList.iterator();cellLineItr.hasNext();) {
								
				cellLine = cellLineItr.next();
				
				//Initialize
				activatorElementList.clear();
				repressorElementList.clear();
				elementNameTag2NumberofTotalRunsMap.clear();
				fillWithDataDrivenExperimentElementNameTypesIncludingAmbigiousElements(
						cellLine,
						activatorElementList,
						repressorElementList);

				//Initialize
				significanceLevelList.clear();
				significanceLevel2ElementNameTag2NumberofEnrichmentMapMap.clear();				
				fillMaps(
						significanceLevelList,
						significanceLevel2ElementNameTag2NumberofEnrichmentMapMap);

			
				for(Iterator<DataDrivenExperimentGeneType> geneTypeItr = geneTypeList.iterator(); geneTypeItr.hasNext();){
					
					geneType = geneTypeItr.next();
					
					if (geneType.isNonExpressingProteinCodingGenes()){
						
						dnaseExclusionTypeList = NE_DnaseExclusionTypeList;	
						tpmTypeList = NE_tpmTypeList;
						
					}else if (geneType.isExpressingProteinCodingGenes()){	
						
						dnaseExclusionTypeList = E_DnaseExclusionTypeList;
						tpmTypeList = E_tpmTypeList;
					}
					
					
					for(Iterator<DataDrivenExperimentTPMType> tpmTypeItr = tpmTypeList.iterator();tpmTypeItr.hasNext();){
						
						tpmType = tpmTypeItr.next();
						
						for(Iterator<DataDrivenExperimentDnaseOverlapExclusionType> dnaseExclusionTypeItr = dnaseExclusionTypeList.iterator(); dnaseExclusionTypeItr.hasNext();){
													
							exclusionType = dnaseExclusionTypeItr.next();											
							
							for(Iterator<AssociationMeasureType> measureTypeItr = associationMeasureTypeList.iterator();measureTypeItr.hasNext();){
								
								measureType = measureTypeItr.next();
																
								switch(toolType){
								
									case GLANET:
										/*******************************************************************************/
										/*******************************GLANET starts***********************************/
										/*******************************************************************************/
										for(Iterator<IsochoreFamilyMode> isochoreFamilyModeItr =  isochoreFamilyModeList.iterator();isochoreFamilyModeItr.hasNext();) {
											

											isochoreFamilyMode = isochoreFamilyModeItr.next();
											
											for(Iterator<GenerateRandomDataMode> generateRandomDataModeItr = generateRandomDataModeList.iterator(); generateRandomDataModeItr.hasNext();){
																								
												generateRandomDataMode =generateRandomDataModeItr.next();
																																	
												runDirectory= "Output" + System.getProperty("file.separator") + cellLine + "_" + geneType.convertEnumtoString() + "_" + tpmType.convertEnumtoString() + "_" + exclusionType.convertEnumtoString() + "_" + generateRandomDataMode.convertEnumtoShortString() + "_" + isochoreFamilyMode.convertEnumtoShortString() + "_" + measureType.convertEnumtoShortString() + "Run";
												
												prepareTabDelimitedDataFileForGLANETTypeIErrorAndPowerFigures(
														mainDirectory, 
														runDirectory, 
														numberofRuns, 
														exclusionType,
														tpmType,
														measureType,
														isochoreFamilyMode,
														generateRandomDataMode,
														significanceLevelList,
														activatorElementList, 
														repressorElementList,
														elementNameTag2NumberofTotalRunsMap,
														significanceLevel2ElementNameTag2NumberofEnrichmentMapMap);																									
											
											}//End of for each generateRandomDataMode														
											
										}//End of each isochoreFamilyMode
										/*******************************************************************************/
										/*******************************GLANET ends*************************************/
										/*******************************************************************************/									
										break;
										
									case GAT:
										/*******************************************************************************/
										/*********************************GAT starts************************************/
										/*******************************************************************************/									
										//Initialization
										//For GAT I wrote wIF and woIF by myself in createTag method
										//Because GAT output filenames does not has wIF woIF in their names.
										generateRandomDataModeList.clear();
										
										//In fact GAT only achieves wIF and woIF
										//For wIF I have named the files with wGC									
										//For woIF I have named the files with woGCM
										generateRandomDataModeList.add(GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITH_GC_CONTENT);
										generateRandomDataModeList.add(GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITHOUT_MAPPABILITY_AND_GC_CONTENT);
										
										for(Iterator<GenerateRandomDataMode> generateRandomDataModeItr = generateRandomDataModeList.iterator(); generateRandomDataModeItr.hasNext();){
											
											generateRandomDataMode =generateRandomDataModeItr.next();
																																	
											//GAT DDE Output EOO uses wGC and woGCM in filenames for wGC and woGCM, respectively.
											//GAT DDE Output NOOB uses wGCM and woGCM in filenames for wGC and woGCM, respectively.			
									
											//My naming --> Consistent Naming with GLANET
											//EOO_wGC --> EOO_wIF_woGCM
											//EOO_woGCM --> EOO_woIF_woGCM											
											//NOOB_wGCM --> NOOB_wIF_woGCM
											//NOOB_woGCM --> NOOB_woIF_woGCM
											
											if (measureType.isAssociationMeasureNumberOfOverlappingBases() && generateRandomDataMode.isGenerateRandomDataModeWithGC()){
												runDirectory= "Output" + System.getProperty("file.separator") + toolType.convertEnumtoString() + "_" +  cellLine + "_" + geneType.convertEnumtoString() + "_" + tpmType.convertEnumtoString() + "_" + exclusionType.convertEnumtoString() + "_" + GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITH_MAPPABILITY_AND_GC_CONTENT.convertEnumtoShortString() + "_" + measureType.convertEnumtoShortString() + "_" + Commons.DDE_RUN;												
											}else{
												runDirectory= "Output" + System.getProperty("file.separator") + toolType.convertEnumtoString() + "_" +  cellLine + "_" + geneType.convertEnumtoString() + "_" + tpmType.convertEnumtoString() + "_" + exclusionType.convertEnumtoString() + "_" + generateRandomDataMode.convertEnumtoShortString() + "_" + measureType.convertEnumtoShortString() + "_" + Commons.DDE_RUN;							
											}
											
											prepareTabDelimitedDataFileForGATTypeIErrorAndPowerFigures(
													mainDirectory, 
													runDirectory, 
													numberofRuns, 
													exclusionType,
													tpmType,
													measureType,
													generateRandomDataMode,
													significanceLevelList,
													activatorElementList, 
													repressorElementList,
													elementNameTag2NumberofTotalRunsMap,
													significanceLevel2ElementNameTag2NumberofEnrichmentMapMap);								
										
										}//End of for each generateRandomDataMode											
										/*******************************************************************************/
										/*********************************GAT ends**************************************/
										/*******************************************************************************/
										break;
										
									default:
										break;
								
								}//End of switch

							}//End of each associationMeasureType							
							
						}//End of for each dnaseExclusionType
						
					}//End of for each tpmType	
					
					
					
					
									
					//Calculate TypeIError Rates
					calculateTypeIErrorAndPower(
							toolType, 
							cellLine,
							geneType, 
							significanceLevelList,
							activatorElementList,
							repressorElementList,
							elementNameTag2NumberofTotalRunsMap,
							significanceLevel2ElementNameTag2NumberofEnrichmentMapMap,
							bufferedWriterTypeIError,
							bufferedWriterPower);

				}//End of for each geneType
				
			}//End of for each cellLine
		
			
		
			//Close
			bufferedWriterTypeIError.close();
			bufferedWriterPower.close();
			
		}catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}

	//My naming --> Consistent Naming with GLANET
	//EOO_wGC --> EOO_wIF_woGCM
	//EOO_woGCM --> EOO_woIF_woGCM											
	//NOOB_wGCM --> NOOB_wIF_woGCM
	//NOOB_woGCM --> NOOB_woIF_woGCM	
	public static String getGLANETConsistentTag(
			AssociationMeasureType measureType,
			GenerateRandomDataMode generateRandomDataMode){
		
		String tag = null;
		
		switch(generateRandomDataMode){
		
			//wGC
			//wGCM			
			case GENERATE_RANDOM_DATA_WITH_GC_CONTENT:
			case GENERATE_RANDOM_DATA_WITH_MAPPABILITY_AND_GC_CONTENT:					
				tag = 	measureType.convertEnumtoShortString() + "_" +
						IsochoreFamilyMode.DO_USE_ISOCHORE_FAMILY.convertEnumtoShortString() + "_" +
						GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITHOUT_MAPPABILITY_AND_GC_CONTENT.convertEnumtoShortString();
				break;
				
			
			//woGCM
			case GENERATE_RANDOM_DATA_WITHOUT_MAPPABILITY_AND_GC_CONTENT:
				tag = 	measureType.convertEnumtoShortString() + "_" +
						IsochoreFamilyMode.DO_NOT_USE_ISOCHORE_FAMILY.convertEnumtoShortString() + "_" +
						generateRandomDataMode.convertEnumtoShortString();
				break;
			
			default:
				break;
		
		}//End of switch
		
	
		return tag;
		
	}
	

	
	//For collecting data for GAT TypeIErrors and Power Figures
	public static void prepareTabDelimitedDataFileForGATTypeIErrorAndPowerFigures(
			String mainDirectory,
			String runDirectory,
			int numberofRuns, 
			DataDrivenExperimentDnaseOverlapExclusionType exclusionType,
			DataDrivenExperimentTPMType tpmType,
			AssociationMeasureType measureType,
			GenerateRandomDataMode generateRandomDataMode, 
			TFloatList significanceLevelList,
			List<DataDrivenExperimentElementNameType> activatorElementList, 
			List<DataDrivenExperimentElementNameType> repressorElementList,
			TObjectIntMap<String> elementNameTag2NumberofTotalRunsMap,
			TFloatObjectMap<TObjectIntMap<String>> significanceLevel2ElementNameTag2NumberofEnrichmentMapMap){
						
		String strLine = null;
		File gatTSVFile = null;
		
		FileReader gatTSVFileReader = null;
		BufferedReader gatTSVBufferedReader = null;
		
		//My naming --> Consistent Naming with GLANET
		//EOO_wGC --> EOO_wIF_woGCM
		//EOO_woGCM --> EOO_woIF_woGCM											
		//NOOB_wGCM --> NOOB_wIF_woGCM
		//NOOB_woGCM --> NOOB_woIF_woGCM		
		String tag = getGLANETConsistentTag(measureType,generateRandomDataMode);
		
		//Add exclusionType
		//Add tpmType
		tag = exclusionType.convertEnumtoString() + "_" + tpmType.convertEnumtoString() + "_" +  tag; 
		
		//Add GAT
		tag = tag + "_" + ToolType.GAT.convertEnumtoString();
		
		String elementNameTag = null;
		
		DataDrivenExperimentElementNameType elementName = null;
		Float pValue = 0f;
		
		Float significanceLevel = null;
		TObjectIntMap<String> elementNameTag2NumberofEnrichmentMap;
	
			
		try{

			// For each run
			for(int i = 0; i <numberofRuns; i++){
				
				//Create Instances
				TObjectFloatMap<DataDrivenExperimentElementNameType> activator2EmpPValueMap  =  new TObjectFloatHashMap<DataDrivenExperimentElementNameType>();
				TObjectFloatMap<DataDrivenExperimentElementNameType> repressor2EmpPValueMap  =  new TObjectFloatHashMap<DataDrivenExperimentElementNameType>();
				
								
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
						
						processLine(
								strLine,
								activatorElementList, 
								repressorElementList, 
								activator2EmpPValueMap,
								repressor2EmpPValueMap);
						
					}// End of WHILE
					
					
					// Close
					gatTSVBufferedReader.close();
					
				}//End of IF gatTSVFile exists
				
				else{
					//Write down this file does not exists
					//Exit from loop
					System.out.println("There is no file for this run: " + i + "\t" +  gatTSVFile.getAbsolutePath());				
				}
				
				//In fact not necessary for GAT
				//For all missing activator elements put a pValue of 1.
				for(Iterator<DataDrivenExperimentElementNameType> itr =activatorElementList.iterator();itr.hasNext();){
					
					elementName = itr.next();
					
					if (!activator2EmpPValueMap.containsKey(elementName)){
						activator2EmpPValueMap.put(elementName,1.0f);					
					}
					
				}//End of FOR
				
				
				//In fact not necessary for GAT
				//For all missing repressor elements put a pValue of 1.
				for(Iterator<DataDrivenExperimentElementNameType> itr =repressorElementList.iterator();itr.hasNext();){
					
					elementName = itr.next();
					
					if (!repressor2EmpPValueMap.containsKey(elementName)){
						repressor2EmpPValueMap.put(elementName,1.0f);
					}
					
				}//End of FOR
				
				
				//Fill significanceLevel2ElementName2NumberofEnrichmentMapMap
				for(int j=0; j<activatorElementList.size();j++){
					
					elementName = activatorElementList.get(j);	
					elementNameTag = elementName.convertEnumtoString() + "_" + tag;
					pValue = activator2EmpPValueMap.get(elementName);
					
					//For each value of significanceLevel
					for(int k=0; k<significanceLevelList.size();k++){
						
						significanceLevel = significanceLevelList.get(k);
						
						elementNameTag2NumberofEnrichmentMap = significanceLevel2ElementNameTag2NumberofEnrichmentMapMap.get(significanceLevel);
						
						//You have to initialize
						if (!elementNameTag2NumberofEnrichmentMap.containsKey(elementNameTag)){
							elementNameTag2NumberofEnrichmentMap.put(elementNameTag, 0);
						}
						
						
						if (pValue <= significanceLevel){							
							if (elementNameTag2NumberofEnrichmentMap.get(elementNameTag)==0){
								elementNameTag2NumberofEnrichmentMap.put(elementNameTag, 1);
							}else{
								elementNameTag2NumberofEnrichmentMap.put(elementNameTag, elementNameTag2NumberofEnrichmentMap.get(elementNameTag)+1);
							}
						}//End of if pValue is less than significanceLevel
						
					}//End of for each significanceLevel
					
					//update elementName2NumberofTotalRunsMap
					elementNameTag2NumberofTotalRunsMap.put(elementNameTag,elementNameTag2NumberofTotalRunsMap.get(elementNameTag)+1);
					
										
				}//End of for each activator element
				
				for(int j=0; j<repressorElementList.size();j++){
					
					elementName = repressorElementList.get(j);					
					elementNameTag = elementName.convertEnumtoString() + "_" + tag;
					pValue = repressor2EmpPValueMap.get(elementName);
					
					//For each value of significanceLevel
					for(int k=0; k<significanceLevelList.size();k++){
						
						significanceLevel = significanceLevelList.get(k);
						
						elementNameTag2NumberofEnrichmentMap = significanceLevel2ElementNameTag2NumberofEnrichmentMapMap.get(significanceLevel);
						
						//You have to initialize
						if (!elementNameTag2NumberofEnrichmentMap.containsKey(elementNameTag)){
							elementNameTag2NumberofEnrichmentMap.put(elementNameTag, 0);
						}
						
						if (pValue <= significanceLevel){							
							if (elementNameTag2NumberofEnrichmentMap.get(elementNameTag)==0){
								elementNameTag2NumberofEnrichmentMap.put(elementNameTag, 1);
							}else{
								elementNameTag2NumberofEnrichmentMap.put(elementNameTag, elementNameTag2NumberofEnrichmentMap.get(elementNameTag)+1);
							}
						}//End of if pValue is less than significanceLevel
						
					}//End of for each significanceLevel
					
					//update elementName2NumberofTotalRunsMap
					elementNameTag2NumberofTotalRunsMap.put(elementNameTag,elementNameTag2NumberofTotalRunsMap.get(elementNameTag)+1);

				}//End of for each repressor element

				
				//Free space
				activator2EmpPValueMap  =  null;
				repressor2EmpPValueMap  =  null;
				
			}// End of FOR each run
			
			
		}catch( IOException e){
			e.printStackTrace();
		}
		
	}
	
	
	//For collecting GAT TypeI Error and Power Data
	//Fill activator2EmpPValueMap
	//Fill repressor2EmpPValueMap
	public static void processLine(
			String strLine,
			List<DataDrivenExperimentElementNameType> activatorElementList, 
			List<DataDrivenExperimentElementNameType> repressorElementList,
			TObjectFloatMap<DataDrivenExperimentElementNameType> activator2EmpPValueMap,
			TObjectFloatMap<DataDrivenExperimentElementNameType> repressor2EmpPValueMap) throws IOException{
						
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
			
			//Question: Shall I consider ln2Fold?
			//Decision: Yes
			//If ln2fold is positive take empiricalPValue as it is
			//If ln2fold is negative take 1-empiricalPValue
			if(ln2Fold<0){
				empiricalPValue = 1-empiricalPValue;
			}						
			
			if (activatorElementList.contains(elementNameType)){				
				activator2EmpPValueMap.put(elementNameType,empiricalPValue);															
			}else if (repressorElementList.contains(elementNameType)){				
				repressor2EmpPValueMap.put(elementNameType,empiricalPValue);				
			}		
			
		
		}//End of IF valid strLine control
		
	}


	
	public static void main(String[] args) {
		
		//ToolType GLANET or GAT
		ToolType toolType = ToolType.convertStringtoEnum(args[0]);
						
		//MainDirectory excluding Output directory
		//Main directory must be the parent of Output Directory
		String mainDirectory=args[1];
		
		//numberofRuns in each case
		int numberofRuns = Integer.parseInt(args[2]);
		
		//cellLine
		DataDrivenExperimentCellLineType cellLine = DataDrivenExperimentCellLineType.convertStringtoEnum(args[3]);
		
		//geneType 
		//NonExpressed Genes
		//Expressed Genes
		DataDrivenExperimentGeneType geneType = DataDrivenExperimentGeneType.convertStringtoEnum(args[4]);
				
		//dnaseExclusionType 
		//CompletelyDiscard
		//TakeTheLongest
		//NoDiscard
		DataDrivenExperimentDnaseOverlapExclusionType exclusionType = DataDrivenExperimentDnaseOverlapExclusionType.convertStringtoEnum(args[5]);
		
		//TopUnknown
		//Top5
		//Top20
		DataDrivenExperimentTPMType tpmType = DataDrivenExperimentTPMType.convertStringtoEnum(args[6]);
		
		//EOO
		//NOOB
		AssociationMeasureType measureType = AssociationMeasureType.convertStringtoEnum(args[7]);
		
		//isochoreMode
		//wIF
		//woIF
		IsochoreFamilyMode isochoreFamilyMode = IsochoreFamilyMode.convertStringtoEnum(args[8]);
		
		//wGC
		//wM
		//wGCM
		//woGCM
		GenerateRandomDataMode generateRandomDataMode = GenerateRandomDataMode.convertStringtoEnum(args[9]);
		
		prepareDataFileForTypeIErrorFigures(
				toolType,
				mainDirectory,
				numberofRuns,
				cellLine,
				geneType,
				exclusionType,
				tpmType,
				measureType,
				isochoreFamilyMode,
				generateRandomDataMode);

	}

}
