/**
 * 
 */
package wilcoxon;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
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
 * @date Mar 28, 2017
 * @project GLANETBioinformatics 
 * 
 * This class will generate data for wilcoxon rank sum tests.
 * 
 * 
 *
 */
public class DataPreparationforWilcoxonTests {
	

	public static void write(
			DataDrivenExperimentGeneType scenario,
			AssociationMeasureType measure,
			IsochoreFamilyMode isochoreFamily,
			GenerateRandomDataMode generateRandomDataMode,
			TFloatList scenario_measure_isochoreFamily_generateRandomDataMode_typeIErrorList,
			BufferedWriter bufferedWriterWilcoxonTestData) throws IOException{
		
		//E_EOO_woIF_wGC	<- 	c(0,	0,	0,	0.002,	0,	0,	0.103,	0)
		
		DecimalFormat df = new DecimalFormat("#.####");
		int i;
		
		bufferedWriterWilcoxonTestData.write(scenario.convertEnumtoShortString() + "_" + 
												measure.convertEnumtoShortString() + "_" +
												isochoreFamily.convertEnumtoShortString() + "_" +
												generateRandomDataMode.convertEnumtoShortString() + " <- c(");
		
		for(i=0; i<scenario_measure_isochoreFamily_generateRandomDataMode_typeIErrorList.size()-1;i++){			
			bufferedWriterWilcoxonTestData.write(df.format(scenario_measure_isochoreFamily_generateRandomDataMode_typeIErrorList.get(i)) + ",");
		}
		//No comma after the last one
		bufferedWriterWilcoxonTestData.write(df.format(scenario_measure_isochoreFamily_generateRandomDataMode_typeIErrorList.get(i)));

		bufferedWriterWilcoxonTestData.write(")" + System.getProperty("line.separator"));
		
	}

	
	public static void fillRepressorElements(List<DataDrivenExperimentElementNameType> repressorElementList){
		//Repressors 
		repressorElementList.add(DataDrivenExperimentElementNameType.H3K27ME3);				
		repressorElementList.add(DataDrivenExperimentElementNameType.H3K9ME3); 		// <--- Ambigious		
				
	}
	
	
	public static void fillActivatorElements(DataDrivenExperimentCellLineType cell, List<DataDrivenExperimentElementNameType> activatorElementList){
		
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
			
		switch (cell){
		
			case K562:
				//Activators
				activatorElementList.add(DataDrivenExperimentElementNameType.H3K9ACB);		
				break;
				
			default:
				break;
	
		}//End  of switch
		
	}
	
	
	
	
	public static void initializeMaps(
			TFloatList significanceLevelList,
			TFloatObjectMap<TObjectIntMap<DataDrivenExperimentElementNameType>> significanceLevel2ElementName2NumberofEnrichmentMapMap){
		
		Float significanceLevel= null;
		
		for(int i=0; i<=25; i++){
			significanceLevelList.add(0.01f*i);
		}
				
		//Initialize
		for(int i=0; i< significanceLevelList.size(); i++){			
			significanceLevel = significanceLevelList.get(i); 			
			significanceLevel2ElementName2NumberofEnrichmentMapMap.put(significanceLevel,new TObjectIntHashMap<DataDrivenExperimentElementNameType>());			
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
			List<DataDrivenExperimentElementNameType> elementList,
			TObjectFloatMap<DataDrivenExperimentElementNameType> element2EmpPValueMap){
		
		
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
					
					if (elementList.contains(element)){
						
						//gets its empirical p value
						empiricalPValue = getEmpiricalPValue(strLine);
						element2EmpPValueMap.put(element,empiricalPValue);
																	
					}				
					
				}//End of IF not an comment line
				
			}//End of while
			
			//Close
			bufferedReader.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}	

		
	}

	
	public static void prepareDataForTypeIError(
			String mainDirectory, 
			String runDirectory, 
			int numberofRuns, 								
			TFloatList significanceLevelList,
			List<DataDrivenExperimentElementNameType> elementList,
			TObjectIntMap<DataDrivenExperimentElementNameType>  elementName2NumberofTotalRunsMap,
			TFloatObjectMap<TObjectIntMap<DataDrivenExperimentElementNameType>> significanceLevel2ElementName2NumberofEnrichmentMapMap){
		
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
		TObjectIntMap<DataDrivenExperimentElementNameType> elementName2NumberofEnrichmentMap;
		
		
		//For each run
		for(int i=0; i< numberofRuns; i++){
			
			//Create Instances
			TObjectFloatMap<DataDrivenExperimentElementNameType> element2EmpPValueMap  =  new TObjectFloatHashMap<DataDrivenExperimentElementNameType>();
					
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
						//contains("_2017_")
						
						thereIsSuchAFile = true;
						numberofFilesReadforHistone++;

						resultFileName = eachHistoneFile.getAbsolutePath();		
						
						readResultFileAndFillMap(resultFileName,elementList,element2EmpPValueMap);
						
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
						//contains("_2017_")
						
						thereIsSuchAFile = true;
						numberofFilesReadforTF++;

						resultFileName = eachTFFile.getAbsolutePath();		
						
						readResultFileAndFillMap(resultFileName,elementList,element2EmpPValueMap);
					
						break;

					}// End of IF EnrichmentFile under EnrichmentDirectory
					
				}//End of for each file in histone Directory
				
				//Check
				if(!thereIsSuchAFile){
					System.out.println("There is no result file for directory " + tfDirectory);
				}
				
			}//End of IF TF directory exists
			
			//For all missing elements put a pValue of 1.
			for(Iterator<DataDrivenExperimentElementNameType> itr =elementList.iterator();itr.hasNext();){
				
				elementNameType = itr.next();
				
				if (!element2EmpPValueMap.containsKey(elementNameType)){
					element2EmpPValueMap.put(elementNameType,1.0f);					
				}
				
			}//End of FOR
						
			
			//Now accumulate the results of one unique run
			//Fill significanceLevel2ElementName2NumberofEnrichmentMapMap
			for(int j=0; j<elementList.size();j++){
				
				elementNameType = elementList.get(j);	
				pValue = element2EmpPValueMap.get(elementNameType);
				
				//For each value of significanceLevel
				for(int k=0; k<significanceLevelList.size();k++){
					
					significanceLevel = significanceLevelList.get(k);
					
					elementName2NumberofEnrichmentMap = significanceLevel2ElementName2NumberofEnrichmentMapMap.get(significanceLevel);
					
					//You have to initialize
					if (!elementName2NumberofEnrichmentMap.containsKey(elementNameType)){
						elementName2NumberofEnrichmentMap.put(elementNameType, 0);
					}
					
					if (pValue <= significanceLevel){							
						if (elementName2NumberofEnrichmentMap.get(elementNameType)==0){
							elementName2NumberofEnrichmentMap.put(elementNameType, 1);
						}else{
							elementName2NumberofEnrichmentMap.put(elementNameType, elementName2NumberofEnrichmentMap.get(elementNameType)+1);
						}
					}//End of if pValue is less than significanceLevel
					
					//There is no enrichment
					//Nothing to do
					
				}//End of for each significanceLevel
				
				//update elementName2NumberofTotalRunsMap
				elementName2NumberofTotalRunsMap.put(elementNameType,elementName2NumberofTotalRunsMap.get(elementNameType)+1);
				
									
			}//End of for each element
				
							
			//Free space
			element2EmpPValueMap = null;
	
		}//End of For each run number
		
		//Check
		if (numberofFilesReadforHistone!=numberofRuns){
			System.out.println("numberofFilesReadforHistone: " + numberofFilesReadforHistone);
		}
		
		//Check
		if (numberofFilesReadforTF!=numberofRuns){
			System.out.println("numberofFilesReadforTF: " + numberofFilesReadforTF);
		}

		
	}
	
	public static void fillTypeIErrorList(
			TObjectIntMap<DataDrivenExperimentElementNameType> elementName2NumberofTotalRunsMap,
			List<DataDrivenExperimentElementNameType> elementList,
			TFloatList significanceLevelList,
			TFloatObjectMap<TObjectIntMap<DataDrivenExperimentElementNameType>> significanceLevel2ElementName2NumberofEnrichmentMapMap,
			TFloatList scenario_measure_isochoreFamily_generateRandomDataMode_typeIErrorList){
		
		
		Float significanceLevel = null;
		DataDrivenExperimentElementNameType  element = null;
		int numberofEnrichment;
		int numberofTotalRuns;
		Float typeIError = null;
		
		TObjectIntMap<DataDrivenExperimentElementNameType> elementName2NumberofEnrichmentMap = null;
		
		//For each value of significanceLevel
		for(int i=0; i<significanceLevelList.size();i++){
			
			significanceLevel = significanceLevelList.get(i);
			elementName2NumberofEnrichmentMap = significanceLevel2ElementName2NumberofEnrichmentMapMap.get(significanceLevel);
			
			for(int j=0; j<elementList.size(); j++){
				
				element = elementList.get(j);
				
				numberofEnrichment = elementName2NumberofEnrichmentMap.get(element);
				numberofTotalRuns = elementName2NumberofTotalRunsMap.get(element);
				typeIError = (numberofEnrichment*1.0f)/numberofTotalRuns;
				
				scenario_measure_isochoreFamily_generateRandomDataMode_typeIErrorList.add(typeIError);
				
				
			}//End of for each element
			
		}//End of for each significance level

		
		
	}

	
	public static void read(
			String mainDirectory,
			int numberofRuns,
			DataDrivenExperimentGeneType scenario,
			AssociationMeasureType measure,
			IsochoreFamilyMode isochoreFamily,
			GenerateRandomDataMode generateRandomDataMode,
			TFloatList scenario_measure_isochoreFamily_generateRandomDataMode_typeIErrorList){
		
		String runDirectory = null;
		DataDrivenExperimentTPMType TopUnknown = DataDrivenExperimentTPMType.TOPUNKNOWN;
		DataDrivenExperimentDnaseOverlapExclusionType NoDiscard = DataDrivenExperimentDnaseOverlapExclusionType.NO_DISCARD;
		
		//Create instances
		List<DataDrivenExperimentElementNameType> activatorElementList = new ArrayList<DataDrivenExperimentElementNameType>();
		List<DataDrivenExperimentElementNameType> repressorElementList = new ArrayList<DataDrivenExperimentElementNameType>();
		
		TFloatList significanceLevelList = new TFloatArrayList();		
		TObjectIntMap<DataDrivenExperimentElementNameType> elementName2NumberofTotalRunsMap = new TObjectIntHashMap<DataDrivenExperimentElementNameType>();
		TFloatObjectMap<TObjectIntMap<DataDrivenExperimentElementNameType>> significanceLevel2ElementName2NumberofEnrichmentMapMap = new TFloatObjectHashMap<TObjectIntMap<DataDrivenExperimentElementNameType>>();
		
		List<DataDrivenExperimentCellLineType> cellLineList= new ArrayList<DataDrivenExperimentCellLineType>();
		List<DataDrivenExperimentDnaseOverlapExclusionType> NonExp_DnaseExclusionTypeList= new ArrayList<DataDrivenExperimentDnaseOverlapExclusionType>();
		List<DataDrivenExperimentTPMType> Exp_tpmTypeList= new ArrayList<DataDrivenExperimentTPMType>();		
		
		//Fill the lists
		cellLineList.add(DataDrivenExperimentCellLineType.GM12878);
		cellLineList.add(DataDrivenExperimentCellLineType.K562);	
		
		NonExp_DnaseExclusionTypeList.add(DataDrivenExperimentDnaseOverlapExclusionType.COMPLETELY_DISCARD_INTERVAL);
		NonExp_DnaseExclusionTypeList.add(DataDrivenExperimentDnaseOverlapExclusionType.PARTIALLY_DISCARD_INTERVAL_TAKE_THE_LONGEST_REMAINING_INTERVAL);
		
		Exp_tpmTypeList.add(DataDrivenExperimentTPMType.TOP5);
		Exp_tpmTypeList.add(DataDrivenExperimentTPMType.TOP20);
		
		//Get the type I errors
		for(DataDrivenExperimentCellLineType cell :  cellLineList){
			
			switch(scenario){
			
				case EXPRESSING_PROTEINCODING_GENES:
					
					repressorElementList.clear();
					fillRepressorElements(repressorElementList);
					
					for(DataDrivenExperimentTPMType tpmType : Exp_tpmTypeList){
											
						//Initialize
						significanceLevelList.clear();
						significanceLevel2ElementName2NumberofEnrichmentMapMap.clear();	
						elementName2NumberofTotalRunsMap.clear();
						initializeMaps(
								significanceLevelList,
								significanceLevel2ElementName2NumberofEnrichmentMapMap);
						
						runDirectory= "Output" + System.getProperty("file.separator") + cell + "_" + scenario.convertEnumtoString() + "_" + tpmType.convertEnumtoString() + "_" + NoDiscard.convertEnumtoString() + "_" + generateRandomDataMode.convertEnumtoShortString() + "_" + isochoreFamily.convertEnumtoShortString() + "_" + measure.convertEnumtoShortString() + "Run";
					
						//Read Runs
						//For each run fill repressorElement2EmpiricalPValueMap
						//For each run Fill significanceLevel2ElementName2NumberofEnrichmentMapMap
						//For each run Fill elementName2NumberofTotalRunsMap						
						prepareDataForTypeIError(
								mainDirectory, 
								runDirectory, 
								numberofRuns, 								
								significanceLevelList,
								repressorElementList,
								elementName2NumberofTotalRunsMap,
								significanceLevel2ElementName2NumberofEnrichmentMapMap);	
	
						//Finally
						//Fill scenario_measure_isochoreFamily_generateRandomDataMode_typeIErrorList
						fillTypeIErrorList(
								elementName2NumberofTotalRunsMap,
								repressorElementList,
								significanceLevelList,
								significanceLevel2ElementName2NumberofEnrichmentMapMap,
								scenario_measure_isochoreFamily_generateRandomDataMode_typeIErrorList);
						
					}//End of for each tpmType
					break;
					
				case NONEXPRESSING_PROTEINCODING_GENES:
					
					activatorElementList.clear();
					fillActivatorElements(cell,activatorElementList);
					
					
					for(DataDrivenExperimentDnaseOverlapExclusionType exclusionType: NonExp_DnaseExclusionTypeList){
						
						//Initialize
						significanceLevelList.clear();
						significanceLevel2ElementName2NumberofEnrichmentMapMap.clear();				
						elementName2NumberofTotalRunsMap.clear();
						initializeMaps(
								significanceLevelList,
								significanceLevel2ElementName2NumberofEnrichmentMapMap);
						
						runDirectory= "Output" + System.getProperty("file.separator") + cell + "_" + scenario.convertEnumtoString() + "_" + TopUnknown.convertEnumtoString() + "_" + exclusionType.convertEnumtoString() + "_" + generateRandomDataMode.convertEnumtoShortString() + "_" + isochoreFamily.convertEnumtoShortString() + "_" + measure.convertEnumtoShortString() + "Run";
						
						//Read Runs
						//For each run fill activatorElement2EmpiricalPValueMap
						//For each run Fill significanceLevel2ElementNameTag2NumberofEnrichmentMapMap
						//For each run Fill elementNameTag2NumberofTotalRunsMap
						prepareDataForTypeIError(
								mainDirectory, 
								runDirectory, 
								numberofRuns, 								
								significanceLevelList,
								activatorElementList,
								elementName2NumberofTotalRunsMap,
								significanceLevel2ElementName2NumberofEnrichmentMapMap);	
						//Finally
						//Fill scenario_measure_isochoreFamily_generateRandomDataMode_typeIErrorList
						fillTypeIErrorList(
								elementName2NumberofTotalRunsMap,
								activatorElementList,
								significanceLevelList,
								significanceLevel2ElementName2NumberofEnrichmentMapMap,
								scenario_measure_isochoreFamily_generateRandomDataMode_typeIErrorList);
						
					}//End of for each exclusionType
					break;
				
				default:
					break;
					
			}//End of switch for scenario
			
		}//End of for each cell 
		
	}
	

	public static void prepareDataFileForWilcoxonRankSumTests(
			ToolType toolType,
			String mainDirectory,
			int numberofRuns){
		
		FileWriter fileWriterWilcoxonTestData = null;
		BufferedWriter bufferedWriterWilcoxonTestData = null;
					
		List<DataDrivenExperimentGeneType> geneTypeList= new ArrayList<DataDrivenExperimentGeneType>();
		List<AssociationMeasureType> associationMeasureTypeList= new ArrayList<AssociationMeasureType>();
		List<IsochoreFamilyMode> isochoreFamilyModeList= new ArrayList<IsochoreFamilyMode>();
		List<GenerateRandomDataMode> generateRandomDataModeList= new ArrayList<GenerateRandomDataMode>();
		
		//Fill them
		geneTypeList.add(DataDrivenExperimentGeneType.NONEXPRESSING_PROTEINCODING_GENES);
		geneTypeList.add(DataDrivenExperimentGeneType.EXPRESSING_PROTEINCODING_GENES);	
		
		associationMeasureTypeList.add(AssociationMeasureType.EXISTENCE_OF_OVERLAP);
		associationMeasureTypeList.add(AssociationMeasureType.NUMBER_OF_OVERLAPPING_BASES);
		
		isochoreFamilyModeList.add(IsochoreFamilyMode.DO_USE_ISOCHORE_FAMILY);
		isochoreFamilyModeList.add(IsochoreFamilyMode.DO_NOT_USE_ISOCHORE_FAMILY);
		
		generateRandomDataModeList.add(GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITH_GC_CONTENT);
		generateRandomDataModeList.add(GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITH_MAPPABILITY);
		generateRandomDataModeList.add(GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITH_MAPPABILITY_AND_GC_CONTENT);
		generateRandomDataModeList.add(GenerateRandomDataMode.GENERATE_RANDOM_DATA_WITHOUT_MAPPABILITY_AND_GC_CONTENT);
	
		
		TFloatList scenario_measure_isochoreFamily_generateRandomDataMode_typeIErrorList = null;
		
		try {
			fileWriterWilcoxonTestData = FileOperations.createFileWriter(mainDirectory + System.getProperty("file.separator") + toolType.convertEnumtoString() + "_Wilcoxon_Test_Data" + System.getProperty("file.separator") + toolType.convertEnumtoString() +  "_Wilcoxon_Test_Data.txt");
			bufferedWriterWilcoxonTestData = new BufferedWriter(fileWriterWilcoxonTestData);
			
			
			for(DataDrivenExperimentGeneType scenario : geneTypeList){
				
				for(AssociationMeasureType measure: associationMeasureTypeList){
					
					for (IsochoreFamilyMode isochoreFamily : isochoreFamilyModeList){
						
						for(GenerateRandomDataMode generateRandomDataMode: generateRandomDataModeList){
							
							//Initialize
							scenario_measure_isochoreFamily_generateRandomDataMode_typeIErrorList = new TFloatArrayList();
							
							//Read runs data
							//Fill scenario_measure_isochoreFamily_generateRandomDataMode_typeIErrorList
							read(mainDirectory,numberofRuns,scenario,measure,isochoreFamily,generateRandomDataMode,scenario_measure_isochoreFamily_generateRandomDataMode_typeIErrorList);
							
							//Write to file
							write(scenario,
								measure,
								isochoreFamily,
								generateRandomDataMode,
								scenario_measure_isochoreFamily_generateRandomDataMode_typeIErrorList,
								bufferedWriterWilcoxonTestData);
							
						}//End of for each generateRandomDataMode
						
						bufferedWriterWilcoxonTestData.write(System.getProperty("line.separator"));
						
					}//End of for each isochoreFamily
					
				}//End of for eacjh measure
				
			}//End of for each scenario Exp or NonExp
			

			//Close
			bufferedWriterWilcoxonTestData.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}
	
	public static void main(String[] args) {
		
		//ToolType GLANET 
		ToolType toolType = ToolType.convertStringtoEnum(args[0]);
				
		//Main directory must be the parent of Output Directory
		String mainDirectory=args[1];
		
		//numberofRuns in each case
		int numberofRuns = Integer.parseInt(args[2]);
		
		prepareDataFileForWilcoxonRankSumTests(
			toolType,
			mainDirectory,
			numberofRuns);
				
		

	}

}
