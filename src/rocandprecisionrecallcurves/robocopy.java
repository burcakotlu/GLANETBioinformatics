/**
 * 
 */
package rocandprecisionrecallcurves;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import auxiliary.FileOperations;
import enumtypes.AssociationMeasureType;
import enumtypes.DataDrivenExperimentDnaseOverlapExclusionType;
import enumtypes.DataDrivenExperimentGeneType;
import enumtypes.DataDrivenExperimentTPMType;
import enumtypes.IsochoreFamilyMode;

/**
 * @author Burçak Otlu
 * @date Feb 26, 2017
 * @project GLANETBioinformatics 
 * 
 * 
 * This class will write a batch file using robocopy commands.
 *
 */
public class robocopy {
	
	public static void writeCommands(String fileName){
		
		FileWriter fileWriter = null;
		BufferedWriter bufferedWriter = null;
				
		int numberofRuns = 1000;
		
		try {
			
			fileWriter = FileOperations.createFileWriter(fileName);
			bufferedWriter = new BufferedWriter(fileWriter);
			
			List<DataDrivenExperimentTPMType> topList = new ArrayList<DataDrivenExperimentTPMType>();
			List<DataDrivenExperimentDnaseOverlapExclusionType> discardList = new ArrayList<DataDrivenExperimentDnaseOverlapExclusionType>();
			
			for(DataDrivenExperimentGeneType geneType: DataDrivenExperimentGeneType.values()){
				
				switch(geneType){
				
					case NONEXPRESSING_PROTEINCODING_GENES:
						topList.clear();
						discardList.clear();
						
						topList.add(DataDrivenExperimentTPMType.TOPUNKNOWN);
						discardList.add(DataDrivenExperimentDnaseOverlapExclusionType.COMPLETELY_DISCARD_INTERVAL);
						discardList.add(DataDrivenExperimentDnaseOverlapExclusionType.PARTIALLY_DISCARD_INTERVAL_TAKE_THE_LONGEST_REMAINING_INTERVAL);
						break;
						
					case EXPRESSING_PROTEINCODING_GENES:
						topList.clear();
						discardList.clear();
						
						topList.add(DataDrivenExperimentTPMType.TOP5);
						topList.add(DataDrivenExperimentTPMType.TOP20);
						discardList.add(DataDrivenExperimentDnaseOverlapExclusionType.NO_DISCARD);
						break;
					
				}
				
				for (IsochoreFamilyMode isochore:IsochoreFamilyMode.values()){
					
					for (AssociationMeasureType measure : AssociationMeasureType.values()){
						
						for(DataDrivenExperimentTPMType tpm: topList){
							
							for(DataDrivenExperimentDnaseOverlapExclusionType discard: discardList){
								
								for(int i=0; i<numberofRuns; i++){

									bufferedWriter.write("robocopy F:" + System.getProperty("file.separator") + 
												"GLANET_DDCE" + System.getProperty("file.separator") + 
												"GM12878_DDE11_wGC_wM_wGCM_woGCM_wIF_woIF"  + System.getProperty("file.separator") + 
												"Output" + System.getProperty("file.separator") + 
												"GM12878_" + geneType.convertEnumtoString() + "_" + tpm.convertEnumtoString() + "_" + discard.convertEnumtoString() +
												"_woGCM_" + isochore.convertEnumtoShortString() + "_" + measure.convertEnumtoShortString() + "Run" + i + System.getProperty("file.separator")  + " " + 
												"C:" + System.getProperty("file.separator") + 
												"Users" + System.getProperty("file.separator") + 
												"Burçak" +  System.getProperty("file.separator") + 
												"GLANET_DDE11_GM12878_woGCM_Runs" + System.getProperty("file.separator") + 
												"GM12878_" + geneType.convertEnumtoString() + "_" + tpm.convertEnumtoString() + "_" + discard.convertEnumtoString() +
												"_woGCM_" + isochore.convertEnumtoShortString() + "_" + measure.convertEnumtoShortString() + "Run" + i + System.getProperty("file.separator") +
												"."  + " " + "/e" + System.getProperty("line.separator"));
									
								

								}
								
							}							
							
						}						
						
					}
					
				}
				
			}
			
			
			
			//Close
			bufferedWriter.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		String fileName = "C:\\Users\\Burçak\\GLANET_DDE11_GM12878_woGCM_Runs\\robocopy_commands.bat";
		
		writeCommands(fileName);

	}

}
