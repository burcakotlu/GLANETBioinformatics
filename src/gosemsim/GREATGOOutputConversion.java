/**
 * 
 */
package gosemsim;

import goterms.GOTermsUtility;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import auxiliary.FileOperations;
import common.Commons;

/**
 * @author Burçak Otlu
 * @date Mar 21, 2017
 * @project GLANETBioinformatics 
 * 
 * I have run GREAT with spp.optimal.wgEncodeSydhTfbsK562bGata2UcdAlnRep0_VS_wgEncodeSydhTfbsK562bInputUcdAlnRep1.narrowPeak
 * under C:\Users\Burçak\Google Drive\Data\demo_input_data\RDA2
 * 
 *  GREAT provides the enrcihed GO Term name not their ID
 *  This class provides GO Term ID for each GO Term name.
 *
 */
public class GREATGOOutputConversion {

	public static void readFileAndFill(String greatEnrichedBPGOTermsFile, List<String> great_enriched_BP_GOTermNames){
		
		FileReader fileReader;
		BufferedReader bufferedReader;
		
		String strLine;
		int indexofFirstTab;
		//int indexofSecondTab;
		String goTermName;
		
		try{
			fileReader = FileOperations.createFileReader(greatEnrichedBPGOTermsFile);
			bufferedReader = new BufferedReader(fileReader);
			
			while( ( strLine = bufferedReader.readLine()) != null){
				
				//skip comment lines
				if (!strLine.startsWith("#")){
					
					indexofFirstTab = strLine.indexOf("\t");
					//indexofSecondTab = strLine.indexOf("\t", indexofFirstTab+1);
					
					if (indexofFirstTab!=-1){
						goTermName = strLine.substring(0,indexofFirstTab);
						
						if(!goTermName.isEmpty()){
							great_enriched_BP_GOTermNames.add(goTermName);
						}
						
					}
					
				}//End of if not a comment				
				
			}// End of While

			bufferedReader.close();

		}catch( FileNotFoundException e){
			e.printStackTrace();
		}catch( IOException e){
			e.printStackTrace();
		}
		
		
	}
	
	//GOTermList_GATA2_P= c("GO:0006351","GO:0035854","GO:0045599","GO:0045746","GO:0045766","GO:0045944","GO:0070345","GO:2000178")
	public static String creatGOTermIDString(
			String great_association_rule,
			List<String> great_enriched_BP_GOTermNames_List,
			Map<String,String> GOTermName2GOTermIDMap){
		
		String great_enriched_BP_GOTermIDs = "great_" + great_association_rule + "_Enriched_BP_GOTerms <- c(";
		
		String GOTermName;
		String GOTermID;
		
		boolean forTheFirstTime = true; 
		
		for(int i=0;i<great_enriched_BP_GOTermNames_List.size();i++){
			
			GOTermName = great_enriched_BP_GOTermNames_List.get(i);
			
			GOTermID = GOTermName2GOTermIDMap.get(GOTermName);
			
			if (GOTermID!=null && !GOTermID.isEmpty()){
				
				if (forTheFirstTime){
					forTheFirstTime = false;
					great_enriched_BP_GOTermIDs = great_enriched_BP_GOTermIDs + "\"" + GOTermID + "\"";
				}
				else{
					great_enriched_BP_GOTermIDs = great_enriched_BP_GOTermIDs + ",\"" + GOTermID + "\"";
				}
			}else{
				
				System.out.println("No GO Term ID found for: "  + GOTermName);
			}
			
			
		}//End of for
		
		great_enriched_BP_GOTermIDs = great_enriched_BP_GOTermIDs + ")" ;
		
		return great_enriched_BP_GOTermIDs;
	}

	
	public static void main(String[] args) {
		
		//read the Great output file
		//great_basalplusextension_found_enriched_BP_GO_Terms.txt
		//great_singlenearestgene_found_enriched_BP_GO_Terms.txt
		//great_twonearestgenes_found_enriched_BP_GO_Terms.txt
		
		//String great_association_rule ="twonearestgenes";
		String great_association_rule ="singlenearestgene";
		
		//String filename = "C:\\Users\\Burçak\\Google Drive\\GLANET_Bioinformatics_2ndSubmission\\GREAT\\great_twonearestgenes_found_enriched_BP_GO_Terms.txt";
		String filename = "C:\\Users\\Burçak\\Google Drive\\GLANET_Bioinformatics_2ndSubmission\\GREAT\\great_singlenearestgene_found_enriched_BP_GO_Terms.txt";
		List<String> great_enriched_BP_GOTermNames_List = new ArrayList<String>();
		readFileAndFill(filename,great_enriched_BP_GOTermNames_List);
		
		System.out.println("Number of GO Term Names found enriched by GREAT: " + great_enriched_BP_GOTermNames_List.size());
		
		//Fill GOTermName 2 GOTermID map
		String GOID2TermInputFile = "C:\\Users\\Burçak\\Google Drive\\Data\\" + Commons.GENE_ONTOLOGY_TERMS +  System.getProperty("file.separator") + Commons.GOID2TERM_INPUTFILE;
		Map<String,String> GOTermName2GOTermIDMap = new HashMap<String,String>();
		GOTermsUtility.fillGOTermName2GOTermIDMap(GOID2TermInputFile, GOTermName2GOTermIDMap);
		
		//output the converted GO Term IDs
		//GOTermList_GATA2_P= c("GO:0006351","GO:0035854","GO:0045599","GO:0045746","GO:0045766","GO:0045944","GO:0070345","GO:2000178")
		String great_enriched_BP_GOTermIDs = creatGOTermIDString(great_association_rule,great_enriched_BP_GOTermNames_List,GOTermName2GOTermIDMap);
		
		System.out.println("**************************************************");
		System.out.println(great_enriched_BP_GOTermIDs);

	}

}
