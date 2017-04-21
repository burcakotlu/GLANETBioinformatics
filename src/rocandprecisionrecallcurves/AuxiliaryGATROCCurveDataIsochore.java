/**
 * 
 */
package rocandprecisionrecallcurves;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import auxiliary.FileOperations;
import enumtypes.IsochoreFamilyMode;

/**
 * @author Burçak Otlu
 * @date Mar 5, 2017
 * @project GLANETBioinformatics 
 * 
 * After DataPreparationForROCCurves class has generated the data for ROC curves and precision recall curves
 * 
 * This class uses the files prepared by DataPreparationForROCCurves class
 * It removes CompletelyDiscard and Top5 from the tag.
 * It simplifies and shortens the tag.
 * 
 * Nothing else.
 *
 */
public class AuxiliaryGATROCCurveDataIsochore {
	
	public static String remove_Top5_CompletelyDiscard(String strLine){
		
		/*******************************************/
		int indeofTop5 = strLine.indexOf("Top5");
		
		if (indeofTop5!=-1){	
			strLine = strLine.substring(0, indeofTop5) + strLine.substring(indeofTop5+5);	
		}
		/*******************************************/
		
		/*******************************************/
		int indeofCompletelyDiscard = strLine.indexOf("CompletelyDiscard");
		
		if (indeofCompletelyDiscard!=-1){	
			strLine = strLine.substring(0, indeofCompletelyDiscard) + strLine.substring(indeofCompletelyDiscard+18);	
		}
		/*******************************************/
		
		return strLine;
		
	}
	
	
	public static void readandPrepare(
			String filename, 
			String fileName_wIF, 
			String fileName_woIF){
		
		FileReader filereader = null;
		BufferedReader bufferedReader = null;
		
		FileWriter fileWriter_wIF = null;
		BufferedWriter bufferedWriter_wIF = null;
		
		FileWriter fileWriter_woIF = null;
		BufferedWriter bufferedWriter_woIF = null;
		
		String strLine = null;
		
		try {
			filereader = FileOperations.createFileReader(filename);
			bufferedReader = new BufferedReader(filereader);
			
			fileWriter_wIF = FileOperations.createFileWriter(fileName_wIF);
			bufferedWriter_wIF = new BufferedWriter(fileWriter_wIF);
			
			fileWriter_woIF = FileOperations.createFileWriter(fileName_woIF);
			bufferedWriter_woIF = new BufferedWriter(fileWriter_woIF);
			
			//Skip header line
			strLine = bufferedReader.readLine();
			
			while((strLine = bufferedReader.readLine())!=null){
				
				strLine = remove_Top5_CompletelyDiscard(strLine);
				
				if (strLine.contains("wIF")){
					bufferedWriter_wIF.write(strLine + System.getProperty("line.separator"));					
				}else if (strLine.contains("woIF")){
					bufferedWriter_woIF.write(strLine + System.getProperty("line.separator"));
				}
				
			}//End of while
			
			
			//Close
			bufferedReader.close();
			bufferedWriter_wIF.close();
			bufferedWriter_woIF.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}


	/**
	 * @param args
	 */
	public static void main(String[] args) {

		String filename = args[0];		
		String generateRandomDataMode = "GenerateRandomDataMode" ;
		
		int index = filename.indexOf(generateRandomDataMode);
		
		String fileName_wIF = filename.substring(0, index-1) + "_" + IsochoreFamilyMode.DO_USE_ISOCHORE_FAMILY.convertEnumtoShortString() + "_" + filename.substring(index);
		String fileName_woIF = filename.substring(0, index-1)+ "_" + IsochoreFamilyMode.DO_NOT_USE_ISOCHORE_FAMILY.convertEnumtoShortString() + "_" + filename.substring(index);
	
		readandPrepare(filename,fileName_wIF,fileName_woIF);
	}

}
