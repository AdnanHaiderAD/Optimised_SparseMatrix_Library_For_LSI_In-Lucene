package latentSemanticAnalysis;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.Term;
import org.apache.lucene.index.TermFreqVector;
import org.apache.lucene.search.Similarity;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.FSDirectory;

/**@author Adnan 
 * @function constructs the analytics engine by :
 * 1) mapping documents in the analytics lucene index to vector representations where features correspond to individual words.
 * 2) A second level transformation is then applied on the documents using the powerful matrix library to get a more meaningful and 
 * compact representation of docs by projection them to a latent topic space.  
 *  
 */
public class AnalyseDocs {

	private double[] singularValues;
	private SparseMatrix matrixStruct;
	private String dirpath;
	
	/** the following two hashMaps are required to perform quick lookups when comapring documents using LDA engine **/ 
	private HashMap<String,Integer>magnumDocIdLuceneDocIDMap;
	private HashMap <Integer,String>luceneDocIDMagnumDocIDMAP;
	
	public AnalyseDocs(String dirpath){
		try {
			this.dirpath = dirpath;
			IndexReader reader =readLuceneIndex(dirpath);
			//=============================================================
			long time  = System.currentTimeMillis();
			System.out.println("matrix is being created");
			createDocTermMatrix(reader);
			reader.close();
			reader=null;
			System.out.println(" matrix created in " +(System.currentTimeMillis()-time));
			//===================================================================
			
			time = System.currentTimeMillis();
			System.out.println("No of terms" +matrixStruct.columns());
			System.out.println("No of docs "+matrixStruct.rows());
			/** Perform LSI and  extract the topic by doc matrix**/
			double[][] result = extractSemanticDocMatrix();
			System.out.println(" LSI completed in  " +(System.currentTimeMillis()-time) );
			System.out.println("no of topics "+result.length+ " no of docs"+ result[0].length);
			
			//========================================================================
			/** order the topics from most significant to least significant**/
			result =organisetopics(result);
			System.gc();
			System.out.println("no of topics "+result.length+ " no of docs"+ result[0].length);
			
			//=======================================================================================
			
			time = System.currentTimeMillis();
			/* this is equivalent to (D^(-1)*U^T *A) where A is doc by term matrix A = U D V^T*/
			double[][]topicTermMatrix =MatrixOperations.denseMatrixSparseMatrixMult(result, matrixStruct,null);
			System.out.println("no of topics "+topicTermMatrix.length+ " no of terms"+ topicTermMatrix[0].length);
			int numofTopics = topicTermMatrix.length;
			for (int i = 0; i < numofTopics; i++) topicTermMatrix[i][i]= topicTermMatrix[i][i]/singularValues[i];
			System.out.println("Topics by term matrix generated" +(System.currentTimeMillis()-time));
			
			//==================================================================================================
			//store the two matrices
			storeMatrix(dirpath, topicTermMatrix,"topicTermMatrix");
			/** transpose the topic doc matrix to doc topic and store the result**/
			result = MatrixOperations.transposeMatrix(result);
			storeMatrix(dirpath, result, "doctopicMatrix");
			storeMap(dirpath,matrixStruct.returnDictionary(),"dictionary",null,null);
			storeMap(dirpath,null,null,matrixStruct.returnidByName(),"idByNames");
			//printTopics(topicTermMatrix);
			topicTermMatrix =null;
			System.gc();
			
			System.out.println("success");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
  
  public static IndexReader readLuceneIndex(String dirpath) throws IOException, InterruptedException{
    	Directory dir = FSDirectory.open(new File(dirpath));
    	System.out.println(new File(dirpath).exists());
		if (!IndexReader.indexExists(dir)){
			throw new IOException();
			}	
		return IndexReader.open(dir);
	}
	
	/**This method integrates lucene with colt
	 * Extracts the number of distinct terms and create a document term matrix by:
	 * Step 1: Iterating over each document and extracting its list of terms and their term-doc frequencies.
	 * Step 2: Compute the tf-idf score associated with each term in a document
	 */
	public  void createDocTermMatrix(IndexReader reader) throws IOException, InterruptedException{
		int docNum= reader.maxDoc();
		Similarity measure =Similarity.getDefault();
		
		//filter out deleted docs from lucene index
		ArrayList<Integer> listOfDocs = new ArrayList<Integer>();
		for(int i=0;i<docNum;i++){
			if (!reader.isDeleted(i)) {
				listOfDocs.add(i);
			}
		}
	    //initialise matrix
		int numofdocs = listOfDocs.size();
		//================
		Integer[] listdocs= new Integer[numofdocs];
		listdocs=(Integer[])listOfDocs.toArray(listdocs);
		
		
		//create a doc by term matrix
		matrixStruct = new SparseMatrix(reader.terms(), numofdocs);
		
		for (int j=0;j<listOfDocs.size();j++){
    		try{
    			TermFreqVector v = reader.getTermFreqVector(listOfDocs.get(j), "contents");
    			if( v == null) continue;
    			int [] freq = v.getTermFrequencies();
    			String[] terms = v.getTerms();
    			int[] nonZeroColumns = new int[terms.length];
    			double[] values = new double[terms.length];
    			for (int k=0;k<terms.length;k++){
    				int numOfdocforterm = reader.docFreq(new Term("contents",terms[k]));
    				nonZeroColumns[k] = matrixStruct.returnTermID(terms[k]);
    				values[k] = measure.tf(freq[k])* measure.idf(numOfdocforterm, listOfDocs.size());
    			 }
    			matrixStruct.addRow(j, nonZeroColumns, values);	
    			
    			
    		}
    		catch (Exception e) {
 				System.out.println("doc: "+j+"");
 				e.printStackTrace();
 			}
    		
    	}
		
	}
	/** The method performs Singular Value decomposition of Document by Terms matrix and extracts  Document by Topic Matrix*/
	public  double[][] extractSemanticDocMatrix() throws IOException{
		EigenDecompositionSparse eigen = new EigenDecompositionSparse(matrixStruct,250);
		singularValues=eigen.getSingularValues();
		return (eigen.getV());
	}
	/** order then topics from most significant to least significant*/
	private double[][] organisetopics(double[][] matrix) {
		int [] index = new int[singularValues.length];
		//=================================================================================================
		for(int i=0;i<index.length;i++)index[i]=i;
		for(int i=0;i<singularValues.length;i++){
			int max=i;
			for(int j=i;j<singularValues.length;j++){
				if (singularValues[j]>singularValues[i]) max=j;
			}
			double tmp = singularValues[i];
			int tmpIndex = index[i];
			singularValues[i] = singularValues[max];
			index[i] = index[max];
			singularValues[max] = tmp;
			index[max] = tmpIndex;
			}
		double sum = 0;
		for (double value: singularValues) sum += Math.abs(value);
		int limit = 0;
		double currentSum = 0;
		for (double value:singularValues) {
			currentSum += Math.abs(value);
			if (currentSum/sum < 0.80) limit++;
			else break;
		}
		//====================================================================================================
		/**update the matrix**/
		double[][] truncated_matrix = new double[limit][];
		for (int i=0;i<limit;i++){
			truncated_matrix[i]= matrix[index[i]];
		}
		matrix = truncated_matrix;
		truncated_matrix = null;
		System.gc();
		
		double[] sig_values= new double[limit];
		for (int i=0;i<limit;i++){
			sig_values[i]= singularValues[i];
			
		}
		this.singularValues = sig_values;
		return matrix;
  }
	
	
	public static  void storeMatrix(String dirpath,Object matrix, String matrixtype) throws IOException{
		File file  = new File (dirpath+"/"+matrixtype+".ser");
		if (file.exists()){
			file.delete();
			file.createNewFile();
		}
		FileOutputStream fileStream= new FileOutputStream(file);
		ObjectOutputStream os = new ObjectOutputStream(fileStream);
		os.writeObject(matrix);
		os.close();
	}
	/** store the associations of each term with its corresponding column in the doc by the term matrix**/
	private void storeMap(String dirpath, HashMap<String,Integer> map,String type,HashMap<Integer,String> map2, String type2) throws IOException{
		FileOutputStream fileStream = null;
		//====================================================================
		if (type2==null){
			File file  = new File (dirpath+"/"+type+".ser");
			if (file.exists()){
				file.delete();
				file.createNewFile();
			}
			fileStream= new FileOutputStream(file);
		}
		//=============================================================================
		else {
			File file  = new File (dirpath+"/"+type2+".ser");
			if (file.exists()){
				file.delete();
				file.createNewFile();
			}
			fileStream = new FileOutputStream(file);
		}
		//==========================================================================
		ObjectOutputStream os = new ObjectOutputStream(fileStream);
		
		if(map2==null)os.writeObject(map);
		else os.writeObject(map2);
		os.close();
	}
	
	
	private void printTopics(double[][] topicTermMatrix) {
		for (int i=0 ; i<5;i++){
			double[] ts =topicTermMatrix[i];
			int[] index = new int[ts.length];
			for(int j=0;j<ts.length;j++){
				index[j]=j;
			}
			StringBuilder strbuil= new StringBuilder();
			for (int k=0;k<20;k++){
				int max=k;
				for(int l=k;l<ts.length;l++){
					if((ts[l])>(ts[max])) max=l;
				}
				double tmp;
				tmp=ts[k];
				ts[k]=ts[max];
				ts[max]=tmp;
				
				int temp;
				temp=index[k];
				index[k]=index[max];
				index[max]=temp;
				strbuil.append(matrixStruct.getTerm(index[k]));
				strbuil.append("*");
				strbuil.append(ts[k]);
				strbuil.append(" ");
				
			}
			System.out.println(strbuil.toString());
			System.out.println();
				
			}
	}
	
	
/** assign a AnalyseDoc object to each workspace ID*/	
	public static void main(String[] args) throws IOException {
	if (!args[0].equals("-lucenedirpath") ) throw new IOException("type -lucenedirpath followed by the dirpath");
	String dirpath =args[1];
	if (dirpath==null)throw new IOException(" the lucene dirpath must be specified");
	new AnalyseDocs(dirpath);
	
	}
	
	
}
