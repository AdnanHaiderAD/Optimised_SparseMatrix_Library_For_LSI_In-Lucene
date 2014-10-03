package co.annotate.LatentSemanticAnalysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.Term;
import org.apache.lucene.index.TermEnum;





/**@author  Adnan
 * 
 * @function provides an efficient way to store and perform algebraic operations on BIG sparse matrices in memory 
 * */
public class SparseMatrix {
  /** dictionary encapsulates the mapping between normalised and their integer ids*/
	private HashMap<String,Integer>dictionary;
	private HashMap<Integer,String> idTermRefs;
	private SparseVec[] rows; 
	private int nt;
	
	
	
	/** construct a memory efficient sparse matrix from the lucene index*/
	public SparseMatrix(TermEnum enumer,int docNum) throws IOException{
		this.dictionary = new HashMap<String, Integer>();
		this.idTermRefs= new HashMap<Integer, String>();
		this.rows = new SparseVec[docNum];
		 nt = 0;
		 while(enumer.next()){
			 Term t = enumer.term();
			 if (t.field().equals("contents") ){
				insertTerm(t.text(), nt);
				 nt=nt+1;
			 }
		}
	}
	/**default constructor*/
	public SparseMatrix(int rows, int columns){
		this.rows = new SparseVec[rows];
		this.nt =columns;
	}

	private void insertTerm(String text, Integer nt) {
		this.dictionary.put(text, nt);
		this.idTermRefs.put(nt, text);
	}
	public void addRow(int row,double[] vec){
		this.rows[row] =new SparseVec(vec);
	}
	/** add a row vector using information of only non -zero column elements*/
	public void addRow(int row, int[] columns, double[] values){
		this.rows[row]= new SparseVec(columns, values, nt+1);
	}
	public int returnTermID(String term){
		return dictionary.get(term);
	}
	
	/** sparse matrix-sparse vector multiplication*/
	public double[] matrixVecMult(SparseVec vec) throws IOException{
		int vecDim = vec.returnDimension();
		if (vecDim!=nt)throw new IOException("the number of columns of the matrix doesnt match"
				+ " the number of rows of the vector");
		double[] output = new double[rows.length];
		for (int i=0;i<output.length;i++){
			output[i]= rows[i].innerProd(vec,nt);
		}
		return output;
	}
	/** sparse matrix dense vector multiplication*/
	public double[] matrixVecMult(double[] vec) throws IOException{
		if (nt!=vec.length) throw new IOException("the number of columns of the matrix doesnt match "
				+ "the number of rows of the vector");
		double[] result = new double[rows.length];
		
		for (int j=0;j<rows.length;j++){
			SparseVec row = rows[j];
			int currentCol =row.getCurrentCol();
			while(currentCol!=-1){
				result[j]= result[j]+ vec[currentCol]*row.currentValue();
				row.next();
				currentCol= row.getCurrentCol();
			}
			row.reset();
		}
		return result;
	}
	
	/** multiplication of the transpose of the  matrix with a vector i.e  A' *v*/
	public double[] matrixTransposeVecMult(double[] vec) throws IOException{
		if (rows.length!=vec.length) throw new IOException("the number of columns of the matrix"
				+ " doesnt match the number of rows of the vector");
		double[] result = new double[nt];
		for (int j=0;j<vec.length;j++){
			double value = vec[j];
			SparseVec column = rows[j];
			for (int i=0;i<nt;i++){
				int nonZeroIndex = column.getCurrentCol();
				if (nonZeroIndex==-1) break;
				if (nonZeroIndex==i) {
					result[i]= result[i]+ value*column.currentValue();
					column.next();
				}else continue;
				
			}
			column.reset();
		}
		return result;
	}
	
	/** matrix multiplied by its transpose AA**/
	public double[][] matrixMultTransposeR(){
		double[][] result = new double[rows.length][];
		//store only the diagonal and off diagonal elements
		for(int i=0;i<rows.length;i++) result[i]=new double[i+1];
		
		for(int i=0;i<rows.length;i++){
			result[i][i] = rows[i].selfInnerProd();
		}
		
		for (int i=0;i<rows.length;i++){
			for (int j=i+1;j<rows.length;j++){
				result[j][i]= rows[i].innerProd(rows[j], nt);
				}
		}
		return result;
	}
	
	public HashMap<String, Integer> returnDictionary(){return dictionary;}
	public HashMap< Integer,String> returnidByName(){return this.idTermRefs;}

	public int columns() {
		// TODO Auto-generated method stub
		return nt;
	}
	public int rows() {
		// TODO Auto-generated method stub
		return rows.length;
	}
	public SparseVec getRow(int j){ return rows[j];}
	public String getTerm(int key) {return this.idTermRefs.get(key) ;
	}
	
}
