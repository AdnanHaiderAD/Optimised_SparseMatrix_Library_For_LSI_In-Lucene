package co.annotate.LatentSemanticAnalysis;

import java.io.IOException;


/*Author: Adnan
 * @description : library functions for algebraic operations between dense  matrix/vectors or both.
 * */
public class MatrixOperations {

	/** assign the same value to all elements in a vector*/
    public static void assign(double[] q2, double f) {
	   for (int i=0; i< q2.length;i++) q2[i]=f;
	}
    /** dot product of two dense arrays*/
   public static double dotproductDense(double[] q1,double[] q2) throws IOException{
    	if (q2.length!=q1.length)throw new IOException("the two vectors must be of the same size");
    	double value=0;
    	for (int i=0; i< q2.length;i++) value+=q2[i]*q1[i];
    	return value;
    }
    public static void multiplyVecByScalar(double[]q,double s){
    	for (int i=0; i< q.length;i++) q[i]=q[i]*s;
    }
    /** compute vector = q-aq1*/
    public static void vectorSub(double[]q,double[]q1,double a) throws IOException{
    	if (q.length!=q1.length)throw new IOException("the two vecotrs must be of the same size");
    	for (int i =0; i<q.length;i++)q[i]= q[i]-a*q1[i]; 
    }
    
    
    /** multiply a dense matrix with a dense vector*/
   public  static double[] denseMatrixVecMult(double[][] matrix, double[] vec) throws IOException{
    	if(matrix[0].length!=vec.length)throw new IOException("the number of columns of the matrix doesnt match "
    			+ "the number of rows of the vector");
    	double[] result = new double[matrix.length];
    	for (int i=0; i< matrix.length;i++){
    		double[] row = matrix[i];
    		for (int k =row.length/10*10+row.length%10;--k>=row.length/10*10;){
    			result[i]+=row[k]*vec[k];
    		}
    		for (int k= row.length/10;--k>=0;){
    			result[i]+= row[k*10+9]*vec[k*10+9]+
    					row[k*10+8]*vec[k*10+8]+
    					row[k*10+7]*vec[k*10+7]+
    					row[k*10+6]*vec[k*10+6]+
    					row[k*10+5]*vec[k*10+5]+
    					row[k*10+4]*vec[k*10+4]+
    					row[k*10+3]*vec[k*10+3]+
    					row[k*10+2]*vec[k*10+2]+
    					row[k*10+1]*vec[k*10+1]+
    					row[k*10+0]*vec[k*10+0];
    					
    		}
    	}
    	return result;
    }
    /** multiply the transpose of a dense matrix with a vector*/
    public static double[] denseMatrixTransposeVecMult(double[][] matrix, double[] vec)throws IOException{
    	if(matrix.length!=vec.length)throw new IOException("the number of columns of the transpose matrix doesnt match "
    			+ "the number of rows of the vector");
    	double[] result = new double[matrix[0].length];
    	//System.out.println(result.length);
    	
    	for (int i=0; i< result.length;i++){
    		double[] column = new double[matrix.length];
    		for(int j=0;j<column.length;j++) column[j]= matrix[j][i];
    		for (int k =column.length/10*10+column.length%10;--k>=column.length/10*10;){
    			result[i]+= column[k]*vec[k];
    		}
    		for (int k= column.length/10;--k>=0;){
    			result[i]+= column[k*10+9]*vec[k*10+9]+
    					column[k*10+8]*vec[k*10+8]+
    					column[k*10+7]*vec[k*10+7]+
    					column[k*10+6]*vec[k*10+6]+
    					column[k*10+5]*vec[k*10+5]+
    					column[k*10+4]*vec[k*10+4]+
    					column[k*10+3]*vec[k*10+3]+
    					column[k*10+2]*vec[k*10+2]+
    					column[k*10+1]*vec[k*10+1]+
    					column[k*10+0]*vec[k*10+0];
    		}
    	}
    	return result;
    }
    /**computes the matrix vector multiplication of a dense matrix A with a sparse matrix V i.e A*V. 
     * m denotes the number of rows of A that take part in the computation*/
    public static double[][] denseMatrixSparseMatrixMult(double[][] A , SparseMatrix V, Integer m) throws IOException{
    	if  (V.rows()!=A[0].length) throw new IOException("the number of columns of A is not equal to the number of rows of V");
    	if (m==null) m = A.length;
    	double[][] result = new double[m][V.columns()];
    	for (int i=0;i<m;i++){
    		double[] A_row = A[i];
    		for (int j =0;j<A_row.length;j++){
    			double value = A_row[j];
    			SparseVec v = V.getRow(j);
    			while(v.getCurrentCol()!=-1){
    				result[i][v.getCurrentCol()]= result[i][v.getCurrentCol()]+value* v.currentValue();
    				v.next();
    			}
    			v.reset();
    		}
    		
    	}
    	return result;
    	
    }
    /**multiplication of the transpose of a dense matrix with a sparse matrix**/
    public static double[][] denseMatrixTransposeSparseMatrixMult(double[][] A , SparseMatrix V) throws IOException{
    	if  (V.rows()!=A.length) throw new IOException("the number of columns of the transpose of A is not equal to the number of rows of V");
    	/**construct the transpose of the dense matrix**/
    	double[][] matrix = new double[A[0].length][A.length];
    	for (int i=0;i<matrix.length;i++){
    		for (int j=0;j<matrix[0].length;j++){
    			matrix[i][j]=A[j][i];
    		}
    	}
    	A=null;
    	System.gc();
    	
    	double[][] result = new double[matrix.length][V.columns()];
    	for (int i=0;i<matrix.length;i++){
    		double[] A_row = matrix[i];
    		for (int j =0;j<A_row.length;j++){
    			double value = A_row[j];
    			SparseVec v = V.getRow(j);
    			while(v.getCurrentCol()!=-1){
    				result[i][v.getCurrentCol()]= result[i][v.getCurrentCol()]+value* v.currentValue();
    				v.next();
    			}
    			v.reset();
    		}
    		
    	}
    	return result;
    	
    }
    
    
    /** multiplication of a dense matrix with a sparse vector*/
    public static double[] denseMatrixSparseVecMult(double[][]A, SparseVec V) throws IOException{
    	if (A[0].length!= V.getDimension())throw new IOException("the number of columns of A doesnt match the number"
    			+ " of rows of the Sparse vector V");
    	double[] result= new double[A.length];
    	for (int i=0;i<A.length;i++){
    		double[] A_row=A[i];
    		while(V.getCurrentCol()!=-1){
    			result[i]= result[i]+ V.currentValue()*A_row[V.getCurrentCol()];
    			V.next();
    			
    		}
    		V.reset();
    	}
    	return result;
    }
    /** compute the cosine similarity matrix of a given doc by topic matrix*/
     public static double[][]angularDistanceMatrix(double[][] docTopicMatrix) throws IOException{
     	double[][] similarityMatrix = new double[docTopicMatrix.length][];
     	for (int i=0;i<similarityMatrix.length;i++){
     		double[] scores = new double[similarityMatrix.length - i];
     		for ( int j=i+1; j<similarityMatrix.length;j++){
     			double dotprod=0;
     			double dotprodA=0;
     			double dotprodB=0;
     			double[] vec= docTopicMatrix[i];
     			double[] vec2= docTopicMatrix[j];
     			for (int b=0;b<vec.length;b++){
     				dotprod= dotprod+vec[b]*vec2[b];
     				dotprodA =dotprodA+vec[b]*vec[b] ;
     				dotprodB =dotprodB+vec2[b]*vec2[b];
     			}
     			double d =1-( Math.acos(dotprod/(Math.sqrt(dotprodA*dotprodB)))/Math.PI);
     	        d = roundToSignificantFigures(d,4);
     	        scores[j-i]=d;
     	          // similarityMatrix[i][j]=d;
     	          //similarityMatrix[j][i]=d;
     	         
            }
     		similarityMatrix[i]=scores;
     		
     	}
        return similarityMatrix;
     }
     /** computes the angular distance between a query string and docs*/
     public static double[] angularDistance(double[][] docTopicMatrix,double[] vector){
    	 double[] result = new double[docTopicMatrix.length];
    	 for (int i=0;i<result.length;i++){
    		 double[] vec= docTopicMatrix[i];
    		 double dotprod=0;
  			 double dotprodA=0;
  			 double dotprodB=0;
  			 for (int b=0;b<vec.length;b++){
  				dotprod= dotprod+vec[b]*vector[b];
  				dotprodA =dotprodA+vec[b]*vec[b] ;
  				dotprodB =dotprodB+vector[b]*vector[b];
  			}
  			double d =1-( Math.acos(dotprod/(Math.sqrt(dotprodA*dotprodB)))/Math.PI);
 	        d = roundToSignificantFigures(d,4);
 	        result[i]=d;
  		}
    	return result;
     }
     /** computes the euclidean distance between different docs in a doc by topic matrix**/
     public static double[][] computeEuclideandistanceMatrix(double[][] docTopicMatrix){
    	 double[][] result = new double[docTopicMatrix.length][docTopicMatrix.length];
    	 for (int i =0 ; i < docTopicMatrix.length; i++){
    		double[] vec = docTopicMatrix[i];
    		result[i] = euclideanDistance(docTopicMatrix, vec); 
    	 }
    	 return result;
     }
     
     /** computes the euclidean distance between a query string and docs*/
     public static double[] euclideanDistance (double[][] docTopicMatrix,double[] vector){
    	 double[] result = new double[docTopicMatrix.length];
    	 for (int i=0;i<result.length;i++){
    		 double[] vec= docTopicMatrix[i];
    		 double value=0;
    		 for (int b=0;b<vec.length;b++){
    			 value= value + (vec[b]-vector[b])*(vec[b]-vector[b]);
    		 }
    		 result[i]=Math.sqrt(value);
    	 }
    	 return result;
     }
     static double[] normalise(double[] vec){
    	 double mod_vec=0;
    	 for (int i=0;i<vec.length;i++){
    		 mod_vec= mod_vec+ vec[i]*vec[i];
    	 }
    	 mod_vec = Math.sqrt(mod_vec);
    	 for (int i=0;i<vec.length;i++){
    		 vec[i] = (vec[i]/mod_vec);
    	 }
    	 return vec;
     }
     
     static double[] convertVecttoProbility(double[] vector){
    	 double sum = 0 ;
    	 for (double item : vector)sum+= item;
    	 for (int i =0; i < vector.length;i++ )vector[i] = roundToSignificantFigures(vector[i]/sum,3);
    	 return vector;
     }
     
     public static double roundToSignificantFigures(double d, int i) {
     	if (d==0) return d;
     	double num = Math.ceil(Math.log10(d<0?-d:d));
     	int power = i -(int)num;
     	double magnitude =Math.pow(10, power);
     	long shifted =Math.round(d*magnitude);
     	return shifted/magnitude;
 	}
	public static double[][] transposeMatrix(double[][] result) {
		double[][] transpose = new double[result[0].length][result.length];
		for (int i =0;i<result.length;i++){
			for (int j=0;j<result[0].length;j++){
				transpose[j][i]=result[i][j];
			}
		}
		return transpose;
	}
     
     
	
}
