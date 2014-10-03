package co.annotate.LatentSemanticAnalysis;

import java.io.IOException;
import java.util.Arrays;


/**@author Adnan
* @function  A memory efficient  and computational fast way of extracting a doc-topic matrix from a doc-term matrix;
 * The method takes a rectangular matrix A(m by n) and extracts the eigenvectors of AA using the Arnold algorithm and the implicit QL decomposition using combination of jacobian and given rotations.
 * */
public class EigenDecompositionSparse implements java.io.Serializable {

    private static final long serialVersionUID = 1L;
    /* number of columns*/
    private int n;
    /* d and e: diagonal and off diagonal elements of the  hessenbergh matrix*/ 
    private double[] d;
	  private double[] e;
    /*V: eigen vectors of H*/
    private double[][] V;
    private double[][] Q;
    private int iter;
    
    public EigenDecompositionSparse(SparseMatrix A,int num) throws IOException{
         n= A.rows();
         iter = (n>num)? num :n;
         V = new double[iter][iter];
         for(int i =0;i<iter;i++)V[i][i]=1;
         d = new double[iter];
         e = new double [iter];
         Q = new double[iter][n];
        /* compute H where H is hessenbergh and H = Q'AQ*/
         ArnoldiAlg(A);
         /* implicit QL decomposition using combination of jacobian and given rotations*/
         tql2();
         
         
    }
    
    public EigenDecompositionSparse(SparseMatrix A) throws IOException{
    	this(A,A.columns()>250?250:A.columns());
    }

    private void ArnoldiAlg(SparseMatrix A) throws IOException {
        /*initialize q1*/
        double[] q =  new double[n];
        MatrixOperations.assign(q,1.0);
        double norm_q =(float)Math.sqrt(MatrixOperations.dotproductDense(q, q));
        MatrixOperations.multiplyVecByScalar(q, 1/norm_q);
       
        //breaking the computation of AA'q to A'q followed by A (A'q)
        double[] r =  A.matrixTransposeVecMult(q);
        r= A.matrixVecMult(r);		
        double aj =MatrixOperations.dotproductDense(q,r);
        MatrixOperations.vectorSub(r, q, aj);
        double bj = (float)Math.sqrt(MatrixOperations.dotproductDense(r, r));
       
        
        d[0]=aj;
        e[0]=bj;
        
        updateQ(0,q);
        
        for (int i=1;i<iter;i++){
        	double[] v =q;
        	MatrixOperations.multiplyVecByScalar(r, (double)1/bj);
    		q=r;
    		System.gc();
    		updateQ(i, q);
    		/* computing Aqi- bjv followed by Aqi-ajqi*/
    		r =  A.matrixTransposeVecMult(q);
    		r= A.matrixVecMult(r);	
    		System.gc();
            
    		MatrixOperations.vectorSub(r, v, bj);
    		aj =MatrixOperations.dotproductDense(q, r);
    		d[i]=aj;
    		MatrixOperations.vectorSub(r, q, aj);
            
    		/* re-orthogonalisation*/
    		if(i>1){
    			double[]r_p = MatrixOperations.denseMatrixVecMult(Q, r);
    			r_p = MatrixOperations.denseMatrixTransposeVecMult(Q, r_p);
    			MatrixOperations.vectorSub(r, r_p, 1.0);
    		}
    		bj = (float)Math.sqrt(MatrixOperations.dotproductDense(r, r));
    		if (i!=(iter-1))e[i]=bj;
    		if ((double)(int)(bj*100000000)==0) break;
        }
        e[iter-1]=0;
    }
   
    
    
    private void updateQ(int index,double[] q){
        Q[index]=q;
       }
    
    private void tql2 () {

        //  This is derived from the Algol procedures tql2, by
        //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.
       
          n=iter;
          //e[iter-1] = (float) 0.0;
          double f = 0;
          double tst1 = 0;
          double eps = Math.pow(2.0,-52.0);
          for (int l = 0; l < n; l++) {
             // Find small sub-diagonal element
             tst1 = Math.max(tst1,Math.abs(d[l]) + Math.abs(e[l]));
             int m = l;
             while (m < n) {
                if (Math.abs(e[m]) <= eps*tst1) {
                   break;
                }
                m++;
             }
             // If m == l, d[l] is an eigenvalue,
             // otherwise, iterate.
             
             if (m > l) {
                int iter = 0;
                do {
                   iter = iter + 1;  // (Could check iteration count here.)
                   // Compute implicit shift through Jacobian rotation
                   double g = d[l];
                   /* p is theta/cot2phi*/
                   double p =  ((d[l+1] - g) / (2.0 * e[l]));
                   
                   /*t = 1/sign(theta)+ sqrt(theta^2+1))*/
                   double r = hypot(p,1.0);
                   if (p < 0) {
                      r = -r;
                   }
                   /* d[l]= e[l]*t*/
                   d[l] =  (e[l] / (p + r));
                   d[l+1] =  (e[l] * (p + r));
                   double dl1 = d[l+1];
                   double h = g - d[l];
                   for (int i = l+2; i < n; i++) {
                      d[i] -= h;//d[i]= d[i]-ks
                   }
                   /* f at l=0 is equal to k_s and is the sum of Sum ks for all l*/
                   
                   f = f + h;
       
                   // Implicit QL transformation.
                  
                   p = d[m];
                   
                   double c = 1.0;
                   double c2 = c;
                   double c3 = c;
                   double el1 = e[l+1];
                   double s = 0.0;
                   double s2 = 0.0;
                   /* 1 Jacobian rotation followed by Given rotations*/
                   for (int i = m-1; i >= l; i--) {
                      c3 = c2;
                      c2 = c;
                      s2 = s;
                      /* for i=m-1, c2,c3,c=1 and s=0*/
                      g =  (c * e[i]);
                      h =  (c * p);
                      /* r  = sqrt( e[m-1]^2 + d[m]^2) for i =m-1*/
                      r =  hypot(p,e[i]);
                      /* for i=m-1 e[i+1]= 0  */
                      e[i+1] =  (s * r);
                      s = e[i] / r;
                      c = p / r;

                      p =  (c * d[i] - s * g);
                      d[i+1] =  (h + s * (c * g + s * d[i]));
                        // Accumulate transformation.
                        for (int k = 0; k < n; k++) {
                        	 h = V[k][i+1];
        					 V[k][i+1] = s * V[k][i] + c * h;
        					 V[k][i] = c * V[k][i] - s * h;
                      }
                   }
                   
                   p =  (-s * s2 * c3 * el1 * e[l] / dl1);
                   e[l] =  (s * p);
                   d[l] =  (c * p);
                   // Check for convergence.
       
                } while (Math.abs(e[l]) > eps*tst1);
             }
             d[l] =  (d[l] + f);
             e[l] =  0.0;
          }
        // Sort eigenvalues and corresponding vectors.
          for (int i=0;i<d.length;i++){
              int max= i;
              for(int k =i;k<d.length;k++){
                  if (d[k]>d[max]) max =k;
              }
              double tmp ;
              tmp=d[i];
              d[i]= d[max];
              d[max] =tmp;
              for (int j=0;j<n;j++){
                  tmp =  V[j][i];
                  V[j][i]=V[j][max];
                  V[j][max]=tmp;
              }
          }
          
          
         
    } 
    /**
    * Returns sqrt(a^2 + b^2) without under/overflow.
    */
   protected static double hypot(double a, double b) {
   	double r;
   	if (Math.abs(a) > Math.abs(b)) {
   		r = b/a;
   		r = Math.abs(a)*Math.sqrt(1+r*r);
   	} else if (b != 0) {
   		r = a/b;
   		r = Math.abs(b)*Math.sqrt(1+r*r);
   	} else {
   		r = 0.0;
   	}
   	return r;
   }
    
    public double[][] getV() throws IOException{
    	/* here Q is iter by doc matrix so we need to multiply Q'V where v is iter by iter matrix to get doc by iter matrix  
    	 * but for efficiency we get iter by doc matrix*/
    	double[][]eigenvectors = new double[iter][n];
        for (int i=0;i<iter;i++){
        	double[] column =new double[V[0].length];
        	//extract column from V
        	for (int k=0;k<iter;k++) column[k] = V[k][i];
        	eigenvectors[i] = MatrixOperations.denseMatrixTransposeVecMult(Q,column);
        }
       // System.out.println(new DenseDoubleMatrix2D(eigenvectors));
		return eigenvectors;
    }
    public double[] getSingularValues(){
    	//System.out.println(new DenseDoubleMatrix1D(d));
        return d;
        
    }
    
   
	
    
}
