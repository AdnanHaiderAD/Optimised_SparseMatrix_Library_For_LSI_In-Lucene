package co.annotate.LatentSemanticAnalysis;

import java.util.ArrayList;
/**
 * @author : Adnan
 *  @function this class  provides a memory efficient way to store Sparse vectors and perform fast algebraic operations on them. */

public class SparseVec {
	private int[] columns;
	private double[] values;
	private int index;
	private int dimension;
	
	
	public SparseVec(int[] columns,double[] values,int dimension){
		this.values=values;
		this.columns=columns;
		this.index=0;
		this.dimension=dimension;
	}
	
	public int returnDimension(){ return this.dimension;}
	
	public SparseVec(double[] rowVec){
		ArrayList<Double> vals= new ArrayList<Double>();
		ArrayList<Integer> cols = new ArrayList<Integer>();
		for(int i =0; i<rowVec.length;i++){
			if (rowVec[i]!=0){
				vals.add(rowVec[i]);
				cols.add(i);
			}
		}
		if (cols.size()==0){
			this.values=new double[1];
			this.columns=new int[1];
			this.values[0]=0;
			this.columns[0]=0;
		}else{
			this.values =new double[cols.size()];
			this.columns = new int[cols.size()];
			for (int i=0;i<cols.size();i++){
				this.values[i]=vals.get(i);
				this.columns[i]=cols.get(i);
			}
			
		}
		
		vals=null;
		cols=null;
		System.gc();
		this.dimension= rowVec.length;
		
	}
	
	public double currentValue(){
		return values[index];
	}
	
	/* this method allows the class to behave like a iterator*/
	public void next(){
		index+=1;
		if (index==columns.length){
			index=-1;
			}
		}
	public void reset(){ this.index=0;}
	
	public int getCurrentCol(){return index==-1?-1:columns[index];}
	public int getDimension(){return dimension;}
	
	/** the dot product of the sparse vector with itself*/
	public double selfInnerProd(){
		double result=0;
		for (double value : values) result+=value*value;
		return result;
	}
	/**computes  the inner product ofthe sparse vector with the sparse vector v2**/
	public double innerProd( SparseVec v2, int nt){
		double value=0;
		for( int i=0; i<nt;i++){
			int v1_col = getCurrentCol();
			int v2_col=v2.getCurrentCol();
			//System.out.println("v1 "+ v1_col+ "v2 "+v2_col);
			if(v1_col ==-1|| v2_col==-1)break;
			if (v1_col>i){
				i=v1_col-1;
				if (v2_col<v1_col){
					while(v2_col<v1_col){
						v2.next();
						v2_col= v2.getCurrentCol();
						//System.out.println( "v2 "+v2_col);
						if (v2_col==-1)  break;
					}
				}
				continue;
			}
			double prod_value = currentValue();
			
			if(v2_col>i){
				i=v1_col-1;
				if (v1_col<v2_col){
					while(v1_col<v2_col){
						next();
						v1_col= getCurrentCol();
						//System.out.println( "v1 "+v1_col);
						if (v1_col==-1)break;
					}
				}
				continue;
			}
			value= value+ prod_value*v2.currentValue();
			next();
			v2.next();
			
		}
		index=0;
		v2.index=0;
		return value;
	}
	
}
