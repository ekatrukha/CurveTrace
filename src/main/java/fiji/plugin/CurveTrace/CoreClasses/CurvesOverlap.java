package fiji.plugin.CurveTrace.CoreClasses;

import java.util.ArrayList;
import java.util.List;

public class CurvesOverlap {

	/** overlap types 0 = FULL, 1 = EXTENSION, 2 = BRANCHING, 3 = CROSS, 4= NONE*/
	public static final int FULL=0, EXTENSION=1, BRANCHING=2, CROSS=3, NONE=4;

	/** direction of overlap with respect to curve1 numbering of points
	 *  0 = SAME, 1 = OPPOSITE*/
	public static final int SAME=0, OPPOSITE=1;

	/** indexes of overlapping points in Curve array**/
	public ArrayList<int []> indexes;
	/** averaged points**/
	public ArrayList<Point> averaged_points;
	/** type of overlap **/
	public int nOvelapType;
	/** direction of overlap (points numbering of curve2 with respect to curve1 **/
	public int nOverlapDirection;
	/** fraction of overlap for curve1**/
	public double overlapFraction1;
	/** fraction of overlap for curve2**/
	public double overlapFraction2;
	/** number of points in the overlap **/
	public int nOverlapSize;
	/** array containing min and max indexes of overlap for curves 
	 * [i][j] -> i = curve number, j = min (0) or max (1) **/
	public int [][] indminmax;
	
	/**indexes of curves**/
	public int [] curves_ind;
	
	/** if one of curve changes direction of numbering at overlap **/
	public boolean bDirectConflict=false;	
	/** array containing indexes of overlap where one of curves 
	 * changes direction of numbering points **/
	public List<Integer> dirChange ;
	
	/** maximum curve index jump during overlay  (any of two curves) **/
	public double dMaxJumpFraction;
	
	
	public CurvesOverlap()
	{
		indminmax= new int [2][2];
		indexes = new ArrayList<int[]>();
		averaged_points = new ArrayList<Point>();
		curves_ind = new int[2];
		dirChange = new ArrayList<>();
	}
	
	public void calculateMinMaxIndexes()
	{
		int i,j;
		int [] ind = new int [2];
		
		//initialization
		for (j=0;j<2;j++)
		{
			indminmax[j][0]=Integer.MAX_VALUE;
			indminmax[j][1]=Integer.MIN_VALUE;
		}

		for (i=0;i<nOverlapSize;i++)
		{
			ind=this.indexes.get(i);
			for (j=0;j<2;j++)
			{
				
				if(ind[j]<indminmax[j][0])
					indminmax[j][0]=ind[j];
				
				if(ind[j]>indminmax[j][1])
					indminmax[j][1]=ind[j];
			
			}	

		}
		return;
		
	}
	
	public void calculateOverlapDirection()
	{
		if(indexes.get(0)[1]<indexes.get(nOverlapSize-1)[1])
		{
			nOverlapDirection = SAME;
		}
		else
		{
			nOverlapDirection = OPPOSITE;
		}
		
	}
	/**  Method analyzes overlap consistency:
	 *  1) if one of curves changes numbering direction
	 *  during overlap, there is multiple crossings between lines,
	 *  so we need to extract biggest one
	 *  2) if there is a big jump in points along the curve inside overlap,
	 *  we should record it and later eliminate this kind of overlap **/
	public void analyzeOverlap()
	{
		// ANALYZE OVERLAP
		int i,j;
		
		//analyzing sequence of points
		int [] dC = new int[2];
		int [] dCurr = new int[2];
		boolean [] bInit= new boolean[2];
		bInit[0]=false;
		bInit[1]=false;

		double dMaxFraction = 0.0;
		double dCurrFraction=0.0;
		int nMaxInd = 0;
		
		//just in case
		nOverlapSize=this.indexes.size();
		dirChange = new ArrayList<>();
		dirChange.add(0);
		for(i=1;i<nOverlapSize;i++)
		{
			for (j=0;j<2;j++)
			{
				dCurr[j]=indexes.get(i)[j]-indexes.get(i-1)[j];
				if(Math.abs(dCurr[j])>0 && !bInit[j])
				{
					dC[j]=dCurr[j];
					bInit[j]=true;
				}
				else
				{
					if(Math.abs(dCurr[j])>0)
						if(Math.signum((double)dC[j])!=Math.signum((double)dCurr[j]))
						{
							bDirectConflict=true;
							dirChange.add(i);
							dC[j]=dCurr[j];
						}
				}

			}//end of for (j=0;j<2;j++)
			
	
		}//end of for(i=1;i<nOverlapSize;i++)
		
		//ok, seems like there is a change in the numbering
		if(bDirectConflict)
		{
			//finding maximum stretch where direction is not changing
			dirChange.add(nOverlapSize);
			for (i=1;i<dirChange.size();i++)
			{
				//indexes.
				dCurrFraction = (double)(dirChange.get(i)-dirChange.get(i-1))/(double)nOverlapSize;
				if(dCurrFraction>dMaxFraction)
				{
					dMaxFraction=dCurrFraction;
					nMaxInd=i;
				}
			}
			
			//got it 
			//not let's shrink overlay arrays to this value
			ArrayList<int []> indexes_new = new ArrayList<int []>();
			ArrayList<Point> points_new = new ArrayList<Point>();
			int [] new_ind; 
			for (i=dirChange.get(nMaxInd-1);i<dirChange.get(nMaxInd);i++)
			{
				new_ind = new int [2];
				new_ind[0]=indexes.get(i)[0];
				new_ind[1]=indexes.get(i)[1];
				indexes_new.add(new_ind);
				points_new.add(averaged_points.get(i));
			}
			
			indexes=indexes_new;
			averaged_points=points_new;
			nOverlapSize=indexes.size();
		}
		//calculate maximum jump fraction
		calculateMaxJumpFraction();
	}
	
	/** Function assessing how continuous overlap is,
	 * i.e. are points of both curves staying close to each other
	 * by looking at the indices of overlap.
	 *  If there is a big jump in points indices along the curve inside overlap,
	 *  we should record it in dMaxJumpFraction and later in code eliminate this kind of overlap **/
	public void calculateMaxJumpFraction()
	{
		int i,j;
		int nStepDiff=0;
		double nMaxStep=0;
		int [] dCurr = new int[2];

		for(i=1;i<nOverlapSize;i++)
		{
			for (j=0;j<2;j++)
			{
				dCurr[j]=indexes.get(i)[j]-indexes.get(i-1)[j];
			}
			nMaxStep=Math.max(Math.abs((double)dCurr[0]), Math.abs((double)dCurr[1]));
			if((int)nMaxStep>nStepDiff)
			{
				nStepDiff=(int)nMaxStep;
			}
		}
		dMaxJumpFraction =nStepDiff/(double)nOverlapSize;
		
	}
}
