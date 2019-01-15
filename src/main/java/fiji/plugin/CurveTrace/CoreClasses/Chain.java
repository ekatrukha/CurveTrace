package fiji.plugin.CurveTrace.CoreClasses;

import java.util.ArrayList;

/** contains chains of curves connectivity,
 * in one int array, where
 * 0 - index of curve from set 1 
 * 1 - index of curve from set 2 
 * 2 - overlap index from alloverlaps array 
 * 3 - inverse flag of whether this element was processed is Breadth-first search (BFS) 0=processed 1=not processed **/

public class Chain extends ArrayList<int []>{
	

	public Chain()
	{
		super();
	}
	public Chain(int i, int j, int overlay_index)
	{
		
		super();
		int [] pair = new int [4];
		pair[0] = i;
		pair[1] = j;
		pair[2] = overlay_index;
		pair[3] = 1;
		this.add(pair);
		
	}
	public void addPair(int i, int j, int overlay_index)
	{
		int [] pair = new int [4];
		pair[0] = i;
		pair[1] = j;
		pair[2] = overlay_index;
		pair[3] = 1;
		this.add(pair);
	}
	public boolean isProcessed()
	{
		//boolean bOut=true;
		int i;
		
		for (i=0;i<this.size();i++)
		{
			if(this.get(i)[3]>0)
			{
				return false; 
			}
		}
		return true;
		
	}
}
