package fiji.plugin.CurveTrace.CoreClasses;

import java.util.ArrayList;
/**
 *	Curve class, contains array of points
 */
public class Curve extends ArrayList<Point> {

	public Curve()
	{
		super();
	}
	
	public Curve duplicate()
	{
		Curve finCurve = new Curve();
		int i;
		
		for (i=0;i<this.size();i++)
		{
			finCurve.add(this.get(i).clone());
		}
		
		
		return finCurve;
	}
	
}
