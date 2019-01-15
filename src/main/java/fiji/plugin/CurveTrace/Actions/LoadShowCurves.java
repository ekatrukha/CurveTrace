package fiji.plugin.CurveTrace.Actions;

import fiji.plugin.CurveTrace.CoreClasses.CurveStack;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;


public class LoadShowCurves implements PlugIn
{
	

	/** main image window */
	public ImagePlus imp; 	
	/** currently active processor */
	public ImageProcessor ip; 
	
	@Override
	public void run(String arg0) 
	{
				
		/** main object storing all curves **/
		CurveStack curvestack = new CurveStack();
		
		if(!curvestack.loadCurvesFromRT())
			return;
		
		if(!curvestack.showCurvesOverlay(1,false))
			return;
		
		
	}

}
