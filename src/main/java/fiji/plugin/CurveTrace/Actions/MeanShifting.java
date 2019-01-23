package fiji.plugin.CurveTrace.Actions;

import java.util.ArrayList;

import fiji.plugin.CurveTrace.CoreClasses.Curve;
import fiji.plugin.CurveTrace.CoreClasses.CurveStack;
import fiji.plugin.CurveTrace.CoreClasses.Normale;
import fiji.plugin.CurveTrace.CoreClasses.Point;
import fiji.plugin.CurveTrace.KDTree.kdtree.TwoDTree;
import fiji.plugin.CurveTrace.KDTree.twod.TwoDRectangle;
import fiji.plugin.CurveTrace.KDTree.IPoint;
import fiji.plugin.CurveTrace.KDTree.kdtree.TwoDFactory;
import ij.IJ;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

public class MeanShifting implements PlugIn {
	
	/**2D tree to store points positions
	 * for acceleration of nearby area search **/
	public TwoDTree tree;
	/**averaging kernel radius in pixels**/
	public double dRadius;
	/**Angle of normales "coalignment" between two points during averaging (in radians).
	 * If the difference in normale angles is higher that that,
	 *  the point will not be included in the averaging (in radians) **/
	public double dAngle;	
	/** whether to render all curves on one image
	 * **/	
	public boolean bMSIgnoreFrame;

	/** whether to update (rewrite) Results table with new curves
	 * after Mean shifting**/
	public boolean bMSUpdateResults;
	
	@Override
	public void run(String arg) 
	{
		int i,j,k;
		int nCount;
		Curve curve;
		Point point;
		float x,y;
		float [] new_coords;
		ArrayList<IPoint> locatedP;
		TwoDRectangle area;
				
		/** main object storing all curves **/
		CurveStack curvestack = new CurveStack();
		
		if(!curvestack.loadCurvesFromRT())
			return;
		
		// check if there is an image	
		if(null == IJ.getImage())
		{
			IJ.noImage();	    
			return;
		}
		
		if(!showMeanShiftingDialog())
			return;
		
		

		//generate 2d tree from all the points
		IJ.showStatus("Generating 2D tree of point...");
		tree=TwoDFactory.generate(curvestack.getAllPointsInArray());
		
		IJ.showStatus("Performing mean shifting...");
		IJ.showProgress(0, (int)curvestack.nPoints);
		
		nCount=0;
		
		//looping through frames		
		for (i=0;i<curvestack.size();i++)
		{
			//looping through curves
			for (j=0;j<curvestack.get(i).size();j++)
			{
				curve=curvestack.get(i).get(j);
				//looping through points
				for (k=0;k<curve.size();k++)				 
				{
					point=curve.get(k);
					x=point.coords[0];
					y=point.coords[1];
					area=new TwoDRectangle(x-3.*dRadius,y-3.*dRadius,x+3.*dRadius,y+3.*dRadius); 
					locatedP = tree.search(area);
					if(locatedP.size()>1)
						{ averagePosition(locatedP,point);}
					

					
					nCount++;
					IJ.showProgress(nCount,(int)curvestack.nPoints);
				}
			}
		}
		IJ.showStatus("Performing mean shifting...Done.");
		IJ.showProgress((int)curvestack.nPoints, (int)curvestack.nPoints);
		if(bMSUpdateResults)
			curvestack.toResultsTable();
		
		if(!curvestack.showCurvesOverlay(1,bMSIgnoreFrame))
			return;
		
		
	}
	public boolean showMeanShiftingDialog()
	{
		GenericDialog shiftingD = new GenericDialog("Mean-shifting parameters");
		shiftingD.addNumericField("Kernel's radius:", Prefs.get("CurveTrace.meanShiftR", 5), 0, 3," pixels");
		shiftingD.addNumericField("Angle cone (0-90):", Prefs.get("CurveTrace.meanShiftAngle", 45), 0, 3," degrees");
		shiftingD.addCheckbox("Ignore frame number?", Prefs.get("CurveTrace.bMSIgnoreFrame", false));
		shiftingD.addCheckbox("Update Results Table", Prefs.get("CurveTrace.bMSUpdateResults", false));
		shiftingD.setResizable(false);
		shiftingD.showDialog();						
		if (shiftingD.wasCanceled())
            return false;
		
		dRadius =  shiftingD.getNextNumber();
		Prefs.set("CurveTrace.meanShiftR", dRadius);
		dAngle = shiftingD.getNextNumber();
		Prefs.set("CurveTrace.meanShiftAngle", dAngle);
		dAngle*=Math.PI/180;
		bMSIgnoreFrame = shiftingD.getNextBoolean();
		Prefs.set("CurveTrace.bMSIgnoreFrame", bMSIgnoreFrame);
		bMSUpdateResults = shiftingD.getNextBoolean();
		Prefs.set("CurveTrace.bMSUpdateResults", bMSUpdateResults);
		return true;
	}
	/** function modifies coordinates of point_in by calculating mean shift
	 * from locatedP array (weighted average with kernel of specified radius and angle cone**/
	public void averagePosition(ArrayList<IPoint> locatedP, Point point_in)
	{
		float xaver=0;
		float yaver=0;
		float weightaver=0;
		double anglex,lenx,weightx;
		double x,y;
		Normale angle;
		int i;
		Point locPoint;
		
		angle=point_in.angle;
		x=point_in.coords[0];
		y=point_in.coords[1];
		for (i=0;i<locatedP.size();i++)
		{
			/*
			 * TODO correct angle calculations,
			 * insert Normale class version
			*/
			locPoint =(Point)locatedP.get(i);
			/*
			anglex=Math.cos(angle-locPoint.angle.val);
			anglex=Math.acos(anglex);
			anglex=Math.min(anglex, Math.PI-anglex);
			*/
			anglex=Normale.Sdist(angle, locPoint.angle);
			anglex=Math.exp(-0.5*Math.pow(anglex/dAngle,2));
			lenx=Math.sqrt(Math.pow(x-locPoint.coords[0],2)+Math.pow(y-locPoint.coords[1],2));
			lenx = Math.exp(-0.5*Math.pow(lenx/dRadius,2));
			weightx=anglex*lenx;
			weightaver+=weightx;
			xaver+=locPoint.coords[0]*weightx;
			yaver+=locPoint.coords[1]*weightx;
		}
		point_in.coords[0]=xaver/weightaver;
		point_in.coords[1]=yaver/weightaver;
		
	}
}
