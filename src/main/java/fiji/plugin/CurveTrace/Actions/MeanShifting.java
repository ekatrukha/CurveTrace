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
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class MeanShifting implements PlugIn {
	
	/**2D tree to store points positions
	 * for acceleration of nearby area search **/
	public TwoDTree tree;
	
	/** what to use for mean shifting 0 = curves 1 = image **/	
	int nMSType;
	
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
	
	/** main image window */
	ImagePlus imp; 
	
	ImageProcessor ip;
	int nW, nH;
	
	Normale angles[][];
	
	double lenxarr[][];
	
	float imvals[][];
	
	@Override
	public void run(String arg) 
	{
		int i,j,k;
		int nCount;
		Curve curve;
		Point point;
		float x,y;
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
		
		imp = IJ.getImage();
		ip=imp.getProcessor().convertToFloat();
		nW=ip.getWidth();
		nH=ip.getHeight();
		
		if(!showMeanShiftingDialog())
			return;
		
		

		//generate 2d tree from all the points
		IJ.showStatus("Generating 2D tree of point...");
		tree=TwoDFactory.generate(curvestack.getAllPointsInArray());
		
		IJ.showStatus("Performing mean shifting...");
		IJ.showProgress(0, (int)curvestack.nPoints);
		
		nCount=0;
		
		if(nMSType==1)
			fillArrays();
		
		//looping through frames		
		for (i=0;i<curvestack.size();i++)
		{
			if(!bMSIgnoreFrame && nMSType==1)
			{
				imp.setSliceWithoutUpdate(curvestack.get(i).nFrame);
				ip=imp.getProcessor().convertToFloat();
			}
			//looping through curves
			for (j=0;j<curvestack.get(i).size();j++)
			{
				curve=curvestack.get(i).get(j);
				//looping through points
				for (k=0;k<curve.size();k++)				 
				{
					point=curve.get(k);
					
					if(nMSType==0)
					{
						x=point.coords[0];
						y=point.coords[1];
						area=new TwoDRectangle(x-3.*dRadius,y-3.*dRadius,x+3.*dRadius,y+3.*dRadius); 
						locatedP = tree.search(area);
						if(locatedP.size()>1)
						{ 
							averagePosition(locatedP,point);
						}
					}
					else
					{
						averagePositionImage(point);
					}

					
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

	/** function modifies coordinates of point_in by calculating mean shift
	 * from locatedP array (weighted average with kernel of specified radius and angle cone**/
	public void averagePosition(ArrayList<IPoint> locatedP, Point point_in)
	{
		float xaver=0;
		float yaver=0;
		float weightaver=0;
		double anglex,lenx, weightx;
		double x,y;
		//double dx,dy, directionx;
		Normale angle;//, direction;
		int i;
		Point locPoint;
		
		angle=point_in.angle;
		x=point_in.coords[0];
		y=point_in.coords[1];
		for (i=0;i<locatedP.size();i++)
		{

			locPoint =(Point)locatedP.get(i);
			//see if it is the same point
//			dx=locPoint.coords[0]-x;
//			dy=locPoint.coords[1]-y;
//			
//			if(Math.abs(dx)<0.0001 && Math.abs(dy)<0.0001)
//			{
//				weightx=1.0;
//				
//			}
//			else
//			{
				//contribution of angle difference between curves normales
				anglex=Normale.Sdist(angle, locPoint.angle);
				anglex=Math.exp(-0.5*Math.pow(anglex/dAngle,2));
				
				//contribution of direction	
//				direction = new Normale(Math.atan2(dy, dx));
//				directionx = Normale.Sdist(angle, direction);
//				directionx = Math.exp(-0.5*Math.pow(directionx/dAngle,2));
				
				//contribution of distance between points
				lenx=Math.sqrt(Math.pow(x-locPoint.coords[0],2)+Math.pow(y-locPoint.coords[1],2));
				lenx = Math.exp(-0.5*Math.pow(lenx/dRadius,2));
				
				weightx=anglex*lenx;
//				weightx=anglex*lenx*directionx;
////			}
			
				//angle=Normale.averageWeighted(angle,locPoint.angle,anglex);
			
			weightaver+=weightx;
			xaver+=locPoint.coords[0]*weightx;
			yaver+=locPoint.coords[1]*weightx;
		}
		point_in.coords[0]=xaver/weightaver;
		point_in.coords[1]=yaver/weightaver;
		//point_in.angle=angle;
		
	}
	
	public void fillArrays()
	{
		double dx,dy;
		double lenx;
		lenxarr = new double[(int)(6.0*dRadius+1.0)][(int)(6.0*dRadius+1.0)];
		angles = new Normale[(int)(6.0*dRadius+1.0)][(int)(6.0*dRadius+1.0)];
		imvals = new float[(int)(6.0*dRadius+1.0)][(int)(6.0*dRadius+1.0)];
		int xi, yi;
		xi=0;
		
		for (dx=-3.0*dRadius;dx<3.0*dRadius+1.0;dx++)
		{
			yi=0;
			for (dy=-3.0*dRadius;dy<3.0*dRadius+1.0;dy++)
			{
				
				lenx=Math.sqrt(Math.pow(dx,2)+Math.pow(dy,2));
				lenxarr[xi][yi]=Math.exp(-0.5*Math.pow(lenx/dRadius,2));
				angles[xi][yi]=new Normale(Math.atan2(dy, dx));
				yi++;
			}
			xi++;
		}
	}
	
	public void averagePositionImage(Point point_in)
	{
		double x,y;
		double dx,dy;
		Normale angle,direction;
		double anglex, lenx,weightx;
		float weightaver=0;
		float xaver=0;
		float yaver=0;
		float intMin = Float.MAX_VALUE;
		
		int xi, yi;
		xi=0;
		
		
		x=point_in.coords[0];
		y=point_in.coords[1];
		angle=point_in.angle;
		
		
		for (dx=x-3.0*dRadius;dx<x+3.0*dRadius;dx++)
		{
			yi=0;
			for (dy=y-3.0*dRadius;dy<y+3.0*dRadius;dy++)
			{
				if(!(dx<0 || dy<0 || dx>(nW-1)||dy>(nH-1)))	
				{
					imvals[xi][yi]=ip.getf((int)dx,(int)dy);
					if(imvals[xi][yi]<intMin)
					{
						intMin=imvals[xi][yi];
					}
				}
				yi++;
			}
			xi++;
		}
		
		xi=0;
		for (dx=x-3.0*dRadius;dx<x+3.0*dRadius;dx++)
		{
			yi=0;
			for (dy=y-3.0*dRadius;dy<y+3.0*dRadius;dy++)
			{
				if(!(dx<0 || dy<0 || dx>(nW-1)||dy>(nH-1)))	
				{
					
					if(Math.abs(dx-x)<0.0001 && Math.abs(dy-x)<0.0001)
					{
						anglex=1.0;
						
					}
					else
					{
						//dLineSum+=ip.getf(dx,dy);
						//direction = new Normale(Math.atan2(dy-y, dx-x));
						direction =angles[xi][yi];
						anglex=Normale.Sdist(angle, direction);
						anglex=Math.exp(-0.5*Math.pow(anglex/dAngle,2));
					}
					//lenx=Math.sqrt(Math.pow(x-dx,2)+Math.pow(y-dy,2));
					//lenx = Math.exp(-0.5*Math.pow(lenx/dRadius,2));
					lenx = lenxarr[xi][yi];
					//weightx=anglex*lenx*((double)ip.getf((int)dx,(int)dy));
					weightx=anglex*lenx*(imvals[xi][yi]-intMin);
					weightaver+=weightx;
					xaver+=dx*weightx;
					yaver+=dy*weightx;
				}
				yi++;
			}
			xi++;
		}
		point_in.coords[0]=xaver/weightaver;
		point_in.coords[1]=yaver/weightaver;
	}
	
	
	
	public boolean showMeanShiftingDialog()
	{
		GenericDialog shiftingD = new GenericDialog("Mean-shifting parameters");
		
		String [] sMSType = new String [] {"Curves","Image"};
		shiftingD.addChoice("Foe MS calculations use: ",sMSType, Prefs.get("CurveTrace.nMSType", "Curves"));
		shiftingD.addNumericField("Kernel's radius:", Prefs.get("CurveTrace.meanShiftR", 5), 0, 3," pixels");
		shiftingD.addNumericField("Angle cone (0-90):", Prefs.get("CurveTrace.meanShiftAngle", 25), 0, 3," degrees");
		shiftingD.addCheckbox("Ignore frame number?", Prefs.get("CurveTrace.bMSIgnoreFrame", false));
		shiftingD.addCheckbox("Update Results Table", Prefs.get("CurveTrace.bMSUpdateResults", false));
		shiftingD.setResizable(false);
		shiftingD.showDialog();						
		if (shiftingD.wasCanceled())
            return false;
		
		nMSType = shiftingD.getNextChoiceIndex();
		Prefs.set("CurveTrace.nMSType", sMSType[nMSType]);
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
}
