package fiji.plugin.CurveTrace.Actions;

import fiji.plugin.CurveTrace.Fit.OneDGaussianBG;
import fiji.plugin.CurveTrace.CoreClasses.Curve;
import fiji.plugin.CurveTrace.CoreClasses.CurveSet;
import fiji.plugin.CurveTrace.CoreClasses.CurveStack;
import fiji.plugin.CurveTrace.CoreClasses.Point;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import jaolho.data.lma.LMA;

public class FitGaussian implements PlugIn {
	/** main image window */
	public ImagePlus imp; 	
	/** currently active processor */
	public ImageProcessor ip;
	
	/** main object storing all curves **/
	CurveStack curvestack;
	
	/** initial value of SD of line profile**/
	double dSD;
	/** Length of tangent fitting line, pixels**/
	int nLineLength;
	/** Sampling step **/
	double dSamplingStep;
	/** Line width, pixels**/
	int nLineWidth;
	/** intensity interpolation method: 0 "None", 1 "Bilinear", 2 "Bicubic" **/
	int nIterpolationMethod;
	/** type of curves added to overlay 
	 * 0 "Nothing", 1 "Only lines", 2 "Lines of average width", 3 "Lines with width border"**/
	int nOverlayLinesType;
	
	boolean bIgnoreFrame;
	
	/** if current image is float (32-bit) **/	
	boolean bImageFloat=false;
	/** Maximum displacement after fit. Fits displaced further than that will be discarded **/
	float dDispThresh;
	
	
	@Override
	public void run(String arg) 
	{
		
		int nFrame;
		int nCurve;
		int nPoint;
		long nPointCount;
		CurveSet curveset;
		Curve curve;
		Point point;
		int i,j;
		float x,y,nx,ny,tx,ty,angle,xp,yp;
		double intP,maxI, minI, meanI, sstot;
		double dFitMean=0., dFitSD=0., dFitAmp=0.,dFitBG=0., dFitR2=0.;
		boolean bFitException;
		
		LMA fitlma;
		
		curvestack = new CurveStack();

		//load curves from results table
		curvestack.loadCurvesFromRT();

		// get active image
		imp = IJ.getImage();
		if(null == imp)
		{
			IJ.noImage();	    
			return;
		}
		else if (imp.getType() != ImagePlus.GRAY8 && imp.getType() != ImagePlus.GRAY16 && imp.getType() != ImagePlus.GRAY32) 
		{
		    IJ.error("8, 16 or 32-bit grayscale image required");
		    return;
		}	
		
		if(imp.getType() == ImagePlus.GRAY32)
		{
			bImageFloat=true;
		}
		
		//ask for fitting parameters
		if(!showFittingDialog())
			return;
		//check total frame number
		if(curvestack.nFramesTotal>imp.getStackSize() && !bIgnoreFrame)
		{
			IJ.error("Number of frames/slices in the image is less than the values in the Results table");
			return;
		}
		
		IJ.showStatus("Fitting lines...");
		IJ.showProgress(0, (int)curvestack.nPoints);

		//setting up averaging/profile arrays
		double [] averagew = new double[nLineWidth];
		double nHalfWidth = ((float)nLineWidth-1.0)*0.5;
		for(i=0;i<averagew.length;i++)
			averagew[i]=(float)i-nHalfWidth;
		
		double [] linew = new double[nLineLength];
		double [][] dFitData = new double [2][nLineLength];
		double nHalfLength = ((float)nLineLength-1.0)*0.5;
		for(i=0;i<linew.length;i++)
		{
			linew[i]=(float)i-nHalfLength;
			dFitData[0][i]=(double)i-nHalfLength;
		}
		//
		
		nPointCount=0;
		//frame cycle
		for (nFrame=1;nFrame<=curvestack.size();nFrame++)
		{
			curveset=curvestack.get(nFrame-1);
			if(!bIgnoreFrame)
			{
				imp.setSliceWithoutUpdate(curveset.nFrame);
			}
			ip = imp.getProcessor().duplicate();
			ip.setInterpolationMethod(nIterpolationMethod);
			//curve cycle
			for (nCurve=0; nCurve<curveset.size(); nCurve++)
			{
				curve = curveset.get(nCurve);
				
				//point cycle
				for (nPoint =0;nPoint<curve.size();nPoint++)
				{
					point=curve.get(nPoint);
					x=point.coords[0];
					y=point.coords[1];
					angle=(float) point.angle.val;
					
					//normal vector
					nx = (float)Math.sin(angle);
					ny = (float)Math.cos(angle);
					
					//tangent vector (turned 90 degrees)
					tx = -ny;
					ty = nx;
					maxI=Double.MIN_VALUE;
					minI=Double.MAX_VALUE;
					meanI=0;
					for(i=0;i<nLineLength;i++)
					{
						//calculating average intensity across the line
						intP=0;
						xp=(float) (x+nx*linew[i]);
						yp=(float) (y+ny*linew[i]);
						for (j=0;j<nLineWidth;j++)
						{
							
							if (bImageFloat)
								intP+=Float.intBitsToFloat(ip.getPixelInterpolated(xp+tx*averagew[j], yp+ty*averagew[j]));
							else							
								intP+=ip.getPixelInterpolated(xp+tx*averagew[j], yp+ty*averagew[j]);
						}
						intP=intP/(float)nLineWidth;
						if(intP>maxI)
							maxI=intP;
						if(intP<minI)
							minI=intP;
						dFitData[1][i]=intP;
						meanI+=intP;
					}
					fitlma = new LMA(new OneDGaussianBG(), new double[] {maxI, 0, dSD,minI}, dFitData);
					fitlma.maxIterations=10;
					bFitException=false;
					try
					{
						fitlma.fit();
					}
					catch(Exception e)
					{
						bFitException=true;
					}
					if(!bFitException)
					{
						dFitAmp  = fitlma.parameters[0];
						dFitMean = fitlma.parameters[1];
						dFitSD   = fitlma.parameters[2];
						dFitBG   = Math.abs(fitlma.parameters[3]);
					}
					else 
						dFitMean=dDispThresh*2.0;
					
					
					if(Math.abs(dFitMean)<dDispThresh && dFitAmp>0)
					{
						meanI=meanI/(double)nLineLength;
						//calculate R2
						sstot=0;
						for(i=0;i<nLineLength;i++)
						{
							sstot+=Math.pow(dFitData[1][i]-meanI, 2);
						}
						dFitR2=1.0 - (fitlma.chi2/sstot);
						point.coords[0]=(float) (point.coords[0]+nx*dFitMean);
						point.coords[1]=(float) (point.coords[1]+ny*dFitMean);
						point.bg=(float) dFitBG;
						point.amp=(float) dFitAmp;
						point.width=(float) dFitSD;
						point.R2=(float) dFitR2;
						point.integr_int =(float) ((dFitAmp*dFitSD)*Math.sqrt(2.0*Math.PI));
					}
					else
					{
						//point.amp=(float) maxI;						
						//point.bg=(float) minI;
						//point.width = (float) 0.0;
						//point.R2=(float) 0.0;
						point.amp=(float) (-1.0*(maxI-minI));						
						point.bg=(float) (-1.0*minI);
						point.width = (float) -1.0;
						point.R2=(float) -1.0;

					}
					nPointCount++;
					IJ.showProgress((int)nPointCount, (int)curvestack.nPoints);
				}//point cycle end
				
			}//curve cycle end
			
		}//frame cycle end
		
		IJ.showStatus("Fitting lines...Done.");
		IJ.showProgress((int)curvestack.nPoints, (int)curvestack.nPoints);

		curvestack.showCurvesOverlay(nOverlayLinesType,bIgnoreFrame);
		curvestack.toResultsTable();
	} 
	public boolean showFittingDialog()
	{
		
		GenericDialog fittingD = new GenericDialog("Fitting parameters");
		String [] InterpolationMethods = new String [] {
				"None", "Bilinear","Bicubic"};
		String [] OverlayLineType = new String [] {
				"Nothing", "Only lines", "Lines of average width", "Lines with width border"};
		
		fittingD.addNumericField("Approximate SD of line profile:", Prefs.get("CurveTrace.fitGsigma", 1.6), 2, 4," pixels");
		fittingD.addNumericField("Length of tangent fitting line:", Prefs.get("CurveTrace.fitGlinelength", 8), 0, 3," pixels");
		fittingD.addNumericField("Sampling step:", Prefs.get("CurveTrace.fitGsamplestep", 0.5), 2, 4," pixels");
		fittingD.addNumericField("Line width:", Prefs.get("CurveTrace.fitGlinewidth", 3), 0, 3," pixels");
		fittingD.addChoice("Intensity interpolation method:", InterpolationMethods, Prefs.get("CurveTrace.fitGinterpolation", "Bilinear"));
		fittingD.addNumericField("Maximum displacement threshold (fitting is bad):", Prefs.get("CurveTrace.fitDispThresh", 1.6), 2, 4," pixels");
		
		fittingD.addChoice("After fitting add to overlay:", OverlayLineType, Prefs.get("CurveTrace.OverlayLineType", "Only lines"));
		fittingD.addCheckbox("Ignore frame number", Prefs.get("CurveTrace.bIgnoreFrame", false));
		fittingD.setResizable(false);
		fittingD.showDialog();						
		if (fittingD.wasCanceled())
            return false;
		
		dSD = fittingD.getNextNumber();
		Prefs.set("CurveTrace.fitGsigma", dSD);
		nLineLength = (int)fittingD.getNextNumber();
		Prefs.set("CurveTrace.fitGlinelength", nLineLength);
		dSamplingStep = fittingD.getNextNumber();
		Prefs.set("CurveTrace.fitGsamplestep", dSamplingStep);
		nLineWidth = (int)fittingD.getNextNumber();
		Prefs.set("CurveTrace.fitGlinewidth", nLineWidth);
		
		nIterpolationMethod = fittingD.getNextChoiceIndex();
		Prefs.set("CurveTrace.fitGinterpolation", InterpolationMethods[nIterpolationMethod]);
		
		dDispThresh = (float) fittingD.getNextNumber();
		Prefs.set("CurveTrace.fitDispThresh", dDispThresh);
		
		nOverlayLinesType = fittingD.getNextChoiceIndex();
		Prefs.set("CurveTrace.OverlayLineType", OverlayLineType[nOverlayLinesType]);
		bIgnoreFrame = fittingD.getNextBoolean();
		Prefs.set("CurveTrace.bIgnoreFrame", bIgnoreFrame);

		return true;
	}
}
