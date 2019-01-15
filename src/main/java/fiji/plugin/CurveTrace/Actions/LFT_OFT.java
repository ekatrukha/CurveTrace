package fiji.plugin.CurveTrace.Actions;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
/** plugin returns LFT and OFT transforms of an image **/

public class LFT_OFT implements PlugIn {
	/** main image window */
	public ImagePlus imp; 	
	/** currently active processor */
	public ImageProcessor ip;
	
	/** radius of filter **/
	double dRadius;
	/** number of filter orientations **/
	double dOrientN;
	
	@Override
	public void run(String arg) {
		
		int nW,nH;
		double dAngleInterval;
		double dMaxAngleIntensity, dLineSum, dMaxIntensityAngle;
			
		int i,j,q;
		
		int x=0,y=0;
		
		double k;
		double halfPi;
		int R;
		float dRho, dTheta;
		ImagePlus impLFT, impLFT_Orient, impOFT;
		
		FloatProcessor ipLFT, ipLFT_Orient, ipOFT;
		
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
		
		//show parameters dialog
		if(!showLFTOFTDialog())
			return;
		
		
		
		nW=imp.getWidth();
		nH=imp.getHeight();
		
		// get active imageprocessor
		ip = imp.getProcessor();
		//for now
		//convert it to float
		ip =(FloatProcessor) ip.duplicate().convertToFloat();
		ipLFT = new FloatProcessor(nW,nH);
		ipLFT_Orient=new FloatProcessor(nW,nH);
		ipOFT =new FloatProcessor(nW,nH);
		
		dAngleInterval = Math.PI/dOrientN;
		//just int value of dRadius;
		R=(int)dRadius;
		halfPi=Math.PI*0.5;

		/*let's make a filter and put it into array,
		 * so we don't have to calculate it each time */
		int [][][] filter;
		
		int ki, qi;
		filter = new int[(int)dOrientN][2*R+1][2];
		
		ki=0;
		
		for(k=-halfPi; k<halfPi; k=k+dAngleInterval)
        {
			qi=0;
        	for(q=-R; q<R+1; q++)
            {
        		filter[ki][qi][0]=(int)Math.floor(q *Math.cos(k)+0.5);
        		filter[ki][qi][1]=(int)Math.floor(q *Math.sin(k)+0.5);

                qi++;
            }
            ki++;
        }
		
		IJ.showStatus("Calculating LFT...");
		/* LFT */
		for(i=0; i<nW; i++)
		    {
			
			IJ.showProgress(i, nW);
			for(j=0; j<nH; j++)
		        { 
		        	dMaxAngleIntensity = 0;
		            dMaxIntensityAngle = 0;
		            ki=0;
		    		
		            for(k=-halfPi; k<halfPi; k=k+dAngleInterval)
		            {
		            	dLineSum = 0;
		            	qi=0;
		            	
		            	for(q=-R; q<R+1; q++)
	                    {
	                        x = i + filter[ki][qi][0];
	                        y = j - filter[ki][qi][1];
	                        
	                        if(!(x<0 || y<0 || x>(nW-1)||y>(nH-1)))	                    	                        	
	                        	dLineSum+=ip.getf(x,y);

	                        
	                        qi++;
	                    }
		            	ki++;
		   
	                    if(dMaxAngleIntensity<dLineSum) 
	                    {
	                        dMaxAngleIntensity = dLineSum;
	                        dMaxIntensityAngle = k;
	                    }
		            }
		            ipLFT.setf(i, j, (float)(dMaxAngleIntensity/(double)(R+R+1)));
		            ipLFT_Orient.setf(i, j, (float)dMaxIntensityAngle);

		        }//j cycle
		    }//i cycle
		
		
		
		 /* OFT */
		IJ.showStatus("Calculating OFT...");
		for(i=0; i<nW; i++)
	    {
			IJ.showProgress(i, nW);
			for(j=0; j<nH; j++)
	        {
					dRho = ipLFT.getf(i,j);
	                dMaxAngleIntensity = 0;
	                ki=0;
	                for(k=-halfPi; k<halfPi; k=k+dAngleInterval)
	                {
	                    dLineSum = 0;
	                    qi=0;
	                    for(q=-R; q<R+1; q++)
	                    {
	                        
	                        x = i + filter[ki][qi][0];
	                        y = j - filter[ki][qi][1];
	                        if(!(x<0 || y<0 || x>(nW-1)||y>(nH-1))) 
	                        {
	                        	dTheta = ipLFT_Orient.getf(x, y);   	
	                        	dLineSum = dLineSum + dRho * Math.cos(2.0*(dTheta-k));
	                     
	                        }
	                        qi++;
	                    }
	                    ki++;
	                    if(dMaxAngleIntensity<dLineSum) 
	                    	dMaxAngleIntensity = dLineSum;
	                }
	                ipOFT.setf(i,j,(float)dMaxAngleIntensity);
	            
	        }
	    }
		String params= "_R="+Integer.toString((int)dRadius)+"_Norient="+Integer.toString((int)dOrientN);
	    impLFT = new ImagePlus(imp.getTitle()+"_LFT"+params,ipLFT);
	    impLFT_Orient = new ImagePlus(imp.getTitle()+"_LFT_orient"+params, ipLFT_Orient);
	    impOFT = new ImagePlus(imp.getTitle()+"_OFT"+params, ipOFT);
		IJ.showStatus("Calculating LFT/OFT...Done.");
		IJ.showProgress(1, 1);
	    impLFT.show();
	    impLFT_Orient.show();
	    impOFT.show();
	}
	
	public boolean showLFTOFTDialog()
	{
		
		GenericDialog lftoftD = new GenericDialog("LFT/OFT parameters");
		
		lftoftD.addNumericField("Filter radius:", Prefs.get("CurveTrace.lftoftRadius", 10), 0, 3,"pixels");
		lftoftD.addNumericField("Number of filter orientations (10-30):", Prefs.get("CurveTrace.lftoftNOrient", 20), 0, 3,"");
		lftoftD.setResizable(false);
		lftoftD.showDialog();	
		if (lftoftD.wasCanceled())
            return false;
		
		dRadius = lftoftD.getNextNumber();
		Prefs.set("CurveTrace.lftoftRadius", dRadius);
		dOrientN = Math.round(lftoftD.getNextNumber());
		Prefs.set("CurveTrace.lftoftNOrient", dOrientN);
	
		
		return true;
	}
	

}
