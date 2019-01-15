package fiji.plugin.CurveTrace.Actions;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class OFT_vote implements PlugIn {
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
		double  dLineSum;
		double dMeanResponse;
		int i,j,q;
		
		int x=0,y=0;
		
		double k;
		double halfPi;
		int R;
		float dRho, dTheta;
		ImagePlus impOFT;
		float [][] iLFTall;
		FloatProcessor ipOFT;
		
		
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
		if(!showOFTvDialog())
			return;
		
		
		
		nW=imp.getWidth();
		nH=imp.getHeight();
		
		// get active imageprocessor
		ip = imp.getProcessor();
		//for now
		//convert it to float
		ip =(FloatProcessor) ip.duplicate().convertToFloat();
	
		ipOFT =new FloatProcessor(nW,nH);
		iLFTall = new float[nW*nH][(int)dOrientN];
		
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
		            ki=0;	
		            dMeanResponse=0;
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
		            	
		            	dLineSum = dLineSum/(2.0*(double)R+1.0);
		            	//y*width+x
		            	iLFTall[j*nW+i][ki]=(float) dLineSum;
		            	dMeanResponse+=dLineSum;
		            	
		            	ki++;
		            }// end for(k=-halfPi; k<halfPi; k=k+dAngleInterval)
		            
		            //subtracting average
		            dMeanResponse = dMeanResponse/(double)ki;
		            ki=0;			            
		            for(k=-halfPi; k<halfPi; k=k+dAngleInterval)
		            {
		            	iLFTall[j*nW+i][ki]-=dMeanResponse;
		            }

		        }//j cycle
		    }//i cycle
		
		//calculated OFT at each point, now to 
		// OFT voting 
		IJ.showStatus("Calculating OFT...");
		
		double k2;
		int ki2;
		//float dPixVal;
		double dPixTotVal;
		double dTetVal;
		double dAnglDiff;
		double dN;
		
		for(i=0; i<nW; i++)
	    {
			IJ.showProgress(i, nW);
			for(j=0; j<nH; j++)
	        {
					
	                ki=0;
	                dPixTotVal=0;
	                for(k=-halfPi; k<halfPi; k=k+dAngleInterval)
	                {
	                    qi=0;
	                    for(q=-R; q<R+1; q++)
	                    {
	                    	x = i + filter[ki][qi][0];
	                        y = j - filter[ki][qi][1];
	                        
	                        if(!(x<0 || y<0 || x>(nW-1)||y>(nH-1)))	  
	                        {
	                        	//ok, let's sum up stuff
	                        	ki2=0;
	                        	for(k2=-halfPi; k2<halfPi; k2=k2+dAngleInterval)
	        	                {
	                        		dTetVal = iLFTall[j*nW+i][ki2];
	                        		//only positive values
	                        		if(dTetVal>0 && q!=0)
	                        		{
	                        			dAnglDiff=Math.cos(k-k2);
	                					dAnglDiff=Math.acos(dAnglDiff);
	                					dAnglDiff=Math.min(dAnglDiff, Math.PI-dAnglDiff);
	                					dN = Math.pow(Math.abs(q), -0.1)*Math.exp(-Math.pow(dAnglDiff/(1.41*0.26),2));
	                					dPixTotVal+=dTetVal*dN*iLFTall[y*nW+x][ki2];
	                        		}
	                        		
	                        		
	                        		ki2++;	        	                
	        	                }
	                        	
	                        }
	                        
	                       
	                        qi++;
	                    }
	                    ki++;
	                    
	                }
	                ipOFT.setf(i,j,(float)dPixTotVal);
	            
	        }
	    }
		
		String params= "_R="+Integer.toString((int)dRadius)+"_Norient="+Integer.toString((int)dOrientN);
	    //impLFT = new ImagePlus(imp.getTitle()+"_LFT"+params,ipLFT);
	    //impLFT_Orient = new ImagePlus(imp.getTitle()+"_LFT_orient"+params, ipLFT_Orient);
	    impOFT = new ImagePlus(imp.getTitle()+"_OFT"+params, ipOFT);
		IJ.showStatus("Calculating OFT vote...Done.");
		IJ.showProgress(1, 1);
	    //impLFT.show();
	    //impLFT_Orient.show();
	    impOFT.show();
	}
	
	public boolean showOFTvDialog()
	{
		
		GenericDialog otfvD = new GenericDialog("OFT voting filter parameters");
		
		otfvD.addNumericField("Filter radius:", Prefs.get("CurveTrace.otfvRadius", 10), 0, 3,"pixels");
		otfvD.addNumericField("Number of filter orientations (10-30):", Prefs.get("CurveTrace.otfvNOrient", 20), 0, 3,"");
		otfvD.setResizable(false);
		otfvD.showDialog();	
		if (otfvD.wasCanceled())
            return false;
		
		dRadius = otfvD.getNextNumber();
		Prefs.set("CurveTrace.otfvRadius", dRadius);
		dOrientN = Math.round(otfvD.getNextNumber());
		Prefs.set("CurveTrace.otfvNOrient", dOrientN);
	
		
		return true;
	}
	

}
