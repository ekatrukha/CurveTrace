package fiji.plugin.CurveTrace.CoreClasses;

import java.awt.Color;
import java.awt.Font;
import java.util.ArrayList;


import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.measure.ResultsTable;
import ij.plugin.filter.Analyzer;

/**
 *	CurveStack class, contains array of CurveSets
 *  (corresponding to different frames)
 */
public class CurveStack extends ArrayList<CurveSet> 
{
	/** Results table **/
	ResultsTable ptable;
	/** total number of frames **/
	public int nFramesTotal;
	/** total number of points **/
	public long nPoints;
	
	public boolean bPointsCounted = false; 
	
	/** total number of curves **/
	public long nCurves;

	/** lock on results table **/
	java.util.concurrent.locks.Lock ptable_lock = new java.util.concurrent.locks.ReentrantLock();
	public CurveStack()
	{
		super();
	}
	
	/** function loads curve coordinates and angles from Results table
	 * to the current object **/
	public boolean loadCurvesFromRT()
	{
		/** x coordinates of curve points **/
		double [] x;
		/** y coordinates of curve points **/
		double [] y;
		/** normal angles **/
		double [] angles;
		/** curves numbers **/
		double [] curveN;	
		/** frames numbers **/
		double [] frameN;
		/** frames numbers **/
		double [] response;

		/** column index**/
		int col_ind;
		/** current frame index **/
		int nFindex=-1;
		/** current curve index **/
		int nCindex=-1;
	
		CurveSet curveset = new CurveSet ();
		Curve curve = new Curve();
		Point point;
		
		ptable = ResultsTable.getResultsTable();
		if (ptable.getCounter()==0 )
		{
			IJ.error("There is no ResultsTable with curves!");
			return false;
		}
		
		// lock results table 
		ptable_lock.lock();
		
		//read results table
		col_ind=ptable.getColumnIndex("X_(px)");
		if(col_ind==-1)
		{
			IJ.error("Table does not contain curves in a proper format (cannot find X_(px) column)");
			return false;
		}
		else
			x   = ptable.getColumnAsDoubles(col_ind);
		
		col_ind=ptable.getColumnIndex("Y_(px)");		
		if(col_ind==-1)
		{
			IJ.error("Table does not contain curves in a proper format (cannot find Y_(px) column)");
			return false;
		}
		else
			y   = ptable.getColumnAsDoubles(col_ind);
		
		col_ind=ptable.getColumnIndex("Angle_of_normal_(radians)");		
		if(col_ind==-1)
		{
			
			IJ.error("Table does not contain curves in a proper format (cannot find Angle_of_normal_(radians)) column");
			return false;
		}
		else
			angles   = ptable.getColumnAsDoubles(col_ind);
		
		col_ind=ptable.getColumnIndex("Contour Number");
		if(col_ind==-1)
		{
			IJ.error("Table does not contain curves in a proper format (cannot find Contour Number column)");
			return false;
		}
		else
			curveN   = ptable.getColumnAsDoubles(col_ind);
		
		col_ind=ptable.getColumnIndex("Frame Number");
		if(col_ind==-1)
		{
			IJ.error("Table does not contain curves in a proper format (cannot find Frame Number column)");
			return false;
		}
		else
			frameN   = ptable.getColumnAsDoubles(col_ind);
		
		col_ind=ptable.getColumnIndex("Response");
		if(col_ind==-1)
		{
			IJ.error("Table does not contain curves in a proper format (cannot find Response column)");
			return false;
		}
		else
			response = ptable.getColumnAsDoubles(col_ind);
		
		nFramesTotal = (int)frameN[frameN.length-1];
		nPoints=x.length;
		bPointsCounted = true;
		nCurves=0;
		for (int i=0;i<nPoints;i++)
		{
			
			//new frame
			if(nFindex!=(int)frameN[i])
			{
				//add previous curveset
				if(i!=0)
				{					
					//add last curve
					curveset.add(curve);
					nCurves++;
					this.add(curveset);
				}
				//reset indexes
				nFindex=(int)frameN[i];
				curveset = new CurveSet();
				curveset.nFrame =(int)frameN[i];
				//IJ.log("frame"+Double.toString(frameN[i]));
				nCindex=-1;				
			}
			if (nCindex!=(int)curveN[i])
			{
				if(nCindex!=-1)
				{
					//add previous curve
					curveset.add(curve);
					nCurves++;
				}
				nCindex=(int)curveN[i];
				curve = new Curve();	
				//IJ.log("Curve"+Double.toString(curveN[i]));
			}
			point = new Point((float)x[i],(float)y[i],(float)angles[i], (float)response[i]);
			curve.add(point);
		}
		//end of the table
		curveset.add(curve);
		nCurves++;
		this.add(curveset);
		
	
		// unlock results table 
		ptable_lock.unlock();
		return true;
	}
	
	/** function adds curve stored inside this object to Overlay 
	 * @param nOverlayLinesType how lines are rendered 
	 * 0 -"Nothing", 1 - "Only lines", 2 - "Lines of average width", 3 - "Lines with width border", 4 - "Lines with normales" 
	 * @param bIgnoreFrame true=render all lines in one frame, false=render at corresponding stack's frame**/
	public boolean showCurvesOverlay(int nOverlayLinesType, boolean bIgnoreFrame)
	{
		/** main image window */
		ImagePlus imp; 	

		/** overlay of main image */
		Overlay image_overlay; 
		
		Curve curve;
		Point point;
		int nCurrFrame;
		Font font;
		
		int i,j,k, pointsN;
		
		PolygonRoi polyline_p;
		TextRoi textRoiN;
		Color[] colors;
		Color[] colors_bord;
		double nx,ny;

		colors = new Color[8];		  
		colors_bord = new Color[8];
		colors[0] = Color.BLUE;
		colors[1] = Color.CYAN;
		colors[2] = Color.GREEN;
		colors[3] = Color.MAGENTA;
		colors[4] = Color.ORANGE;
		colors[5] = Color.PINK;
		colors[6] = Color.RED;
		colors[7] = Color.YELLOW;
		colors_bord[0] = new Color(0,0,139);//dark blue
		colors_bord[1] = new Color(0,139,139); //dark cyan
		colors_bord[2] = new Color(0,100,0); //dark green
		colors_bord[3] = new Color(139,0,139); //dark magenta
		colors_bord[4] = new Color(255,140,0); //dark orange
		colors_bord[5] = new Color(199,21,133); //mediumvioletred
		colors_bord[6] = new Color(139,0,0); //dark red
		colors_bord[7] = new Color(153,153,0); //dark yellow
		
		image_overlay = new Overlay();
		
		font = new Font("SansSerif",Font.PLAIN,6);
		// get active image
		imp = IJ.getImage();
		if(null == imp)
		{
			IJ.noImage();	    
			return false;
		}	
		
		if(this.nFramesTotal>imp.getStackSize()&& !bIgnoreFrame)
		{
			IJ.error("Number of frames/slices in the image is less than the values in the Results table");
			return false;
		}
		//cycling through curvestack
		for (i=0;i<this.size();i++)
		{
			if(!bIgnoreFrame)
				nCurrFrame=this.get(i).nFrame;
			else
				nCurrFrame=imp.getCurrentSlice();
			//cycling through curveset
			for (j=0;j<this.get(i).size();j++)
			{
				 curve=this.get(i).get(j);
				 pointsN=curve.size();
				 float[] px = new float[pointsN];
				 float[] py = new float[pointsN];
				 
				 
				 for (k=0;k<pointsN;k++)
				 {
					 px[k]=curve.get(k).coords[0];
					 py[k]=curve.get(k).coords[1];
				 }
				 textRoiN = new TextRoi(px[0],py[1],Integer.toString(j+1),font);
				 textRoiN.setStrokeColor(colors[j%8]);
				 textRoiN.setPosition(nCurrFrame);
				 image_overlay.add(textRoiN);
				 
				 polyline_p = new PolygonRoi(px, py, Roi.POLYLINE);
				 polyline_p.setPosition(nCurrFrame);
				 polyline_p.setStrokeColor(colors[j%8]);
				 polyline_p.setStrokeWidth(0.0);
				 
				 if(nOverlayLinesType==1 || nOverlayLinesType >2)
					 image_overlay.add(polyline_p);
				 if(nOverlayLinesType ==2)
				 {
					 double line_avg_width = 0;
					 int num_pnt=0;
					 for (k=0;k<pointsN;k++)
					 {	
						 if (curve.get(k).R2>0.0)
						 {
							 line_avg_width +=curve.get(k).width;
							 num_pnt++;
						 }
					 }  
					 line_avg_width =line_avg_width /num_pnt;

					 polyline_p.setStrokeWidth(line_avg_width);
					 polyline_p.updateWideLine((float) line_avg_width);
					 image_overlay.add(polyline_p);
				  }
				 if(nOverlayLinesType==3)
				 {
					  for(k=0;k<pointsN;k++)
					  {
						point=curve.get(k);
						nx = Math.sin(point.angle.val);
				      	ny = Math.cos(point.angle.val);		
					  	px[k]=(float) (point.coords[0] + nx*point.width*2.0);
					  	py[k]=(float) (point.coords[1] + ny*point.width*2.0);
					  }
					  polyline_p = new PolygonRoi(px, py, Roi.POLYLINE);
					  polyline_p.setStrokeColor(colors_bord[j%8]);
					  polyline_p.setPosition(nCurrFrame);
					  polyline_p.setStrokeWidth(0.1);
					  image_overlay.add(polyline_p);				  
	
					  for(k=0;k<pointsN;k++)
					  {
						point=curve.get(k);
						nx = Math.sin(point.angle.val);
				      	ny = Math.cos(point.angle.val);		
					  	px[k]=(float) (point.coords[0] - nx*point.width*2.0);
					  	py[k]=(float) (point.coords[1] - ny*point.width*2.0);
					  }
					  polyline_p = new PolygonRoi(px, py, Roi.POLYLINE);
					  polyline_p.setStrokeColor(colors_bord[j%8]);		
					  polyline_p.setStrokeWidth(0.1);
					  polyline_p.setPosition(nCurrFrame);
					  image_overlay.add(polyline_p);	
					 
				 }
				 if(nOverlayLinesType==4)
				 {
					 px = new float[2];
					 py = new float[2];

					 for(k=0;k<pointsN;k++)
					  {
						point=curve.get(k);
						nx = Math.sin(point.angle.val);
				      	ny = Math.cos(point.angle.val);		
					  	px[0]=(float) (point.coords[0] + nx*3.0);
					  	py[0]=(float) (point.coords[1] + ny*3.0);
					  	px[1]=(float) (point.coords[0] - nx*3.0);
					  	py[1]=(float) (point.coords[1] - ny*3.0);
					  	polyline_p = new PolygonRoi(px, py, Roi.POLYLINE);
						polyline_p.setPosition(nCurrFrame);
						polyline_p.setStrokeColor(colors[j%8]);
						polyline_p.setStrokeWidth(0.0);
						image_overlay.add(polyline_p);
					  } 
				 }
				 /*
				 curve=this.get(i).get(j);
				 pointsN=curve.size();
				 
				  px = new float[3];
				  py = new float[3];
				  for (k=0;k<pointsN;k++)
				  {
					 px[0]=curve.get(k).coords[0]+(float)Math.sin(curve.get(k).angle);
					 py[0]=curve.get(k).coords[1]+(float)Math.cos(curve.get(k).angle);
					 px[1]=curve.get(k).coords[0];
					 py[1]=curve.get(k).coords[1];
					 px[2]=curve.get(k).coords[0]-(float)Math.sin(curve.get(k).angle);
					 py[2]=curve.get(k).coords[1]-(float)Math.cos(curve.get(k).angle);
				 polyline_p = new PolygonRoi(px, py, Roi.POLYLINE);
				 polyline_p.setPosition(this.get(i).nFrame);
				 polyline_p.setStrokeColor(colors[i%8]);
				 polyline_p.setStrokeWidth(0.0);
				 image_overlay.add(polyline_p);
				  }
				  */
			}
			
		}
		imp.setOverlay(image_overlay);
		imp.updateAndRepaintWindow();
		imp.show();
		return true;
	}
	
	/** function adds curve stored inside object to Overlay **/
	public void toResultsTable()
	{
		int nFrame;
		int nCurve;
		int nPoint;
		CurveSet curveset;
		Curve curve;
		Point point;
		
		ptable = Analyzer.getResultsTable();
		ptable.reset();
		ptable_lock.lock();
		
		//frame cycle
		for (nFrame=1;nFrame<=this.size();nFrame++)
		{
			curveset=this.get(nFrame-1);

			//curve cycle
			for (nCurve=0; nCurve<curveset.size(); nCurve++)
			{
				curve = curveset.get(nCurve);
				//point cycle
				for (nPoint =0;nPoint<curve.size();nPoint++)
				{
					point=curve.get(nPoint);
					ptable.incrementCounter();									
					ptable.addValue("Frame Number", curveset.nFrame);
					ptable.addValue("Contour Number", nCurve+1);
					ptable.addValue("Point Number", nPoint+1);
					ptable.addValue("X_(px)",point.coords[0]);	
					ptable.addValue("Y_(px)",point.coords[1]);
					//ptable.addValue("Frame_Number", nFrame+1);
					ptable.addValue("Angle_of_normal_(radians)", point.angle.val);
					ptable.addValue("Response", point.response);
					ptable.addValue("Amplitude_Fit", point.amp);
					ptable.addValue("SD_fit", point.width);		
					ptable.addValue("BG_fit", point.bg);
					ptable.addValue("IntegrInt_fit", point.integr_int);
					ptable.addValue("R2_fit", point.R2);
	
				}//point cycle end
				
			}//curve cycle end
			
		}//frame cycle end
				
		ptable_lock.unlock();
		ptable.show("Results");
	}
	
	/** function counts all points in all curves of CurveSet
	 * and stores them in nPoints field **/
	public void countPoints()
	{
		int i,j;
	
		nPoints =0;
		for (i=0;i<this.size();i++)
		{
			//cycling through curveset
			for (j=0;j<this.get(i).size();j++)
			{
				nPoints += this.get(i).get(j).size();
			}
		}
	}
	/** returns copies (clones) of all curve points in one array 
	 * (needed to build KD tree) **/
	public Point [] getAllPointsInArray()
	{
		Point [] allPoints;
		Curve curve;
		int i,j,k;
		int nCount =0;
		
		if (!bPointsCounted)
			this.countPoints();
		allPoints = new Point[(int) nPoints];
		
		
		
		//looping through frames		
		for (i=0;i<this.size();i++)
		{
			//looping through curves
			for (j=0;j<this.get(i).size();j++)
			{
				curve=this.get(i).get(j);
				//looping through points
				for (k=0;k<curve.size();k++)				 
				{
					allPoints[nCount]=curve.get(k).clone();
					nCount++;
				}
			}
		}
		
		return allPoints;
		
	}

}

