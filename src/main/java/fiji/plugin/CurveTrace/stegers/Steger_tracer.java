/*  detect-lines, extract lines and their width from images.
    Copyright (C) 1996-1998 Carsten Steger
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2, or (at your option)
    any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

/* 	Changes Made by R. Balasubramanian for incorporating the the detect lines code to incorporate
   	within GRASP (May 10th 1999) */

/*	Port to ImageJ plugin Eugene Katrukha August 2015 */

package fiji.plugin.CurveTrace.stegers;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import org.apache.commons.math3.analysis.MultivariateFunction;

import fiji.plugin.CurveTrace.Fit.OneDGaussian;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.ImageCalculator;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.GaussianBlur;
import ij.process.AutoThresholder.Method;
import ij.process.AutoThresholder;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.process.LUT;
import jaolho.data.lma.LMA;



public class Steger_tracer implements MultivariateFunction
{
	
	/* Constants */
	
	/** The pixel boundaries need to be enlarged slightly since in practice it
	   frequently happens for neighboring pixels a and b that pixel a says a
	   maximum lies within pixel b and vice versa.  This presents no problem since
	   linking algoritm will take care of this. */
	public static double PIXEL_BOUNDARY = 0.6;
	public double SIGMA = 1.8;
	 /** Extract bright lines if true and extract dart lines otherwise*/
	public boolean MODE_LIGHT = true;
	public boolean extend_lines = true;
	public boolean split_lines = true;
	public boolean correct_pos = true;
	public boolean compute_width = true;

	
		

	/** Minimum number of nodes that line can have. Shorter lines will be removed */
	public int nMinNumberOfNodes = 0;
	
	//int nRoiColor;
	/** Maximum angular difference of neighboring line points allowed during
	   linking.  If all feasible neighbors have a larger angular difference the
	   line is stopped at that point. */
	public double MAX_ANGLE_DIFFERENCE = Math.PI/6.0;
	/** Maximum length by which a line is possibly extended in order to find a
	   junction with another line. */
	public double MAX_LINE_EXTENSION = SIGMA*2.5;
	
	/** Maximum search line width. */
	public double MAX_LINE_WIDTH = (2.5*SIGMA);
	
	/** This constant is introduced because for very narrow lines the facet model
	   width detection scheme sometimes extracts the line width too narrow.  Since
	   the correction function has a very steep slope in that area, this will lead
	   to lines of almost zero width, especially since the bilinear interpolation
	   in correct.c will tend to overcorrect.  Therefore it is wise to make the
	   extracted line width slightly larger before correction.  */
	public double LINE_WIDTH_COMPENSATION =1.05;

	/** Minimum line width allowed (used for outlier check in fix_locations()) */
	public double MIN_LINE_WIDTH = 0.1;
	
	/** Maximum contrast allowed (used for outlier check in fix_locations()) */
	public double MAX_CONTRAST=275.0;
	
	
	/** Image width */
	public long width;
	/** Image height */
	public long height;
	
	/** original image minus average value, for MSE estimation**/
	FloatProcessor origMinusMean= null; 
	public FloatProcessor ipMSE;
	public FloatProcessor ipCurveN;
	
	public int nCurveN;
		/* Constant arrays */

		/** This table contains the three appropriate neighbor pixels that the linking
		   algorithm must examine.  It is indexed by the octant the current line
		   angle lies in, e.g., 0 if the angle in degrees lies within [-22.5,22.5]. */
		public static int dirtab[][][] = {
		  { {  1, 0 }, {  1,-1 }, {  1, 1 } },
		  { {  1, 1 }, {  1, 0 }, {  0, 1 } },
		  { {  0, 1 }, {  1, 1 }, { -1, 1 } },
		  { { -1, 1 }, {  0, 1 }, { -1, 0 } },
		  { { -1, 0 }, { -1, 1 }, { -1,-1 } },
		  { { -1,-1 }, { -1, 0 }, {  0,-1 } },
		  { {  0,-1 }, { -1,-1 }, {  1,-1 } },
		  { {  1,-1 }, {  0,-1 }, {  1, 0 } }
		};

		/** This table contains the two neighbor pixels that the linking algorithm
		   should examine and mark as processed in case there are double responses. */
		public static int cleartab[][][] = {
		  { {  0, 1 }, {  0,-1 } },
		  { { -1, 1 }, {  1,-1 } },
		  { { -1, 0 }, {  1, 0 } },
		  { { -1,-1 }, {  1, 1 } },
		  { {  0,-1 }, {  0, 1 } },
		  { {  1,-1 }, { -1, 1 } },
		  { {  1, 0 }, { -1, 0 } },
		  { {  1, 1 }, { -1,-1 } }
		};
		
	ImagePlus impx;
	ResultsTable rt;
	Analyzer anLoc;
	
	//public ResultsTable ptable = new ResultsTable();
	
	/** currently active processor */
	public ImageProcessor ip; 
	/** image calibration**/
	Calibration calIm;
	
	/** Array of all found contours */
	public ArrayList<contour> cont;
	/** Array containing junctions */
	public ArrayList<junction> junc; 
	/** saliency thresholds**/
	double [] salth;
	/** saliency lower threshold */
	public double low; 
	/** saliency upper threshold */
	public double high;
	
	/**Largest eigenvalue of image */
	public double[] eigval; 
	/** eigenvector x component of image */
	public double[] normx; 
	/** eigenvector y component of image */
	public double[] normy; 
	/** x derivative of image */
	public double[] gradx; 
	/** y derivative of image */
	public double[] grady; 
	/** subpixel resolved position of maxima in x*/
	public double[] posx; 
	/** subpixel resolved position of maxima in y*/
	public double[] posy;
		
	/** thresholded array */
	int[] ismax;
	
	
	public double [] threshold_est_min;
	
	public double threshold_est_max;
	
	/** min of saliency map **/
	public double threshold_tmp_ip_min;
	/** max of saliency map **/
	public double threshold_tmp_ip_max;
	
	/** estimated SD of current imageprocessor**/
	public float imgSD;
	
	/** maximum intensity of current imageprocessor**/
	public float imgMAX;
	
	public int nEvalCount=0;
	
	public ArrayList<Double[]> ptable = new ArrayList<Double[]>();
	

	/**   This function links the line points into lines.  The input to this function
	 * are the response of the filter, i.e., the second directional derivative
	 * along (nx[l],ny[l]), contained in eigval[l], and the sub-pixel position of
	 * each line point, contained in (px[l],py[l]).  The parameters low and high
	 * are the hysteresis thresholds for the linking, while width and height are
	 * the dimensions of the five float-images.  The linked lines are returned in
	 * result, and the number of lines detected is returned in num_result. 
	 * */	
	public void compute_contours()
	{
		long area;		
		int i, k, l, pos, nexti;
		int nextpos;
		int it;
		long num_pnt, num_cont, num_junc;
		long x, y;
		long begin, end;
		long indx_max;
		double max;
		long maxx, maxy, nextx, nexty;
		double nx, ny, mx, my;
		double px, py, nextpx, nextpy;
		double alpha, nextalpha, beta, last_beta;
		int octant, last_octant;
		//float tmp;
		boolean nextismax;
		double diff, mindiff, diff1, diff2, dist, mindist;
		double dx, dy;
		double s, t, gx, gy;
		double length, response;
		int num_line, num_add;
		int m=0;
		int j=0;
		double end_angle = 0;
		double end_resp = 0;
		contour tmp_cont;
		boolean add_ext;
		doublepoint closestpnt;
		
		contour_class cls;
		ArrayList<Float> row = new ArrayList<Float>();              
		ArrayList<Float> col = new ArrayList<Float>();             
		ArrayList<Float> angle = new ArrayList<Float>();
		ArrayList<Float> resp = new ArrayList<Float>();  
				
		
		ArrayList<Float> extx;
		ArrayList<Float> exty;
		ArrayList<offset> line;
		ArrayList<Float> trow;              
		ArrayList<Float> tcol;             
		ArrayList<Float> tangle;
		ArrayList<Float> tresp;  
		
		ArrayList<chord> rl = new ArrayList<chord>();
		
		/** The image label contains information on the pixels that have been
	     processed by the linking algorithm. */	
		short label[] = new short[(int) (width*height)];
		/** The image indx is an index into the table of all pixels that possibly
	     could be starting points for new lines.  It is used to quickly determine
	     the next starting point of a line. */
		long indx[] = new long[(int) (width*height)];
		
				
		ArrayList<crossref> cross;
					

		region seg = new region();
		
		
		// Select all pixels that can be starting points for lines. 		  
		seg = threshold(ismax, 2);
		
		// Count the number of possible starting points. 
		area = 0;
		for (i=0; i<seg.num; i++)
			area += seg.rl.get(i).ce - seg.rl.get(i).cb + 1;
		
		// Create the index of possible starting points. 
		cross = new ArrayList<crossref>();					
		rl = seg.rl;
		for (i=0; i<seg.num; i++) 
		{
			x = rl.get(i).r;
			for (y=rl.get(i).cb; y<=rl.get(i).ce; y++) 
			{
				pos = (int) LINCOOR(x,y,width);
				cross.add(new crossref((short)x,(short)y,Math.abs(eigval[pos]), false));
			}
		}
		
		Collections.sort(cross);
		
		for (i=0;i<area;i++)
		    indx[(int) LINCOOR(cross.get(i).x,cross.get(i).y,width)] = i+1;
		
		//reset everything
		num_cont = 0;
		num_junc = 0;
		cont = new ArrayList<contour>();
		junc = new ArrayList<junction>();
		
		// Link lines points. 
		indx_max = 0;
		for (;;) 
		{
			// Contour class unknown at this point; therefore assume both ends free. 
			cls = contour_class.cont_no_junc;
		    while (indx_max < area && cross.get((int) indx_max).done)
		        indx_max++;
		    // Stop if no feasible starting point exists. 
		    if (indx_max == area)
		    	break;
		    max = cross.get((int) indx_max).value;
		    maxx = cross.get((int) indx_max).x;
		    maxy = cross.get((int) indx_max).y;
		    if (max == 0.0)
		    	break;
		    row = new ArrayList<Float>();
		    col = new ArrayList<Float>();
		    resp = new ArrayList<Float>();
		    angle = new ArrayList<Float>();
		    // Add starting point to the line. 
		    num_pnt = 0;
		    pos = (int) LINCOOR(maxx,maxy,width);
		    label[pos] = (short) (num_cont+1);		    
		    if (Math.abs(indx[pos])>0)
		        cross.get((int) (indx[pos]-1)).done = true;
		    row.add((float) posx[pos]);
		    col.add((float) posy[pos]);
		    /*if(DEBUG_show_tracing)
		    {
	    		resolved_p = new OvalRoi(posy[pos]-0.25,posx[pos]-0.25, 0.5, 0.5);
	    		resolved_p.setStrokeColor(Color.YELLOW);
	    		resolved_p.setPosition(nFrame+1);
	    		image_overlay.add(resolved_p);			    
	    		imp.setOverlay(image_overlay);
	    		imp.updateAndRepaintWindow();
	    		imp.show();
		    }*/
		    
		    // Select line direction. 
		    
		    nx = -normy[pos];
		    ny = normx[pos];
		    alpha = Math.atan2(ny,nx);
		    if (alpha < 0.0)
		      alpha += 2.0*Math.PI;
		    if (alpha >= Math.PI)
		      alpha -= Math.PI;
		    octant = (int)(Math.floor(4.0/Math.PI*alpha+0.5))%4;		    
		    // Select normal to the line.  The normal points to the right of the line
		    // as the line is traversed from 0 to num-1.  Since the points are sorted
		    // in reverse order before the second iteration, the first beta actually
		    // has to point to the left of the line! 
		    beta = alpha+Math.PI/2.0;
		    if (beta >= 2.0*Math.PI)
		      beta -= 2.0*Math.PI;
		    angle.add((float) beta);		    
		    resp.add((float)interpolate_response(eigval,maxx,maxy,posx[pos],posy[pos],width,height));
		    num_pnt++;
		    
		    // Mark double responses as processed. 
		    for (i=0;i<2;i++) 
		    {
		      nextx = maxx+cleartab[octant][i][0];
		      nexty = maxy+cleartab[octant][i][1];
		      if (nextx < 0 || nextx >= height || nexty < 0 || nexty >= width)
		        continue;
		      nextpos = (int) LINCOOR(nextx,nexty,width);
		      if (ismax[nextpos] > 0) 
		      {
		        nx = -normy[nextpos];
		        ny = normx[nextpos];
		        nextalpha = Math.atan2(ny,nx);
		        if (nextalpha < 0.0)
		          nextalpha += 2.0*Math.PI;
		        if (nextalpha >= Math.PI)
		          nextalpha -= Math.PI;
		        diff = Math.abs(alpha-nextalpha);
		        if (diff >= Math.PI/2.0)
		          diff = Math.PI-diff;
		        if (diff < MAX_ANGLE_DIFFERENCE) 
		        {
		          label[nextpos] = (short) (num_cont+1);
		          if (Math.abs(indx[nextpos])>0)
		            cross.get((int) (indx[nextpos]-1)).done = true;
		        }
		      }
		    }
		    
		    for (it=1;it<=2;it++) 
		    {
		        if (it == 1) 
		        {
		          // Search along the initial line direction in the first iteration. 
		          x = maxx;
		          y = maxy;
		          pos = (int) LINCOOR(x,y,width);
		          nx = -normy[pos];
		          ny = normx[pos];
		          alpha = Math.atan2(ny,nx);
		          if (alpha < 0.0)
		            alpha += 2.0*Math.PI;
		          if (alpha >= Math.PI)
		            alpha -= Math.PI;
		          last_octant = (int)(Math.floor(4.0/Math.PI*alpha+0.5))%4;
		          last_beta = alpha+Math.PI/2.0;
		          if (last_beta >= 2.0*Math.PI)
		            last_beta -= 2.0*Math.PI;
		        } 
		        else 
		        {
		          // Search in the opposite direction in the second iteration. 
		          x = maxx;
		          y = maxy;
		          pos = (int) LINCOOR(x,y,width);
		          nx = -normy[pos];
		          ny = normx[pos];
		          alpha = Math.atan2(ny,nx);
		          if (alpha < 0.0)
		            alpha += 2.0*Math.PI;
		          if (alpha >= Math.PI)
		            alpha -= Math.PI;
		          last_octant = (int)(Math.floor(4.0/Math.PI*alpha+0.5))%4+4;
		          last_beta = alpha+Math.PI/2.0;
		          if (last_beta >= 2.0*Math.PI)
		            last_beta -= 2.0*Math.PI;
		        }
		        if (it == 2) 
		        {
		          // Sort the points found in the first iteration in reverse. 
                    Collections.reverse(row);
                    Collections.reverse(col);
                    Collections.reverse(angle);
                    Collections.reverse(resp);
		        }		    
		    

		        // Now start adding appropriate neighbors to the line. 
		        for (;;) {
		          pos = (int) LINCOOR(x,y,width);
		          nx = -normy[pos];
		          ny = normx[pos];
		          px = posx[pos];
		          py = posy[pos];
		          // Orient line direction w.r.t. the last line direction. 
		          alpha = Math.atan2(ny,nx);
		          if (alpha < 0.0)
		            alpha += 2.0*Math.PI;
		          if (alpha >= Math.PI)
		            alpha -= Math.PI;
		          octant = (int)(Math.floor(4.0/Math.PI*alpha+0.5))%4;
		          switch(octant) {
		            case 0:
		              if (last_octant >= 3 && last_octant <= 5)
		                octant = 4;
		              break;
		            case 1:
		              if (last_octant >= 4 && last_octant <= 6)
		                octant = 5;
		              break;
		            case 2:
		              if (last_octant >= 4 && last_octant <= 7)
		                octant = 6;
		              break;
		            case 3:
		              if (last_octant == 0 || last_octant >= 6)
		                octant = 7;
		              break;
		          }
		          last_octant = octant;	    

		          // Determine appropriate neighbor. 
		          nextismax = false;
		          nexti = 1;
		          mindiff = Double.MAX_VALUE;
		          for (i=0;i<3;i++) {
		            nextx = x+dirtab[octant][i][0];
		            nexty = y+dirtab[octant][i][1];
		            if (nextx < 0 || nextx >= height || nexty < 0 || nexty >= width)
		              continue;
		            nextpos = (int) LINCOOR(nextx,nexty,width);
		            if (ismax[nextpos] == 0)
		              continue;
		            nextpx = posx[nextpos];
		            nextpy = posy[nextpos];
		            dx = nextpx-px;
		            dy = nextpy-py;
		            dist = Math.sqrt(dx*dx+dy*dy);
		            nx = -normy[nextpos];
		            ny = normx[nextpos];
		            nextalpha = Math.atan2(ny,nx);
		            if (nextalpha < 0.0)
		              nextalpha += 2.0*Math.PI;
		            if (nextalpha >= Math.PI)
		              nextalpha -= Math.PI;
		            diff = Math.abs(alpha-nextalpha);
		            if (diff >= Math.PI/2.0)
		              diff = Math.PI-diff;
		            diff = dist+diff;
		            if (diff < mindiff) {
		              mindiff = diff;
		              nexti = i;
		            }
		            if (Math.abs(ismax[nextpos])>0)
		              nextismax = true;
		          }		          
		          // Mark double responses as processed. 
		          for (i=0;i<2;i++) {
		            nextx = x+cleartab[octant][i][0];
		            nexty = y+cleartab[octant][i][1];
		            if (nextx < 0 || nextx >= height || nexty < 0 || nexty >= width)
		              continue;
		            nextpos = (int) LINCOOR(nextx,nexty,width);
		            if (ismax[nextpos] > 0) {
		              nx = -normy[nextpos];
		              ny = normx[nextpos];
		              nextalpha = Math.atan2(ny,nx);
		              if (nextalpha < 0.0)
		                nextalpha += 2.0*Math.PI;
		              if (nextalpha >= Math.PI)
		                nextalpha -= Math.PI;
		              diff = Math.abs(alpha-nextalpha);
		              if (diff >= Math.PI/2.0)
		                diff = Math.PI-diff;
		              if (diff < MAX_ANGLE_DIFFERENCE) {
		                label[nextpos] = (short) (num_cont+1);
		                if (Math.abs(indx[nextpos])>0)
		                  cross.get((int) (indx[nextpos]-1)).done = true;
		              }
		            }
		          }

		          // Have we found the end of the line? 
		          if (!nextismax)
		            break;
		          // If not, add the neighbor to the line. 
		          x += dirtab[octant][nexti][0];
		          y += dirtab[octant][nexti][1];

		          pos = (int) LINCOOR(x,y,width);		          	          
		          row.add((float) posx[pos]);
		          col.add((float) posy[pos]);
		          
		          // Orient normal to the line direction w.r.t. the last normal. 
		          nx = normx[pos];
		          ny = normy[pos];
		          beta = Math.atan2(ny,nx);
		          if (beta < 0.0)
		            beta += 2.0*Math.PI;
		          if (beta >= Math.PI)
		            beta -= Math.PI;
		          diff1 = Math.abs(beta-last_beta);
		          if (diff1 >= Math.PI)
		            diff1 = 2.0*Math.PI-diff1;
		          diff2 = Math.abs(beta+Math.PI-last_beta);
		          if (diff2 >= Math.PI)
		            diff2 = 2.0*Math.PI-diff2;
		          if (diff1 < diff2) {
		            angle.add((float) beta);
		            last_beta = beta;
		          } else {
		            angle.add((float) (beta+Math.PI));
		            last_beta = beta+Math.PI;
		          }
		          
		          resp.add((float) interpolate_response(eigval,x,y,posx[pos],posy[pos],width,height));
		          num_pnt++;
		          
		          /*
		          if(DEBUG_show_tracing)
		          {
		    		resolved_p = new OvalRoi(posy[pos]-0.25,posx[pos]-0.25, 0.5, 0.5);
		    		resolved_p.setStrokeColor(Color.YELLOW);
		    		resolved_p.setPosition(nFrame+1);
		    		image_overlay.add(resolved_p);			    
		    		imp.setOverlay(image_overlay);
		    		imp.updateAndRepaintWindow();
		    		imp.show();
		          }*/
		          
		          // If the appropriate neighbor is already processed a junction point is found. 
		        if (label[pos] > 0) {
		          
		            // Look for the junction point in the other line. 
		            k = label[pos]-1;
		            if (k == num_cont) {
		              // Line intersects itself. 
		              for (j=0;j<num_pnt-1;j++) {
		                if (row.get(j) == posx[pos] && col.get(j) == posy[pos]) {
		                  if (j == 0) {
		                    // Contour is closed. 
		                    cls = contour_class.cont_closed;
		                    Collections.reverse(row);
		                    Collections.reverse(col);
		                    Collections.reverse(angle);
		                    Collections.reverse(resp);
		                    it = 2;		        	
		                  } else {
		                      if (it == 2) {
		                        // Determine contour class. 
		                        if (cls == contour_class.cont_start_junc)
		                          cls = contour_class.cont_both_junc;
		                        else
		                          cls = contour_class.cont_end_junc;
		                        // Index j is the correct index. 
		                        junc.add(new junction(num_cont,num_cont,j,(float)posx[pos],(float)posy[pos]));
		                        num_junc++;
		                      } else {
		                        // Determine contour class. 
		                        cls = contour_class.cont_start_junc;
		                        // Index num_pnt-1-j is the correct index since the line
		                         //  is going to be sorted in reverse. 
		                        junc.add(new junction(num_cont,num_cont,num_pnt-1-j,(float)posx[pos],(float)posy[pos]));

		                        num_junc++;
		                      }
		                    }
		                    break;
		                  }
		                }		        	
		              // Mark this case as being processed for the algorithm below. 
		              j = -1;
		            } else {
		              for (j=0;j<cont.get(k).num;j++) 
		              {
		            	mindist = Math.sqrt(Math.pow(cont.get(k).row.get(j) -posx[pos],2)+Math.pow(cont.get(k).col.get(j) - posy[pos],2));  
		                if (mindist<0.00001)
		                  break;
		              }
		              // If no point can be found on the other line a double response
		              //   must have occured.  In this case, find the nearest point on
		              //   the other line and add it to the current line. 
		              if (j == cont.get(k).num) {
		                mindist = Double.MAX_VALUE;
		                j = -1;
		                for (l=0;l<cont.get(k).num;l++) {
		                  dx = posx[pos]-cont.get(k).row.get(l);
		                  dy = posy[pos]-cont.get(k).col.get(l);
		                  dist = Math.sqrt(dx*dx+dy*dy);
		                  if (dist < mindist) {
		                    mindist = dist;
		                    j = l;
		                  }
		                }
		                // Add the point with index j to the current line. 

		                row.add(cont.get(k).row.get(j));
		                col.add(cont.get(k).col.get(j));
		                beta = cont.get(k).angle.get(j);
		                if (beta >= Math.PI)
		                  beta -= Math.PI;
		                diff1 = Math.abs(beta-last_beta);
		                if (diff1 >= Math.PI)
		                  diff1 = 2.0*Math.PI-diff1;
		                diff2 = Math.abs(beta+Math.PI-last_beta);
		                if (diff2 >= Math.PI)
		                  diff2 = 2.0*Math.PI-diff2;
		                if (diff1 < diff2)
		                  angle.add((float) beta);
		                else
		                  angle.add((float) (beta+Math.PI));
		                resp.add(cont.get(k).response.get(j));
		                num_pnt++;
		              }
		            }

		              // Add the junction point only if it is not one of the other line's
		              // endpoints. 
		            if (j > 0 && j < cont.get(k).num-1) {
		              // Determine contour class. 
		              if (it == 1)
		                cls = contour_class.cont_start_junc;
		              else if (cls == contour_class.cont_start_junc)
		                cls = contour_class.cont_both_junc;
		              else
		                cls = contour_class.cont_end_junc;
		              // Add the new junction. 
		              junc.add(new junction(k, num_cont,j,row.get((int) (num_pnt-1)), col.get((int) (num_pnt-1))));
		              num_junc++;
		            }
		            break;
		          }
		          label[pos] = (short) (num_cont+1);
		          if (Math.abs(indx[pos])>0)
		            cross.get((int) (indx[pos]-1)).done = true;
		        }
		      }
		    
		    
		    if (num_pnt > 1) {
		        // Only add lines with at least two points. 	        
		        cont.add(new contour(num_pnt,row,col,angle,resp,cls));		        
		        num_cont++;
		      } else {
		        // Delete the point from the label image; we can use maxx and maxy
		        // as the coordinates in the label image in this case. 
		        for (i=-1;i<=1;i++) {
		          for (j=-1;j<=1;j++) {
		            pos = (int) LINCOOR(BR(maxx+i),BC(maxy+j),width);
		            if (label[pos] == num_cont+1)
		              label[pos] = 0;
		          }
		        }
		      }
		    }

		  // Now try to extend the lines at their ends to find additional junctions. 
		  if (extend_lines) {
		    // Sign by which the gradient has to be multiplied below. 
		    if (MODE_LIGHT)
		      s = 1;
		    else
		      s = -1;
		    length = MAX_LINE_EXTENSION;
		    

		    for (i=0; i<num_cont; i++) {		            
			    
		    	tmp_cont = cont.get(i);
		        num_pnt = tmp_cont.num;
		        if (num_pnt == 1)
		          continue;
		        if (tmp_cont.cont_class == contour_class.cont_closed)
		          continue;
		        trow = tmp_cont.row;
		        tcol = tmp_cont.col;
		        tangle = tmp_cont.angle;
		        tresp = tmp_cont.response;
		        // Check both ends of the line (it==-1: start, it==1: end). 
		        for (it=-1; it<=1; it+=2) {

		        // Determine the direction of the search line.  This is done by using 
		          // the normal to the line (angle).  Since this normal may point to
		          // the left of the line (see below) we have to check for this case by
		          // comparing the normal to the direction of the line at its respective
		          // end point. 
		            if (it==-1) {
		                // Start point of the line. 
		                if (tmp_cont.cont_class == contour_class.cont_start_junc ||
		                    tmp_cont.cont_class == contour_class.cont_both_junc)
		                  continue;
		                dx = trow.get(1)-trow.get(0);
		                dy = tcol.get(1)-tcol.get(0);
		                alpha = tangle.get(0);
		                nx = Math.cos(alpha);
		                ny = Math.sin(alpha);
		                if (nx*dy-ny*dx < 0) {
		                  // Turn the normal by +90 degrees. 
		                  mx = -ny;
		                  my = nx;
		                } else {
		                  // Turn the normal by -90 degrees. 
		                  mx = ny;
		                  my = -nx;
		                }
		                px = trow.get(0);
		                py = tcol.get(0);
		                response = tresp.get(0);
		              } else {		    	
		                  // End point of the line. 
		                  if (tmp_cont.cont_class == contour_class.cont_end_junc ||
		                      tmp_cont.cont_class == contour_class.cont_both_junc)
		                    continue;
		                  dx = trow.get((int) (num_pnt-1))-trow.get((int) (num_pnt-2));
		                  dy = tcol.get((int) (num_pnt-1))-tcol.get((int) (num_pnt-2));
		                  alpha = tangle.get((int) (num_pnt-1));
		                  nx = Math.cos(alpha);
		                  ny = Math.sin(alpha);
		                  if (nx*dy-ny*dx < 0) {
		                    // Turn the normal by -90 degrees. 
		                    mx = ny;
		                    my = -nx;
		                  } else {
		                    // Turn the normal by +90 degrees. 
		                    mx = -ny;
		                    my = nx;
		                  }
		                  px = trow.get((int) (num_pnt-1));
		                  py = tcol.get((int) (num_pnt-1));
		                  response = tresp.get((int) (num_pnt-1));
		                }	
		            // Determine the current pixel and calculate the pixels on the search line. 
		         x = (long) Math.floor(px+0.5);
		         y = (long) Math.floor(py+0.5);
		         dx = px-x;
		         dy = py-y;
		         line = bresenham(mx,my,dx,dy,length);
		         num_line = line.size();
		         // Now determine whether we can go only uphill (bright lines) or
		         // downhill (dark lines) until we hit another line. 
			     
				 extx = new ArrayList<Float>();
				 exty = new ArrayList<Float>();
			     num_add = 0;
		         add_ext = false;
		         for (k=0; k<num_line; k++) {
		           nextx = x+line.get(k).x;
		           nexty = y+line.get(k).y;
		           closestpnt = closest_point(px,py,mx,my,(double)nextx,(double)nexty);
                   nextpx = closestpnt.cx;
                   nextpy = closestpnt.cy;
                   t = closestpnt.t;
                   
		           // Ignore points before or less than half a pixel away from the
		           // true end point of the line. 
		           if (t <= 0.5)
		             continue;
		           // Stop if the gradient can't be interpolated any more or if the
		           // next point lies outside the image. 
		           if (nextpx < 0 || nextpy < 0 ||
		               nextpx >= height-1 || nextpy >= width-1 ||
		               nextx < 0 || nexty < 0 ||
		               nextx >= height || nexty >= width)
		             break;
		           closestpnt = interpolate_gradient(gradx,grady,nextpx,nextpy,(int)width);
		           gx = closestpnt.cx;
		           gy = closestpnt.cy;
		           // Stop if we can't go uphill anymore.  This is determined by the
		           // dot product of the line direction and the gradient.  If it is
		           // smaller than 0 we go downhill (reverse for dark lines). 
		           nextpos = (int) LINCOOR(nextx,nexty,width);		            	  
		           if (s*(mx*gx+my*gy) < 0 && label[nextpos] == 0)
		               break;
		             // Have we hit another line? 
		             if (label[nextpos] > 0) {
		               m = label[nextpos]-1;
		               // Search for the junction point on the other line. 
		               mindist = Double.MAX_VALUE;
		               j = -1;
		               for (l=0; l<cont.get(m).num; l++) {
		                 dx = nextpx-cont.get(m).row.get(l);
		                 dy = nextpy-cont.get(m).col.get(l);
		                 dist = Math.sqrt(dx*dx+dy*dy);
		                 if (dist < mindist) {
		                   mindist = dist;
		                   j = l;
		                 }
		               }		            	  
		               // This should not happen...  But better safe than sorry... 
		               if (mindist > 3.0)
		                 break;
		               extx.add(cont.get(m).row.get(j));
		               exty.add(cont.get(m).col.get(j));
		               end_resp = cont.get(m).response.get(j);
		               end_angle = cont.get(m).angle.get(j);
		               beta = end_angle;
		               if (beta >= Math.PI)
		                 beta -= Math.PI;
		               diff1 = Math.abs(beta-alpha);
		               if (diff1 >= Math.PI)
		                 diff1 = 2.0*Math.PI-diff1;
		               diff2 = Math.abs(beta+Math.PI-alpha);
		               if (diff2 >= Math.PI)
		                 diff2 = 2.0*Math.PI-diff2;
		               if (diff1 < diff2)
		                 end_angle = beta;
		               else
		                 end_angle = beta+Math.PI;
		               num_add++;
		               /*
		               if(DEBUG_show_extensions)
		               {
			    		resolved_p = new OvalRoi(cont.get(m).col.get(j)-0.25,cont.get(m).row.get(j)-0.25, 0.5, 0.5);
			    		resolved_p.setStrokeColor(Color.GREEN);
			    		resolved_p.setPosition(nFrame+1);
			    		image_overlay.add(resolved_p);			    
			    		imp.setOverlay(image_overlay);
			    		imp.updateAndRepaintWindow();
			    		imp.show();
		               }*/
		               add_ext = true;
		               break;
		             } else {
		               extx.add((float) nextpx);
		               exty.add((float) nextpy);
		               /*
		               if(DEBUG_show_extensions)
		               {
			    		resolved_p = new OvalRoi(nextpy-0.25,nextpx-0.25, 0.5, 0.5);
			    		resolved_p.setStrokeColor(Color.GREEN);
			    		resolved_p.setPosition(nFrame+1);
			    		image_overlay.add(resolved_p);			    
			    		imp.setOverlay(image_overlay);
			    		imp.updateAndRepaintWindow();
			    		imp.show();
		               }*/
		               num_add++;
		             }
		           }
		         if (add_ext) {
		             // Make room for the new points. 
		             num_pnt += num_add;

		             tmp_cont.row = trow;
		             tmp_cont.col = tcol;
		             tmp_cont.angle = tangle;
		             tmp_cont.response = tresp;
		             tmp_cont.num = num_pnt;
		             if (it == -1) {
		               // Move points on the line up num_add places. 
	
		               // Insert points at the beginning of the line. 
		               for (k=0; k<num_add; k++) 
		               {
		            	 //cause order of insertion is different
			             trow.add(0,extx.get(k));
			             tcol.add(0,exty.get(k));

  	            	     tangle.add(0,(float) alpha);
		                 tresp.add(0,(float) response);
		               }
		               
		               tangle.set(0,(float)end_angle);
		               tresp.set(0,(float)end_resp);
		               // Adapt indices of the previously found junctions. 
		               for (k=0; k<num_junc; k++) {
		                 if (junc.get(k).cont1 == i)
		                   junc.get(k).pos += num_add;
		               }
		             } else {
		               // Insert points at the end of the line. 
		               for (k=0; k<num_add; k++) {

		            	 trow.add(extx.get(k));
			             tcol.add(exty.get(k));
			             tangle.add((float) alpha);
			             tresp.add((float) response);
		               }
		               tangle.set( (int)(num_pnt-1),(float)end_angle);
		               tresp.set((int)(num_pnt-1),(float)end_resp);
		             }		         

		             // Add the junction point only if it is not one of the other line's
		              //  endpoints. 
		             if (j > 0 && j < cont.get(m).num-1) {
		               if (it == -1) {
		                 if (tmp_cont.cont_class == contour_class.cont_end_junc)
		                   tmp_cont.cont_class = contour_class.cont_both_junc;
		                 else
		                   tmp_cont.cont_class = contour_class.cont_start_junc;
		               } else {
		                 if (tmp_cont.cont_class == contour_class.cont_start_junc)
		                   tmp_cont.cont_class = contour_class.cont_both_junc;
		                 else
		                   tmp_cont.cont_class = contour_class.cont_end_junc;
		               }
		               junc.add(new junction());
		               junc.get((int) num_junc).cont1 = m;
		               junc.get((int) num_junc).cont2 = i;
		               junc.get((int) num_junc).pos = j;
		               if (it == -1) {
		                 junc.get((int) num_junc).x = trow.get(0);
		                 junc.get((int) num_junc).y = tcol.get(0);
		               } else {
		                 junc.get((int) num_junc).x = trow.get((int) (num_pnt-1));
		                 junc.get((int) num_junc).y = tcol.get((int) (num_pnt-1));
		               }
		               num_junc++;
		             }
		           }
		         }
		       }		         
		  }

		  // Done with linking.  Now split the lines at the junction points. 
		  
		  Collections.sort(junc);
		  if(split_lines)
		  {
		  
			  for (i=0;i<num_junc;i+=k) {
			    j = (int) junc.get(i).cont1;
			    tmp_cont = cont.get(j);
			    num_pnt = tmp_cont.num;
			    // Count how often line j needs to be split. 
			    for (k=0;i+k<num_junc;k++)
			    {
			    	if(junc.get(i+k).cont1 != j)
			    		break;
			    }
			    if (k == 1 &&
			        tmp_cont.row.get(0) == tmp_cont.row.get((int) (num_pnt-1)) &&
			        tmp_cont.col.get(0) == tmp_cont.col.get((int) (num_pnt-1))) {
			      // If only one junction point is found and the line is closed it only
			       // needs to be rearranged cyclically, but not split. 
			      begin = junc.get(i).pos;
			      trow = tmp_cont.row;
			      tcol = tmp_cont.col;
			      tangle = tmp_cont.angle;
			      tresp = tmp_cont.response;	    
			      //tmp_cont->row = xcalloc(num_pnt,sizeof(float));
			      //tmp_cont->col = xcalloc(num_pnt,sizeof(float));
			      //tmp_cont->angle = xcalloc(num_pnt,sizeof(float));
			      //tmp_cont->response = xcalloc(num_pnt,sizeof(float));
			      for (l=0;l<num_pnt;l++) {
			        pos = (int) (begin+l);
			        // Skip starting point so that it is not added twice. 
			        if (pos >= num_pnt)
			          pos = (int) (begin+l-num_pnt+1);
			        tmp_cont.row.set(l,trow.get(pos));
			        tmp_cont.col.set(l,tcol.get(pos));
			        tmp_cont.angle.set(l,tangle.get(pos));
			        tmp_cont.response.set(l,tresp.get(pos));
			      }
			      // Modify contour class. 
			      tmp_cont.cont_class = contour_class.cont_both_junc;

			    } else {		    
			    	// Otherwise the line has to be split. 
			        for (l=0;l<=k;l++) {
			          if (l == 0)
			            begin = 0;
			          else
			            begin = junc.get(i+l-1).pos;
			          if (l==k)
			            end = tmp_cont.num-1;
			          else
			            end = junc.get(i+l).pos;
			          num_pnt = end-begin+1;
			          if (num_pnt == 1 && k > 1) {
			            // Do not add one point segments. 
			            continue;
			          }
			          
			          cont.add(new contour());
			          //till begin+num_pnt, since last index is exclusive in subList function
			          cont.get((int) num_cont).row = new ArrayList<Float>(tmp_cont.row.subList((int)begin,(int)( begin+num_pnt)));
			          cont.get((int) num_cont).col = new ArrayList<Float>(tmp_cont.col.subList((int)begin,(int)( begin+num_pnt)));
			          cont.get((int) num_cont).angle = new ArrayList<Float>(tmp_cont.angle.subList((int)begin,(int)( begin+num_pnt)));
			          cont.get((int) num_cont).response = new ArrayList<Float>(tmp_cont.response.subList((int)begin,(int)( begin+num_pnt)));
			          cont.get((int) num_cont).num = num_pnt;

			          // Modify contour class. 
			          if (l == 0) {
			            if (tmp_cont.cont_class == contour_class.cont_start_junc ||
			                tmp_cont.cont_class == contour_class.cont_both_junc)
			              cont.get((int)num_cont).cont_class = contour_class.cont_both_junc;
			            else
			              cont.get((int)num_cont).cont_class = contour_class.cont_end_junc;
			          } else if (l == k) {
			            if (tmp_cont.cont_class == contour_class.cont_end_junc ||
			                tmp_cont.cont_class == contour_class.cont_both_junc)
			              cont.get((int)num_cont).cont_class = contour_class.cont_both_junc;
			            else
			              cont.get((int)num_cont).cont_class = contour_class.cont_start_junc;
			          } else {
			            cont.get((int)num_cont).cont_class = contour_class.cont_both_junc;
			          }
			          num_cont++;
			        }
			        num_cont= num_cont -1;
			        cont.set(j,cont.get((int) num_cont)); //*???? WTF
			        cont.remove((int) num_cont);
			        
			        tmp_cont = null;
			      }
			    }		    

		  }
		  		 
		  		    	
		  
		  // Finally, check whether all angles point to the right of the line. 
		  for (i=0; i<num_cont; i++) 
		  {
		    tmp_cont = cont.get(i);
		    num_pnt = tmp_cont.num;
		    trow = tmp_cont.row;
		    tcol = tmp_cont.col;
		    tangle = tmp_cont.angle;
		    
		    
		    
		    	
		    // One point of the contour is enough to determine the orientation. 
		    k = (int) ((num_pnt-1)/2);
		    // The next few lines are ok because lines have at least two points. 
		    dx = trow.get(k+1)-trow.get(k);
		    dy = tcol.get(k+1)-tcol.get(k);
		    nx = Math.cos(tangle.get(k));
		    ny = Math.sin(tangle.get(k));
		    // If the angles point to the left of the line they have to be adapted.
		    //  The orientation is determined by looking at the z-component of the
		    //  cross-product of (dx,dy,0) and (nx,ny,0). 
		    if (nx*dy-ny*dx < 0) 
		    {
		      for (j=0; j<num_pnt; j++) 
		      {
		        tangle.set(j,(float) (tangle.get(j)+Math.PI));
		        if (tangle.get(j) >= 2*Math.PI)
		          tangle.set(j,(float) (tangle.get(j)-2*Math.PI));
		      }
		    }
		  }
		 
		  //Remove lines with number of points less than threshold
		  if(nMinNumberOfNodes>0)
		  {
		   
			  for (i=0; i<num_cont; i++) 
			  {
				  if(cont.get(i).row.size()<nMinNumberOfNodes)
				  {
					  cont.remove(i);
					  num_cont--;
					  i--;					  
				  }				  
			  }
		  }
		    
		
	}
	
	/** For each point in the image determine whether there is a local maximum of
	 * the second directional derivative in the direction (nx[l],ny[l]) within the
	 * pixels's boundaries.  If so, set ismax[l] to 2 if the eigenvalue ev[l] is
	 * larger than high, to 1 if ev[l] is larger than low, and to 0 otherwise.
	 * Furthermore, put the sub-pixel position of the maximum into (px[l],py[l]).
	 * The parameter mode determines whether maxima (dark lines points) or minima
	 * (bright line points) should be selected.  The partial derivatives of the
	 * image are input as ku[].*/ 
	public void compute_line_points()
	{
		double[][] line_points;
		position poscalc = new position(ip, SIGMA);
		
		//get line points		
		line_points = poscalc.compute_line_points(MODE_LIGHT);
		
		
		//old version
		//line_points = DerivativeOfGaussian.get_line_points_stegers(ip, SIGMA, MODE_LIGHT);
		eigval = line_points[0];
		normx  = line_points[1];
		normy  = line_points[2];
		gradx  = line_points[3];
		grady  = line_points[4];
		posx   = line_points[5];
		posy   = line_points[6];
	}
	

	/** Threshold an image above min and return the result as a run-length encoded
	region in out. */
	public region threshold(int image[],long min)
	{
		
	region outr = new region();
	
	long   grey;
	long   r,c,l,num;
	boolean   inside;
	ArrayList<chord> rl = new ArrayList<chord>();
	
	inside = false;
	num = 0;
	rl.add(new chord());
	
	for (r=0; r<height; r++) {
	    for (c=0; c<width; c++) {
	      l = LINCOOR(r,c,width);
	      grey = (long) image[(int) l];
	      if (grey >= min) {
	        if (!inside) {
	        	inside = true;
	            rl.get((int) num).r = (short) r;
	            rl.get((int) num).cb = (short) c;
	        }
	      } else {
	        if (inside) {
	        	inside = false;
	            rl.get((int) num).ce =(short)(c - 1);
	            num++;
	            rl.add(new chord());
	        }
	      }
	    }
	    if (inside) {
	    	   inside = false;
	    	   rl.get((int)num).ce = (short)(width-1);
	    	   num++;
	    	   rl.add(new chord());
	
	    }
	  }
	
	outr.rl = new ArrayList<chord>(rl);
	outr.num = num;
	
	return outr;
	
	}

	
	/** 
	 * Get valid line points (maxima within boundary of pixel, plus higher than low and high)
	 * */
	public void getismax(double [] arg)
	{
		double low=arg[0];
		double high=arg[1];
		
		ismax = new int [(int)(width*height)];
		int l;
		double val;
		for(int py = 0; py < height; ++py)
		{
			for(int px = 0; px < width; ++px)
			{
				l = (int) LINCOOR(py,px,width);

				val = eigval[l];

				if(!Double.isNaN(posx[l]))
				{
					if(Math.abs(posy[l]-px)<=PIXEL_BOUNDARY && Math.abs(posx[l]-py)<=PIXEL_BOUNDARY)
					{
						if(val>=low)
						{
							if(val>=high)
								ismax[l]=2;
							else
								ismax[l]=1;
						}
					
					}
					else
						ismax[l]=0;
				
				}
				else
					ismax[l]=0;
			}
		
		}
		return;

	}
	
	/** Modified Bresenham algorithm.  It returns in line all pixels that are
	   intersected by a half line less than length away from the point (px,py)
	   along the direction (nx,ny).  The point (px,py) must lie within the pixel
	   of the origin, i.e., fabs(px) <= 0.5 and fabs(py) <= 0.5. 
	   */
	public ArrayList<offset> bresenham(double nx, double ny, double px, double py, double length)
	{
		
	  ArrayList<offset> line = new ArrayList<offset>(); 
	  int i, x, y, s1, s2, xchg, maxit;
	  double e, dx, dy, t;

	  x = 0;
	  y = 0;
	  dx = Math.abs(nx);
	  dy = Math.abs(ny);
	  s1 = (int) Math.signum(nx);
	  s2 = (int) Math.signum(ny);
	  px *= s1;
	  py *= s2;
	  if (dy > dx) {
	    t = dx;
	    dx = dy;
	    dy = t;
	    t = px;
	    px = py;
	    py = t;
	    xchg = 1;
	  } else {
	    xchg = 0;
	  }
	  maxit = (int) Math.ceil(length*dx);
	  e = (0.5-px)*dy/dx-(0.5-py);
	  for (i=0; i<=maxit; i++) {
	    line.add(new offset(x,y));		  	    
	  
	    while (e >= -1e-8) {
	      if (Math.abs(xchg)>0) x += s1;
	      else y += s2;
	      e--;
	      if (e > -1) {
	    	line.add(new offset(x,y));
	   
	      }
	    }
	    if (Math.abs(xchg)>0) y += s2;
	    else x += s1;
	    e += dy/dx;
	  }
//	  *num_points = n;
	  return line;
	}
	 /** Translate row and column coordinates of an image into an index into its
	    one-dimensional array. */
	public static long LINCOOR(long row, long col,long width){ 
		return (long)((row)*(width)+(col));				
	}
	
	/** Mirror the row coordinate at the borders of the image; height must be a
	   defined variable in the calling function containing the image height. */
	public long BR(long row) {
		if (row < 0)  
			return -1*row; 				
		else
		{
			if (row >= height) 
				return height - row + height - 2; 
			else
				return row;
		}
		
	}

	/** Mirror the column coordinate at the borders of the image; width must be a
	   defined variable in the calling function containing the image width. */
	public long BC(long col) { 
	 if (col < 0) 
		 return (-1*col);
	 else
	 {
		 if (col >= width)
			return  width - (col) + width - 2;
		 else 
	     	return col;
		 }
	}
	/** Compute the response of the operator with sub-pixel accuracy by using the
	 *  facet model to interpolate the pixel accurate responses. */
	public double interpolate_response(double resp[],long x, long y, double px, double py, long width, long height)
	{
	  double i1, i2, i3, i4, i5, i6, i7, i8, i9;
	  double t1, t2, t3, t4, t5, t6;
	  double d, dr, dc, drr, drc, dcc;
	  double xx, yy;

	  i1 = resp[(int) LINCOOR(BR(x-1),BC(y-1),width)];
	  i2 = resp[(int) LINCOOR(BR(x-1),y,width)];
	  i3 = resp[(int) LINCOOR(BR(x-1),BC(y+1),width)];
	  i4 = resp[(int) LINCOOR(x,BC(y-1),width)];
	  i5 = resp[(int) LINCOOR(x,y,width)];
	  i6 = resp[(int) LINCOOR(x,BC(y+1),width)];
	  i7 = resp[(int) LINCOOR(BR(x+1),BC(y-1),width)];
	  i8 = resp[(int) LINCOOR(BR(x+1),y,width)];
	  i9 = resp[(int) LINCOOR(BR(x+1),BC(y+1),width)];
	  t1 = i1+i2+i3;
	  t2 = i4+i5+i6;
	  t3 = i7+i8+i9;
	  t4 = i1+i4+i7;
	  t5 = i2+i5+i8;
	  t6 = i3+i6+i9;
	  d = (-i1+2*i2-i3+2*i4+5*i5+2*i6-i7+2*i8-i9)/9;
	  dr = (t3-t1)/6;
	  dc = (t6-t4)/6;
	  drr = (t1-2*t2+t3)/6;
	  dcc = (t4-2*t5+t6)/6;
	  drc = (i1-i3-i7+i9)/4;
	  xx = px-x;
	  yy = py-y;
	  return d+xx*dr+yy*dc+xx*xx*drr+xx*yy*drc+yy*yy*dcc;
	}

	/** Calculate the closest point to (px,py) on the line (lx,ly) + t*(dx,dy)
	 *  and return the result in (cx,cy), plus the parameter in t. */
	public doublepoint closest_point(double lx,double ly,double dx,double dy,double px,double py)		  
	{
		
	  doublepoint cp = new doublepoint(0,0,0); 
	  double mx, my, den, nom, tt;

	  mx = px-lx;
	  my = py-ly;
	  den = dx*dx+dy*dy;
	  nom = mx*dx+my*dy;
	  if (den != 0) 
		  tt = nom/den;
	  else 
		  tt = 0;
	  cp.cx = lx+tt*dx;
	  cp.cy = ly+tt*dy;
	  cp.t = tt;
	  return cp;
	  
	}
	
	/** 
	 * Interpolate the gradient of the gradient images gradx and grady with width
	 * width at the point (px,py) using linear interpolation, and return the
	 * result in (gx,gy) (doublepoint format)
	 * */
	public doublepoint interpolate_gradient(double gradx[], double grady[], double px, double py, int width)		  
	{
	  doublepoint gradxy = new doublepoint(0,0,0); 
	  long   gix, giy, gpos;
	  double gfx, gfy, gx1, gy1, gx2, gy2, gx3, gy3, gx4, gy4;

	  gix = (long) Math.floor(px);
	  giy = (long) Math.floor(py);
	  gfx = px % 1.0; //check whether it works as promised
	  gfy = py % 1.0;
	  gpos = LINCOOR(gix,giy,width);
	  gx1 = gradx[(int) gpos];
	  gy1 = grady[(int) gpos];
	  gpos = LINCOOR(gix+1,giy,width);
	  gx2 = gradx[(int) gpos];
	  gy2 = grady[(int) gpos];
	  gpos = LINCOOR(gix,giy+1,width);
	  gx3 = gradx[(int) gpos];
	  gy3 = grady[(int) gpos];
	  gpos = LINCOOR(gix+1,giy+1,width);
	  gx4 = gradx[(int) gpos];
	  gy4 = grady[(int) gpos];
	  gradxy.cx = (1-gfy)*((1-gfx)*gx1+gfx*gx2)+gfy*((1-gfx)*gx3+gfx*gx4);
	  gradxy.cy = (1-gfy)*((1-gfx)*gy1+gfx*gy2)+gfy*((1-gfx)*gy3+gfx*gy4);
	  return gradxy;
	}

	
	/** Provides user with image of all valid points with saliency value at it
	 * and dialogue to pick thresholds for line detection
	 * 
	 * */
	public boolean get_threshold_levels()
	{
		int l,i;
		double val;
		
		
		
		//new image processor with saliency map 
		//final ImageProcessor threshold_tmp_ip = new FloatProcessor((int)width, (int)height);
		final ImageProcessor threshold_tmp_ip = new FloatProcessor((int)width, (int)height);
		for(int py = 0; py < height; ++py)
		{
			for(int px = 0; px < width; ++px)
			{
				l = (int) LINCOOR(py,px,width);
				val = eigval[l];

				
				if(!Double.isNaN(posx[l]))
				{
					if(Math.abs(posy[l]-px)<=PIXEL_BOUNDARY && Math.abs(posx[l]-py)<=PIXEL_BOUNDARY && val>0)
					{
						threshold_tmp_ip.setf(px,py,(float) val);
					}				
					else
						threshold_tmp_ip.setf(px,py,(float) 0);
				}
			}		
		}//for cycle
		
		
		threshold_tmp_ip.resetMinAndMax();
		threshold_tmp_ip_min = threshold_tmp_ip.getMin();
		threshold_tmp_ip_max = threshold_tmp_ip.getMax();
		AutoThresholder autoT=new AutoThresholder();
		ImageStatistics stats = threshold_tmp_ip.getStats();
		threshold_est_min = new double[4];
		for (i=0;i<4;i++)
		{
			switch (i)
			{
				case 0:
					threshold_est_min[0]=threshold_tmp_ip_min+stats.binSize*(autoT.getThreshold(Method.Huang, stats.histogram)+1);
					break;
				case 1:
					threshold_est_min[1]=threshold_tmp_ip_min+stats.binSize*(autoT.getThreshold(Method.Default, stats.histogram)+1);
					break;
				case 2:
					threshold_est_min[2]=threshold_tmp_ip_min+stats.binSize*(autoT.getThreshold(Method.Otsu, stats.histogram)+1);
					break;
				case 3:
					threshold_est_min[3]=threshold_tmp_ip_min+stats.binSize*(autoT.getThreshold(Method.MinError, stats.histogram)+1);
					break;
			}
		}
		

		
		threshold_est_max=threshold_tmp_ip_max;
		//threshold_est_max=getSaturationLevel(threshold_tmp_ip, 0.35);
		final double range = threshold_tmp_ip_max - threshold_tmp_ip_min;
		//final ImagePlus threshold_imp_fl = new ImagePlus("Threshold map float", threshold_tmp_ip);
		//threshold_imp_fl.show();
		//final ImageProcessor threshold_ip = threshold_tmp_ip.convertToByteProcessor(); // final to make accessible in anonymous inner class
		final ImagePlus threshold_imp = new ImagePlus("Threshold map", threshold_tmp_ip); // final to make accessible in anonymous inner class		
		final double threshold_map_scale_factor = (range) / 256;
		threshold_imp.resetDisplayRange();
		threshold_imp.show();
		
		
		
	
		ThresholdSal thresholds_gd = new ThresholdSal(); // RSLV: make NonBlockingGenericDialog();
	
		thresholds_gd.addMessage("Upper threshold (green) defines where line tracing will start");
		thresholds_gd.addSlider("Upper threshold", threshold_tmp_ip_min, threshold_tmp_ip_max, Prefs.get("stegers_source.high", threshold_tmp_ip_min + 0.75*range)); // RSLV: limit?
		thresholds_gd.addMessage("Lower threshold (blue) defines the area where line tracing will occur");
		thresholds_gd.addSlider("Lower threshold", threshold_tmp_ip_min, threshold_tmp_ip_max, Prefs.get("stegers_source.low",threshold_tmp_ip_min + 0.25*range)); // RSLV: limit?
		
		DialogListener dlg_listener = new DialogListener(){
			@Override
			public boolean dialogItemChanged(GenericDialog gd, java.awt.AWTEvent e)
			{
				// get threshold parameters
				double ut = gd.getNextNumber()-threshold_tmp_ip_min;
				double lt = gd.getNextNumber()-threshold_tmp_ip_min;
							
				// use luts to highlight regions
				// NOTE: LUTs are 8-bit, so map is only an approximation
				int but = (int)(ut / threshold_map_scale_factor);
				int blt = (int)(lt / threshold_map_scale_factor);
				
				byte[] threshold_lut_r = new byte[256];
				byte[] threshold_lut_g = new byte[256];
				byte[] threshold_lut_b = new byte[256];
				
				// create lut
				for(int i = 0; i < 256; ++i)
				{
					if(i < blt)
					{
						// retain gray scale
						threshold_lut_r[i] = (byte)i;
						threshold_lut_g[i] = (byte)i;
						threshold_lut_b[i] = (byte)i;
					}
					else if(i < but)
					{
						// set to lower threshold colour
						threshold_lut_r[i] = (byte)0;
						threshold_lut_g[i] = (byte)0;
						threshold_lut_b[i] = (byte)255;
					}
					else
					{
						// set to upper threshold colour
						threshold_lut_r[i] = (byte)0;
						threshold_lut_g[i] = (byte)255;
						threshold_lut_b[i] = (byte)0;
					}
				}
				LUT threshold_lut = new LUT(threshold_lut_r, threshold_lut_g, threshold_lut_b);
				
				// set LUT and update image
				//threshold_ip.setLut(threshold_lut);
				threshold_tmp_ip.setLut(threshold_lut);
				threshold_imp.updateAndDraw();
				
				return true; // true == accepted values, false == incorrect values
			}
		};
		thresholds_gd.addDialogListener(dlg_listener);
		dlg_listener.dialogItemChanged(thresholds_gd, null); // force update of lut
		
		thresholds_gd.setOKLabel("Continue tracing");
		thresholds_gd.setCancelLabel("Cancel tracing");
		
		// focus on threshold imp window, (sometimes get lost behind other debug image windows)
		threshold_imp.getWindow().setVisible(true);
		threshold_imp.getWindow().toFront();
		
		thresholds_gd.showDialog();
		
		if(thresholds_gd.wasCanceled())
		{
			threshold_imp.close();
			return false;
		}
		
		// get user specified threshold
		high = thresholds_gd.getNextNumber();
		Prefs.set("stegers_source.high",high);
		low = thresholds_gd.getNextNumber();
		Prefs.set("stegers_source.low",low);
		//*/
		threshold_imp.close();
		
		return true;

	}
	
	//////////contents of width.c
	
	/** Extract the line width by using a facet model line detector on an image of
	   the absolute value of the gradient. */
	public void compute_line_width()
	{
	  double[] dx = gradx;
	  double[] dy = grady;
	  double[]   grad;
	  int    i, j, k;
	  int    r, c;
	  int l;
	  long    x, y, dir;
	  ArrayList<offset> line;
	  int    num_line;//max_line, 
	  double  length;
	  contour tmp_cont;
	  long    num_points;// max_num_points;
	  float[]   width_r, width_l;
	  float[]  grad_r, grad_l;
	  float[]   pos_x, pos_y;// correct, asymm, contrast;
	  double  d, dr, dc, drr, drc, dcc;
	  double  i1, i2, i3, i4, i5, i6, i7, i8, i9;
	  double  t1, t2, t3, t4, t5, t6;
	  double[][]  eigvalvect;	  
	  double  a, b, t;
	  //long    num;
	  double  nx, ny;
	  double  n1, n2;
	  double  p1, p2;
	  double  val;
	  double  px, py;
	  int    num_contours;
	  


	  /*  max_num_points = 0;
	  for (i=0; i<num_contours; i++) {
	    num_points = contours[i]->num;
	    if (num_points > max_num_points)
	      max_num_points = num_points;
	  }

	  width_l = xcalloc(max_num_points,sizeof(*width_l));
	  width_r = xcalloc(max_num_points,sizeof(*width_r));
	  grad_l = xcalloc(max_num_points,sizeof(*grad_l));
	  grad_r = xcalloc(max_num_points,sizeof(*grad_r));
	  pos_x = xcalloc(max_num_points,sizeof(*pos_x));
	  pos_y = xcalloc(max_num_points,sizeof(*pos_y));
	  correct = xcalloc(max_num_points,sizeof(*correct));
	  contrast = xcalloc(max_num_points,sizeof(*contrast));
	  asymm = xcalloc(max_num_points,sizeof(*asymm));

	  grad = xcalloc(width*height,sizeof(*grad));
	  memset(grad,0,width*height*sizeof(*grad));*/
	  grad = new double[(int) (width*height)];

	  length = MAX_LINE_WIDTH;
	/*  max_line = ceil(length*3);
	  line = xcalloc(max_line,sizeof(*line));*/

	  /* Compute the gradient image. */
	  for (r=0; r<height; r++) {
	    for (c=0; c<width; c++) {
	      l = (int) LINCOOR(r,c,width);
	      grad[l] = Math.sqrt(dx[l]*dx[l]+dy[l]*dy[l]);
	    }
	  }

	  num_contours = this.cont.size();
	  for (i=0; i<num_contours; i++) {
	    tmp_cont = this.cont.get((int) i);
	    num_points = tmp_cont.num;
	    
	    pos_x = new float[(int) num_points];
	    pos_y = new float[(int) num_points];
	    grad_r = new float[(int) num_points];
	    grad_l = new float[(int) num_points];
	    width_r = new float[(int) num_points];
	    width_l = new float[(int) num_points];

	    
	    
	    for (j=0; j<num_points; j++) {	    		    	
	      px = tmp_cont.row.get(j);
	      py = tmp_cont.col.get(j);
	      pos_x[j] = (float) px;
	      pos_y[j] = (float) py;
	      r = (int) Math.floor(px+0.5);
	      c = (int) Math.floor(py+0.5);
	      nx = Math.cos(tmp_cont.angle.get(j));
	      ny = Math.sin(tmp_cont.angle.get(j));
	      // Compute the search line. 
	      line = bresenham(nx,ny,0.0,0.0,length);
	      width_r[j] = width_l[j] = 0;
	      // Look on both sides of the line. 
	      for (dir=-1; dir<=1; dir+=2) {
	    	  num_line = line.size();
	        for (k=0; k<num_line; k++) {
	          x = BR(r+dir*line.get(k).x);
	          y = BC(c+dir*line.get(k).y);
	          i1 = grad[(int) LINCOOR(BR(x-1),BC(y-1),width)];
	          i2 = grad[(int) LINCOOR(BR(x-1),y,width)];
	          i3 = grad[(int) LINCOOR(BR(x-1),BC(y+1),width)];
	          i4 = grad[(int) LINCOOR(x,BC(y-1),width)];
	          i5 = grad[(int) LINCOOR(x,y,width)];
	          i6 = grad[(int) LINCOOR(x,BC(y+1),width)];
	          i7 = grad[(int) LINCOOR(BR(x+1),BC(y-1),width)];
	          i8 = grad[(int) LINCOOR(BR(x+1),y,width)];
	          i9 = grad[(int) LINCOOR(BR(x+1),BC(y+1),width)];
	          t1 = i1+i2+i3;
	          t2 = i4+i5+i6;
	          t3 = i7+i8+i9;
	          t4 = i1+i4+i7;
	          t5 = i2+i5+i8;
	          t6 = i3+i6+i9;
	          dr = (t3-t1)/6;
	          dc = (t6-t4)/6;
	          drr = (t1-2*t2+t3)/6;
	          dcc = (t4-2*t5+t6)/6;
	          drc = (i1-i3-i7+i9)/4;
	          //eigvalvect = compute_eigenvals(2*drr,drc,2*dcc,eigval,eigvec);
	          eigvalvect = position.compute_eigenvals(2*drr,drc,2*dcc);
	          val = -eigvalvect[0][0];
	          if (val > 0.0) {
	            n1 = eigvalvect[1][0];
	            n2 = eigvalvect[1][1];
	            a = 2.0*(drr*n1*n1+drc*n1*n2+dcc*n2*n2);
	            b = dr*n1+dc*n2;
	            //solve_linear(a,b,&t,&num);
	            t= (-1)*b/a;
	            //if (num != 0) {
	            if (!Double.isNaN(t)) {
	              p1 = t*n1;
	              p2 = t*n2;
	              if (Math.abs(p1) <= 0.5 && Math.abs(p2) <= 0.5) {
	                /* Project the maximum point position perpendicularly onto the
	                   search line. */
	                a = 1;
	                b = nx*(px-(r+dir*line.get(k).x+p1))+ny*(py-(c+dir*line.get(k).y+p2));
	                //solve_linear(a,b,&t,&num);
	                t= (-1)*b/a;
	                d = (-i1+2*i2-i3+2*i4+5*i5+2*i6-i7+2*i8-i9)/9;
	                if (dir == 1) {
	                  grad_r[j] = (float) (d+p1*dr+p2*dc+p1*p1*drr+p1*p2*drc+p2*p2*dcc);
	                  width_r[j] = (float) Math.abs(t);
	                } else {
	                  grad_l[j] = (float) (d+p1*dr+p2*dc+p1*p1*drr+p1*p2*drc+p2*p2*dcc);
	                  width_l[j] = (float) Math.abs(t);
	                }
	                break;
	              }
	            }
	          }
	        }
	      }
	    }

	    fix_locations(width_l,width_r,grad_l,grad_r,pos_x,pos_y,tmp_cont);
	  }
/*
	  free(line);
	  free(grad);
	  free(asymm);
	  free(contrast);
	  free(correct);
	  free(pos_y);
	  free(pos_x);
	  free(grad_r);
	  free(grad_l);
	  free(width_r);
	  free(width_l);*/
	}

	
	
	/** Correct the extracted line positions and widths.  The algorithm first closes
	   gaps in the extracted data width_l, width_r, grad_l, and grad_r to provide
	   meaningful input over the whole line.  Then the correction is calculated.
	   After this, gaps that have been introduced by the width correction are again
	   closed.  Finally, the position correction is applied if correct_pos is set.
	   The results are returned in width_l, width_r, and cont. */
	public void fix_locations(float[] width_l,float[] width_r, float[] grad_l,float[] grad_r, float[] pos_x, float[] pos_y, contour cont)
	{
	  int    i;
	  long    num_points;
	  double  px, py;
	  double  nx, ny;
	  double  w_est, r_est, w_real, h_real, corr;//, w_strong, w_weak;
	  double  correct, asymmetry, response, widthx, contrast;
	  boolean    weak_is_r;
	  boolean    correct_start, correct_end;
	  
	  float[] asymm;
	  float[] contr;
	  float[][] tempd;
	  float[] correction;
	  correctionx corrx = new correctionx(); 

	  tempd = fill_gaps(width_l,grad_l,null,cont);
	  width_l = tempd[0];
	  grad_l = tempd[1];
	  tempd =fill_gaps(width_r,grad_r,null,cont);
	  width_r = tempd[0];
	  grad_r = tempd[1];
	  
	  num_points = cont.num;
	  asymm = new float[(int) num_points];
	  contr = new float[(int) num_points];
	  correction = new float[(int) num_points];

	  // Calculate true line width, asymmetry, and position correction. 
	  if(correct_pos) 
	  {
	    // Do not correct the position of a junction point if its width is found
	    // by interpolation, i.e., if the position could be corrected differently
	    // for each junction point, thereby destroying the junction. */
	    correct_start = ((cont.cont_class == contour_class.cont_no_junc ||
	                      cont.cont_class == contour_class.cont_end_junc ||
	                      cont.cont_class == contour_class.cont_closed) &&
	                     (width_r[0] > 0 && width_l[0] > 0));
	    correct_end = ((cont.cont_class == contour_class.cont_no_junc ||
	                    cont.cont_class == contour_class.cont_start_junc ||
	                    cont.cont_class == contour_class.cont_closed) &&
	                   (width_r[(int) (num_points-1)] > 0 && width_l[(int) (num_points-1)] > 0));
	    // Calculate the true width and assymetry, and its corresponding
	    //   correction for each line point. 
	    for (i=0; i<num_points; i++) 
	    {
	      if (width_r[i] > 0 && width_l[i] > 0) 
	      {
	        w_est = (width_r[i]+width_l[i])*LINE_WIDTH_COMPENSATION;
	        if (grad_r[i] <= grad_l[i]) 
	        {
	          r_est = grad_r[i]/grad_l[i];
	          weak_is_r = true;
	        } 
	        else 
	        {
	          r_est = grad_l[i]/grad_r[i];
	          weak_is_r = false;
	        }
	        //corrx = correct.line_corrections(new Double(SIGMA), new Double(w_est), new Double(r_est));
	        corrx = correctx.line_corrections(SIGMA, w_est, r_est);
	        
	        w_real = corrx.w;
	        h_real = corrx.h;
	        corr = corrx.correction;
	        //w_strong = corrx.w_strong;
	        //w_weak = corrx.w_weak;
	        w_real /= LINE_WIDTH_COMPENSATION;
	        corr /= LINE_WIDTH_COMPENSATION;
	        width_r[i] = (float) w_real;
	        width_l[i] = (float) w_real;
	        if (weak_is_r) 
	        {
	          asymm[i] = (float) h_real;
	          correction[i] = (float) (-corr);
	        } 
	        else 
	        {
	          asymm[i] = (float) (-h_real);
	          correction[i] = (float) corr;
	        }
	      }
	    }

	    tempd = fill_gaps(width_l,correction,asymm,cont);
	    width_l = tempd[0];
	    correction = tempd[1];
	    asymm = tempd[2];
	    
	    for (i=0; i<num_points; i++)
	      width_r[i] = width_l[i];

	    /* Adapt the correction for junction points if necessary. */
	    if (!correct_start)
	      correction[0] = 0;
	    if (!correct_end)
	      correction[(int) (num_points-1)] = 0;

	    for (i=0; i<num_points; i++) 
	    {
	      px = pos_x[i];
	      py = pos_y[i];
	      nx = Math.cos(cont.angle.get(i));
	      ny = Math.sin(cont.angle.get(i));
	      px = px+correction[i]*nx;
	      py = py+correction[i]*ny;
	      pos_x[i] = (float) px;
	      pos_y[i] = (float) py;
	    }
	  }


	  /* Update the position of a line and add the extracted width. */
	  //cont->width_l = xcalloc(num_points,sizeof(float));
	  //cont->width_r = xcalloc(num_points,sizeof(float));
	  
	  cont.width_l = new ArrayList<Float>();
	  cont.width_r = new ArrayList<Float>();
	  for (i=0; i<num_points; i++) 
	  {
	    cont.width_l.add(width_l[i]);
	    cont.width_r.add(width_r[i]);
	    cont.row.set(i,(float) pos_x[i]);
	    cont.col.set(i,(float) pos_y[i]);
	  }

	  /* Now calculate the true contrast. */
	  if (correct_pos) 
	  {
	    //cont->asymmetry = xcalloc(num_points,sizeof(float));
	    //cont->contrast = xcalloc(num_points,sizeof(float));
	    for (i=0; i<num_points; i++) 
	    {
	      response = cont.response.get(i);
	      asymmetry = Math.abs(asymm[i]);
	      correct = Math.abs(correction[i]);
	      widthx = cont.width_l.get(i);
	      if (widthx < MIN_LINE_WIDTH)
	        contrast = 0;
	      else
	        contrast = 
	          (response/Math.abs(convol.phi2(correct+widthx,SIGMA)+
	                         (asymmetry-1)*convol.phi2(correct-widthx,SIGMA)));
	      //seems like it always above threshold in my experience
	      //if (contrast > MAX_CONTRAST)
	      //  contrast = 0;
	      contr[i] = (float) contrast;
	    }
	  
	   tempd= fill_gaps(contr,null,null,cont);
	   contr=tempd[0];
	    
	    cont.asymmetry = new ArrayList<Float>();
  	    cont.contrast = new ArrayList<Float>();
	    
	    for (i=0; i<num_points; i++) 
	    {
	      cont.asymmetry.add(asymm[i]);
	      if (MODE_LIGHT)
	        cont.contrast.add(contr[i]);
	      else
	        cont.contrast.add(-contr[i]);
	    }
	  }
	
	
	}

	
	/** Fill gaps in the arrays master, slave1, and slave2, i.e., points where
	   master=0, by interpolation (interior points) or extrapolation (end points).
	   The array master will usually be the width of the line, while slave1 and
	   slave2 will be values that depend on master[i] being 0, e.g., the gradient
	   at each line point.  The arrays slave1 and slave2 can be NULL. 
	 * @return **/
	public static float[][] fill_gaps(float[] master, float[] slave1, float[] slave2, contour t_cont)
	{
	  int    i, j, k, s, e;
	  long    num_points;
	  double  m_s, m_e, s1_s, s1_e, s2_s, s2_e, d_r, d_c, arc_len, len;
	  
	  num_points = t_cont.num;
	  float [][] result_m_s = new float[3][master.length];

	  //for(i=0;i<master.length;i++)
	  //{
		  result_m_s[0]=master;
		  result_m_s[1]=slave1;
		  result_m_s[2]=slave2;
	  //}
	  
	  for (i=0; i<num_points; i++) {
	    if (result_m_s[0][i] == 0.0) {
	      for (j=i+1; j<num_points; j++) {
	        if (result_m_s[0][j] > 0)
	          break;
	      }
	      m_s = 0;
	      m_e = 0;
	      s1_s = 0;
	      s1_e = 0;
	      s2_s = 0;
	      s2_e = 0;
	      if (i > 0 && j < num_points-1) {
	        s = i;
	        e = j-1;
	        m_s = result_m_s[0][s-1];
	        m_e = result_m_s[0][e+1];
	        if (slave1 != null) {
	          s1_s = result_m_s[1][s-1];
	          s1_e = result_m_s[1][e+1];
	        }
	        if (slave2 != null) {
	          s2_s = result_m_s[2][s-1];
	          s2_e = result_m_s[2][e+1];
	        }
	      } else if (i > 0) {
	        s = i;
	        e = (int) (num_points-2);
	        m_s = result_m_s[0][s-1];
	        m_e = result_m_s[0][s-1];
	        result_m_s[0][e+1] = (float) m_e;
	        if (slave1 != null) {
	          s1_s = result_m_s[1][s-1];
	          s1_e = result_m_s[1][s-1];
	          result_m_s[1][e+1] = (float) s1_e;
	        }
	        if (slave2 != null) {
	          s2_s = result_m_s[2][s-1];
	          s2_e = result_m_s[2][s-1];
	          result_m_s[2][e+1] = (float) s2_e;
	        }
	      } else if (j < num_points-1) {
	        s = 1;
	        e = j-1;
	        m_s = result_m_s[0][e+1];
	        m_e = result_m_s[0][e+1];
	        result_m_s[0][s-1] = (float) m_s;
	        if (slave1 != null) {
	          s1_s = result_m_s[1][e+1];
	          s1_e = result_m_s[1][e+1];
	          result_m_s[1][s-1] = (float) s1_s;
	        }
	        if (slave2 != null) {
	          s2_s = result_m_s[2][e+1];
	          s2_e = result_m_s[2][e+1];
	          result_m_s[2][s-1] = (float) s2_s;
	        }
	      } else {
	        s = 1;
	        e = (int) (num_points-2);
	        m_s = result_m_s[0][s-1];
	        m_e = result_m_s[0][e+1];
	        if (slave1 != null) {
	          s1_s = result_m_s[1][s-1];
	          s1_e = result_m_s[1][e+1];
	        }
	        if (slave2 != null) {
	          s2_s = result_m_s[2][s-1];
	          s2_e = result_m_s[2][e+1];
	        }
	      }
	      arc_len = 0;
	      for (k=s; k<=e+1; k++) {
	        d_r = t_cont.row.get(k)-t_cont.row.get(k-1);
	        d_c = t_cont.col.get(k)-t_cont.col.get(k-1);
	        arc_len += Math.sqrt(d_r*d_r+d_c*d_c);
	      }
	      len = 0;
	      for (k=s; k<=e; k++) {
	        d_r = t_cont.row.get(k)-t_cont.row.get(k-1);
	        d_c = t_cont.col.get(k)-t_cont.col.get(k-1);
	        len += Math.sqrt(d_r*d_r+d_c*d_c);
	        result_m_s[0][k] = (float) ((arc_len-len)/arc_len*m_s+len/arc_len*m_e);
	        if (slave1 != null)
	        	result_m_s[1][k] = (float) ((arc_len-len)/arc_len*s1_s+len/arc_len*s1_e);
	        if (slave2 != null)
	        	result_m_s[2][k] = (float) ((arc_len-len)/arc_len*s2_s+len/arc_len*s2_e);
	      }
	      i = j;
	    }
	  }
	  return result_m_s;
	}

	
	public void prepareOrigMSE()
	{
		//ImageStatistics imgstat;
		float[] pxOrig;
		float[] pxBG;
		int i,j,ind;
		float [] ModeSD;
		FloatProcessor flBG;
		origMinusMean=(FloatProcessor) ip.duplicate().convertToFloat();	
		
		/**/
		//subtract background
		flBG=(FloatProcessor) origMinusMean.duplicate();
		flBG.blurGaussian(SIGMA*5.0);
		pxOrig=(float[])origMinusMean.getPixels();
		pxBG=(float[])flBG.getPixels();
		for(j=0; j<height; j++) 
		{
			for(i=0; i<width; i++) 
			{
				ind=i+j*(int)width;
				origMinusMean.setf(ind, pxOrig[ind]-pxBG[ind]);
			}
		}
		
		/**/
		ModeSD=getThreshold(origMinusMean);
		imgSD=ModeSD[1];
		imgMAX=ModeSD[2]-ModeSD[0];
		//imgstat=ImageStatistics.getStatistics(origMinusMean, Measurements.MIN_MAX, null);
		//origMinusMean.subtract(imgstat.mean);
		origMinusMean.subtract(ModeSD[0]);
		//ImagePlus impx;
		rt=new ResultsTable();
		//Analyzer anLoc;
		impx = new ImagePlus("measure", origMinusMean);
		anLoc=new Analyzer(impx,Measurements.MEAN,rt);
		//anLoc=new Analyzer(impx,Measurements.MODE,rt);
		//impx.show();
		new ImagePlus("corrected", origMinusMean.duplicate()).show();
	}
	@Override
	/** function returning MSE between traced and convoluted image 
	 * and original image**/
	public double value(double[] arg0) {
		// TODO Auto-generated method stub
		PolygonRoi polyline_p;
		contour tmp_cont;
		long num_pnt;
		FloatProcessor conv= null; 
		GaussianBlur gBlur = new GaussianBlur();
		int i,j,ind;
		
		//ImageStatistics imgstat;
		float[] pxTrace;
		float[] pxOrig;
		double dMSE=0;
		//double dImW=0;
		double lineval;
		

		double [] limits = new double [2];
		//this one should be calculated already
		//compute_line_points();
		
		//#of evaluations count
		nEvalCount=nEvalCount+1;
		
		//check for the boundaries
		if(arg0[0]>arg0[1])
		{
			limits[0]=arg0[1];
			limits[1]=arg0[0];
		}
		else
		{
			limits[0]=arg0[0];
			limits[1]=arg0[1];
		}
		if(limits[0]<threshold_tmp_ip_min)
			limits[0]=threshold_tmp_ip_min;
		if(limits[1]>threshold_tmp_ip_max)
			limits[1]=threshold_tmp_ip_max;
		
		// get valid points lines map
		getismax(limits); 
		compute_contours();
		//new imageprocessor
		conv = new FloatProcessor((int)width, (int)height);
		conv.setColor(imgMAX*0.5);
		conv.setLineWidth(1);
		
		//if(cont.size()>0)
		//{
			//draw contours
			for (i=0; i<cont.size(); i++) 
			{
				  tmp_cont = cont.get(i);  
				  num_pnt = tmp_cont.num;
				  float[] pxs = new float[(int) num_pnt];
				  float[] pxy = new float[(int) num_pnt];
				  for(j=0;j<num_pnt;j++)
				  {
					  pxs[j]=tmp_cont.row.get(j);
					  pxy[j]=tmp_cont.col.get(j);
				  }
					  			  
				  polyline_p = new PolygonRoi(pxy, pxs, Roi.POLYLINE);
				  
				  //polyline_p.setStrokeColor(Color.WHITE);
				  polyline_p.setStrokeWidth(1.0);
				  
				  //origMinusMean.setMask(polyline_p.getMask());
				  //origMinusMean.setRoi(polyline_p);
				  impx.setRoi(polyline_p, false);
				  anLoc.measure();
				  lineval=rt.getValue("Mean", rt.getCounter()-1);

				  //if(lineSD/Math.abs(lineval)<0.6)
				  //{
					  conv.setColor(lineval*SIGMA/0.342);
					  //conv.setColor(imgMAX*0.5);
				  //}
				  //else
				  //{
					//  conv.setColor((-1.0)*Math.abs(lineval*SIGMA/0.342));
				  //}
				 // IJ.log(Double.toString(lineval));
				  conv.draw(polyline_p);
			}
			gBlur.blurFloat(conv, SIGMA, SIGMA, 0.0002);

			new ImagePlus("convoluted_"+Integer.toString(nEvalCount), conv.duplicate()).show();

		//}
		//else
		//{
			//penalize absense of lines by multiple SD of original image
		//	conv.add(imgSD*4.0);
		//}
		pxTrace = (float[])conv.getPixels();
		pxOrig = (float[])origMinusMean.getPixels();
		
		//calculate MSE
		for(j=0; j<height; j++) 
		{
			for(i=0; i<width; i++) 
			{
				ind = i+j*(int)width;
				//double dev=Math.pow(pxTrace[ind]-pxOrig[ind],2);
				//dImW+=pxTrace[ind]*pxTrace[ind];
				double dev=Math.abs(pxTrace[ind]-pxOrig[ind]);
				//dImW+=Math.abs(pxTrace[ind]);

				if(!Double.isNaN(dev))
					dMSE+=dev;
				else
					dev=0;
				//dMSE+=Math.abs(pxTrace[i+j*(int)width]-pxOrig[i+j*(int)width]);
				
				
			}
		}
		dMSE=dMSE/((double)(width*height));
		//dMSE=dMSE/dImW;
		IJ.log("Eval: "+Integer.toString(nEvalCount)+" MSE:"+Double.toString(dMSE) + " lb:"+Double.toString(limits[0])+ " ub:"+Double.toString(limits[1]) +" curveN:"+Integer.toString(cont.size()));
		ptable.add(new Double[] { (double)(nEvalCount-5), dMSE,limits[0],limits[1], (double)(cont.size())});
		//ptable.addValue("Eval", nEvalCount-5);
		//ptable.addValue("MSE", dMSE);
		//ptable.addValue("LB", dMSE);
		//ptable.addValue("UB", dMSE);
		//ptable.addValue("CurveN", dMSE);
		nCurveN=cont.size();
		return dMSE;
	}


	/** 
	 *  returns value of mode intensity[0] and SD[1] in float processor based on 
	 *  fitting of image histogram to Gaussian function
	 *  In addition, returns max [2]
	**/
	float [] getThreshold(ImageProcessor thImage)
	{
		ImageStatistics imgstat;

		double  [][] dNoiseFit;
		int nHistSize;
		int nMaxCount;
		int nDownCount, nUpCount;
		int i,nPeakPos,k; 
		double dRightWidth, dLeftWidth;
		double dWidth=0;
		double dMean, dSD;
		double dMAX;
		double [] dFitErrors;
		double dErrCoeff;
		LMA fitlma;
		float [] results;
		int [] nHistgr;
		int nBinSizeEst = 256;
		boolean bOptimal = false;
		int nPeakNew;
		
		
		
		nBinSizeEst =getBinOptimalNumber(thImage);
		
		//searching for the optimal for fitting intensity histogram's bin size
	
		thImage.setHistogramSize(nBinSizeEst);			
		imgstat = ImageStatistics.getStatistics(thImage, Measurements.MODE + Measurements.MEAN+Measurements.STD_DEV+Measurements.MIN_MAX, null);
		nHistSize = imgstat.histogram.length;											
		nPeakPos = imgstat.mode;
		nMaxCount = imgstat.maxCount;
		dMAX=imgstat.max;
		
		
		results = new float [] {(float) imgstat.dmode,(float) imgstat.stdDev, (float)dMAX};
		return results;
		
		/*
		nHistgr = new int [nHistSize];
		for(k=0;k<nHistSize;k++)
			nHistgr[k]=imgstat.histogram[k];
		
		//Plot histplot = new Plot("Histogram","intensity", "count", dHistogram[0], dHistogram[1]);
		//histplot.show();

		while (!bOptimal)
		{
			//estimating width of a peak
			//going to the left
			i = nPeakPos;
			while (i>0 && nHistgr[i]>0.5*nMaxCount)
			{
				i--;			
			}
			if(i<0)
				i=0;
			dLeftWidth = i;
			//going to the right
			i=nPeakPos;
			while (i<nHistSize && nHistgr[i]>0.5*nMaxCount)
			{
				i++;			
			}
			if(i==nHistSize)
				i=nHistSize-1;
			dRightWidth = i;
			//FWHM in bins
			dWidth = (dRightWidth-dLeftWidth);
			//histogram is too narrow for fitting, increase number of bins
			if(dWidth<12)
			{
				nBinSizeEst = nBinSizeEst + 100;
				//bin set is too dense
				if(nBinSizeEst> 1000)
				{
					//ok, seems there is one very high peak/bin, let's remove it					
					//nBinSizeEst = 256;
					thImage.setHistogramSize(nBinSizeEst);	
					imgstat = ImageStatistics.getStatistics(thImage, Measurements.MODE + Measurements.MEAN+Measurements.STD_DEV+Measurements.MIN_MAX, null);
					nHistSize = imgstat.histogram.length;											
					nPeakPos = imgstat.mode;
	
					
					nHistgr = new int [nHistSize];
					nPeakNew=0; nMaxCount=0;
					for(k=0;k<nPeakPos;k++)
					{
						nHistgr[k]=imgstat.histogram[k];
						if(nHistgr[k]>nMaxCount)
						{
							nMaxCount = nHistgr[k];
							nPeakNew = k;
						}
					}
					//no particles or flat image
					if (nPeakPos==0)
					{
						results = new float [] {(float) imgstat.mean,(float) imgstat.stdDev, (float)dMAX};
						return results;
					}
					nHistgr[nPeakPos]=(int) (0.5*(imgstat.histogram[nPeakPos-1]+imgstat.histogram[nPeakPos+1]));
					for(k=nPeakPos+1;k<nHistSize;k++)
					{
						nHistgr[k]=imgstat.histogram[k];
						if(nHistgr[k]>nMaxCount)
						{
							nMaxCount = nHistgr[k];
							nPeakNew = k;
						}
					}
					nPeakPos=nPeakNew;
					//nHistSize = nHistSize;
					
					//estimating width of a peak
					//going to the left
					i = nPeakPos;
					while (i>0 && nHistgr[i]>0.5*nMaxCount)
					{
						i--;			
					}
					if(i<0)
						i=0;
					dLeftWidth = i;
					//going to the right
					i=nPeakPos;
					while (i<nHistSize && nHistgr[i]>0.5*nMaxCount)
					{
						i++;			
					}
					if(i==nHistSize)
						i=nHistSize-1;
					dRightWidth = i;
					//FWHM in bins
					dWidth = (dRightWidth-dLeftWidth);
					bOptimal = true;
				}
				else
				{
					//recalculate parameters
					thImage.setHistogramSize(nBinSizeEst);			
					imgstat = ImageStatistics.getStatistics(thImage, Measurements.MODE + Measurements.MEAN+Measurements.STD_DEV+Measurements.MIN_MAX, null);
					nHistSize = imgstat.histogram.length;											
					nPeakPos = imgstat.mode;
					nMaxCount = imgstat.maxCount;
					nHistgr = new int[nHistSize];
					for(k=0;k<nHistSize;k++)
						nHistgr[k]=imgstat.histogram[k];					
				}
			}
			//histogram is ok, proceed to fitting
			else
			{
				bOptimal = true;				
			}
		}

					
		dMean = imgstat.min + nPeakPos*imgstat.binSize;
		dSD = dWidth*imgstat.binSize/2.35;
		//fitting range +/- 3*SD
		dLeftWidth = nPeakPos - 3*dWidth/2.35;
		if(dLeftWidth<0)
			dLeftWidth=0;
		dRightWidth = nPeakPos + 3*dWidth/2.35;
		if(dRightWidth>nHistSize)
			dRightWidth=nHistSize;
		nUpCount = (int)dRightWidth;
		nDownCount = (int)dLeftWidth;
		//preparing histogram range for fitting
		dNoiseFit = new double [2][nUpCount-nDownCount+1];
		for(i=nDownCount;i<=nUpCount;i++)
		{
			dNoiseFit[0][i-nDownCount] = imgstat.min + i*imgstat.binSize;
			dNoiseFit[1][i-nDownCount] = (double)nHistgr[i];
		}
		
		fitlma = new LMA(new OneDGaussian(), new double[] {(double)nMaxCount, dMean, dSD}, dNoiseFit);
		fitlma.fit();
		dMean = fitlma.parameters[1];
		dSD = fitlma.parameters[2];
		
		dFitErrors = fitlma.getStandardErrorsOfParameters();
		// scaling coefficient for parameters errors estimation 
		// (Standard deviation of residuals)
		dErrCoeff = Math.sqrt(fitlma.chi2/(nUpCount-nDownCount+1-3));
		for (i=0;i<3;i++)
			dFitErrors[i] *= dErrCoeff;
		for (i=0;i<3;i++)
			dFitErrors[i] *= 100/fitlma.parameters[i]; 
		
	
		
			results = new float [] {(float) dMean,(float) dSD, (float)dMAX};
			return results;
	
		*/
		
	}
	
	/** function returns optimal bin number for the image histogram
	 *  according to the Freedman-Diaconis rule (check wiki) **/
	int getBinOptimalNumber(ImageProcessor ip)
	{
		//int nBinSize;
		int width, height;
		int pixelCount;

		
		width=ip.getWidth();
		height=ip.getHeight();
		pixelCount=width*height;
		

		float[] pixels2 = new float[pixelCount];
		//float[] pixels2;
		System.arraycopy((float[])ip.getPixels(),0,pixels2,0,pixelCount);
	
		Arrays.sort(pixels2);
		//int middle = pixels2.length/2;
		int qi25 = Math.round(pixelCount*0.25f);
		int qi75 = Math.round(pixelCount*0.75f);
	
		float IQR = pixels2[qi75]-pixels2[qi25];
		double h= 2*IQR*Math.pow((double)pixelCount, -1.0/3.0);
		
		return (int)Math.round((pixels2[pixelCount-1]-pixels2[0])/h);
			
	}
	
	/** given an image function returns intensity value
	 * corresponding to specific percentile (saturation level) **/
	float getSaturationLevel(ImageProcessor ip, double dSatPercentage)
	{
		int width, height;
		int pixelCount;

		
		width=ip.getWidth();
		height=ip.getHeight();
		pixelCount=width*height;
		

		float[] pixels2 = new float[pixelCount];
		//float[] pixels2;
		System.arraycopy((float[])ip.getPixels(),0,pixels2,0,pixelCount);
	
		Arrays.sort(pixels2);
		//int middle = pixels2.length/2;
		long satind = Math.round(pixelCount*(1.0-(dSatPercentage/100)));
		//int qi75 = Math.round(pixelCount*0.75f);
	
		return pixels2[(int) satind];
		//float IQR = pixels2[qi75]-pixels2[qi25];
		//double h= 2*IQR*Math.pow((double)pixelCount, -1.0/3.0);
		
		//return (int)Math.round((pixels2[pixelCount-1]-pixels2[0])/h);
	
	}
	
	public void calcMap()
	{
		double lb,ub;
		double nMin=0.1;
		double nMax=15;
		double nStep=0.1;
		int i=-1;
		int j=-1;
		int nDim=(int)Math.ceil((nMax-nMin)/nStep)+1;
		ipCurveN=new FloatProcessor(nDim,nDim);
		ipMSE=new FloatProcessor(nDim,nDim);
		
		for( lb=nMin;lb<nMax; lb+=nStep)
		{
			i++; j=i-1;
			for(ub=lb;ub<nMax; ub+=nStep)
			{
				j++;
				ipMSE.setf(i, j, (float)value(new double[] { lb, ub}));
				ipCurveN.setf(i, j, (float)nCurveN);
				ipMSE.setf(j, i, (float)value(new double[] { lb, ub}));
				ipCurveN.setf(j, i, (float)nCurveN);
			}
			IJ.log(Double.toString(lb));
		}
	}


}

