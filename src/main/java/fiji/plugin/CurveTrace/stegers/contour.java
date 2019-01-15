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

import java.util.ArrayList;


/** This type holds one extracted line.  The field num contains the number of
points in the line.  The coordinates of the line points are given in the
arrays row and col.  The array angle contains the direction of the normal
to each line point, as measured from the row-axis.  Some people like to
call the col-axis the x-axis and the row-axis the y-axis, and measure the
angle from the x-axis.  To convert the angle into this convention, subtract
PI/2 from the angle and normalize it to be in the interval [0,2*PI).  The
array response contains the response of the operator, i.e., the second
directional derivative in the direction of angle, at each line point.  The
arrays width_l and width_r contain the width information for each line point
if the algorithm was requested to extract it; otherwise they are NULL.  If
the line position and width correction was applied the contents of width_l
and width_r will be identical.  The arrays asymmetry and contrast contain
the true asymmetry and contrast of each line point if the algorithm was
instructed to apply the width and position correction.  Otherwise, they are
set to NULL.  If the asymmetry, i.e., the weaker gradient, is on the right
side of the line, the asymmetry is set to a positive value, while if it is
on the left side it is set to a negative value. */

public class contour {
	/** number of points */
	public long  num
	/** row coordinates of the line points (Y coordinate in ImageJ) */;                
	public ArrayList<Float> row;
	/** column coordinates of the line points (X coordinate in ImageJ)  */
	public ArrayList<Float> col;
	/** angle of normal (measured from the row (Y) axis) */
	public ArrayList<Float> angle;
	/** response of line point (second derivative) */
	public ArrayList<Float> response;
	/** width to the left of the line */
	public ArrayList<Float> width_l;
	/** width to the right of the line */
	public ArrayList<Float> width_r;
	/** asymmetry of the line point */
	public ArrayList<Float> asymmetry;
	/** contrast of the line point */
	public ArrayList<Float> contrast;
	/** contour class (e.g., closed, no_junc) */
	public contour_class cont_class; 
	
	
	//default constructor
	public contour()
	{
		num = 0;
		cont_class = contour_class.cont_no_junc;
		row = new ArrayList<Float>();
		col = new ArrayList<Float>();
		angle = new ArrayList<Float>();
		response = new ArrayList<Float>();
		width_l = new ArrayList<Float>();
		width_r = new ArrayList<Float>();
		asymmetry = new ArrayList<Float>();
		contrast = new ArrayList<Float>();		
			
	}
	//constructor
	public contour(long nnum, ArrayList<Float> nrow, ArrayList<Float> ncol, ArrayList<Float> nangle, ArrayList<Float> nresponse, contour_class ncont_class)
	{
		num = nnum;
		row = nrow;
		col = ncol;
		angle = nangle;
		response = nresponse;
		width_l = new ArrayList<Float>();
		width_r = new ArrayList<Float>();
		asymmetry = new ArrayList<Float>();
		contrast = new ArrayList<Float>();
		cont_class = ncont_class;
	}
		
}
