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

/** This data structure facilitates the quick search for the next possible
starting point of a line.  An array of crossrefs will be accumulated and
sorted according to its value.  x and y are the coordinates of a point in
the image.  When this point has been processed it will be marked as done. **/
public class crossref implements Comparable <crossref>{
	public short  x;
	public short  y;
	public double value;
	public boolean   done;

	  public crossref()
	  {
			super();
	  }
	  
	  public crossref(short nx, short ny, double nvalue, boolean ndone)
	  {
			x = nx;
			y = ny;
			value = nvalue;
			done = ndone;
	  }

	  /** This function compares two crossrefs according to their value.  It is called
	   by qsort. **/
	public int compareTo(crossref otherref) {
		// TODO Auto-generated method stub
		if (this.value>otherref.value)
			return -1;
		if (this.value<otherref.value)
			return 1;
		return 0;
	}
} ;
