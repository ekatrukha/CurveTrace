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

/** This data structure is used to accumulate junction information.  It is
needed to split lines at junction points. */
	public class junction implements Comparable <junction>{
		  /** Index of line that is already processed */
		  public long  cont1; 
		  /** Index of line that runs into cont1 */
		  public long  cont2; 
		  /** Index of the junction point in cont1 */
		  public long  pos;   
		  /** y-(row-)coordinate of the junction point (corrected for ImageJ)*/
		  public float x;     
		  /** x-(col-)coordinate of the junction point (corrected for ImageJ)*/
		  public float y;     
		  
		  public junction()
		  {
				super();
		  }
		  
		  public junction(long ncont1, long ncont2, long npos, float nx, float ny)
		  {
				cont1=ncont1;
				cont2=ncont2;
				pos = npos;
				x = nx;
				y = ny;
		  }
		  
		  /** This function compares two junctions according to their first line indexes,
		   and, if needed, by the position of the junction within the line.  It is
		   called by qsort. */
		  public int compareTo(junction otherjunction) {
			if(this.cont1 == otherjunction.cont1)
				return (int)(this.pos-otherjunction.pos);
			else
				return (int)(this.cont1-otherjunction.cont1);
	
		  }
	} ;

