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

/** Data structure to store three doubles: x,y and t (distance along line) */
public class doublepoint {
	
	public double cx;
	public double cy;
	public double t;
	
	public doublepoint ()
	  {
			super();
	  }
	  
	  public doublepoint (double ccx, double ccy, double ct)
	  {
			cx = ccx;
			cy = ccy;
			t=ct;
	  }

	  public doublepoint (double ccx, double ccy)
	  {
			cx = ccx;
			cy = ccy;
			t=0;
	  }	  
}
