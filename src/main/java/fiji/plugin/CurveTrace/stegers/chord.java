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

/** A chord in a run-length encoded region */
public class chord {
	
	 /** row coordinate of the chord */
	 public short r;   
	 /** column coordinate of the start of the chord */
	 public short cb;  
	 /** column coordinate of the end of the chord */
	 public short ce;  
	  
		
	  public chord()
	  {
				r=0;
				cb=0;
				ce=0;
	  }
	  
	  public chord (short nr, short ncb, short nce)
	  {
			r = nr;
			cb = ncb;
			ce = nce;
	  }

}
