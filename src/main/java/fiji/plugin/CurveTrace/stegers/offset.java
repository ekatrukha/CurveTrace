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

/** Offsets to a specific location in the image.  An array of this type is
returned by the modified Bresenham algorithm in width.c.  It is also used
in link.c to hold an array of pixel locations to check for appropriate
neighbors. */
public class offset {
	public long x;
	public long y;
	
	  public offset()
	  {
			super();
	  }
	  
	  public offset(long nx, long ny)
	  {
			x = nx;
			y = ny;
	  }
	
}
