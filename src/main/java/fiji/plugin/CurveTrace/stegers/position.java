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

import ij.process.ImageProcessor;

public class position {
	
	
	public double sigma;
	public ImageProcessor ip;
	public double[] image;
	public long width, height;
	
	public position(ImageProcessor cip, double csigma)
	{
		int l;
		ip=cip;
		sigma = csigma;
		width = ip.getWidth();
		height = ip.getHeight();
		
		image = new double[(int) (width*height)];
		for(int py = 0; py < height; ++py)		
			for(int px = 0; px < width; ++px)
			{
				//since x and y are switched
				l = (int) LINCOOR(py,px,width);
				image[l] = ip.getf(px,py);
			}
		
		
	}
	/** Compute the eigenvalues and eigenvectors of the Hessian matrix given by
	   dfdrr, dfdrc, and dfdcc, and sort them in descending order according to
	   their absolute values. */	
	public static double[][] compute_eigenvals(double dfdrr,double dfdrc,double dfdcc)	  
	{
		
		double[][] eigenvalvect = new double[3][2];
		//eigenvalvect[0][0] - first eigenvector
		//eigenvalvect[0][1] - second eigenvector
		//eigenvalvect[1][0] - first eigenvector x
		//eigenvalvect[1][1] - first eigenvector y
		//eigenvalvect[2][0] - second eigenvector x
		//eigenvalvect[2][1] - second eigenvector y
	  double theta, t, c, s, e1, e2, n1, n2; /* , phi; */

	  /* Compute the eigenvalues and eigenvectors of the Hessian matrix. */
	  if (dfdrc != 0.0) {
	    theta = 0.5*(dfdcc-dfdrr)/dfdrc;
	    t = 1.0/(Math.abs(theta)+Math.sqrt(theta*theta+1.0));
	    if (theta < 0.0) t = -t;
	    c = 1.0/Math.sqrt(t*t+1.0);
	    s = t*c;
	    e1 = dfdrr-t*dfdrc;
	    e2 = dfdcc+t*dfdrc;
	  } else {
	    c = 1.0;
	    s = 0.0;
	    e1 = dfdrr;
	    e2 = dfdcc;
	  }
	  n1 = c;
	  n2 = -s;

	  /* If the absolute value of an eigenvalue is larger than the other, put that
	     eigenvalue into first position.  If both are of equal absolute value, put
	     the negative one first. */
	  if (Math.abs(e1) > Math.abs(e2)) {
		  eigenvalvect[0][0] = e1;
		  eigenvalvect[0][1] = e2;
		  eigenvalvect[1][0] = n1;
		  eigenvalvect[1][1] = n2;
		  eigenvalvect[2][0] = -n2;
		  eigenvalvect[2][1] = n1;
	  } else if (Math.abs(e1) < Math.abs(e2)) {
		  eigenvalvect[0][0] = e2;
		  eigenvalvect[0][1] = e1;
		  eigenvalvect[1][0] = -n2;
		  eigenvalvect[1][1] = n1;
		  eigenvalvect[2][0] = n1;
		  eigenvalvect[2][1] = n2;
	  } else {
	    if (e1 < e2) {
	    	eigenvalvect[0][0] = e1;
	    	eigenvalvect[0][1] = e2;
	    	eigenvalvect[1][0] = n1;
	    	eigenvalvect[1][1] = n2;
	    	eigenvalvect[2][0] = -n2;
	    	eigenvalvect[2][1] = n1;
	    } else {
	    	eigenvalvect[0][0] = e2;
	    	eigenvalvect[0][1] = e1;
	    	eigenvalvect[1][0] = -n2;
	    	eigenvalvect[1][1] = n1;
	    	eigenvalvect[2][0] = n1;
	    	eigenvalvect[2][1] = n2;
	    }
	  }
	  return eigenvalvect;
	}

	/* Solve the linear equation a*x+b=0 and return the result in t and the number
	   of solutions in num. 
	public static solve_linear(double a,double b)
	  double a;
	  double b;
	  double *t;
	  long   *num;
	{
	  if (a == 0.0) {
	    *num = 0;
	    return;
	  } else {
	    *num = 1;
	    *t = -b/a;
	    return;
	  }
	}
	*/
	
	/** For each point in the image determine whether there is a local maximum of
	   the second directional derivative in the direction (nx[l],ny[l]) within the
	   pixels's boundaries.  If so, set ismax[l] to 2 if the eigenvalue ev[l] is
	   larger than high, to 1 if ev[l] is larger than low, and to 0 otherwise.
	   Furthermore, put the sub-pixel position of the maximum into (px[l],py[l]).
	   The parameter mode determines whether maxima (dark lines points) or minima
	   (bright line points) should be selected.  The partial derivatives of the
	   image are input as ku[]. */	
	public double[][] compute_line_points(boolean mode)
	{
		
		int r, c;
		// store results of eigendecomposition
		//	[0] = first eigenvalue
		//	[1] = eigenvector coordinate y
		//	[2] = eigenvector coordinate x
		//	[3] = derivative of image in y
		//	[4] = derivative of image in x
		//	[5] = "super-resolved" y		t_y, or dlpy
		//	[6] = "super-resolved" x		t_x, or dlpx		
		double[][] line_points = new double[7][(int) (width*height)];
		double[][] eigenvalvect; 
		
		//all derivatives
		double[][] k;
		int l;
		double a,b,n1,n2,t;
		
		convol conobj = new convol(width, height);
		 
		k = conobj.get_all_derivatives(image,sigma);
		
		line_points[3] = k[0];
		line_points[4] = k[1];
		
		/*compute eigenval and vectors */
		  for (r=0; r<height; r++) 
		  {
			    for (c=0; c<width; c++) 			    	
			    {
			    	 l = (int) LINCOOR(r,c,width);
			    	 eigenvalvect = compute_eigenvals(k[2][l],k[3][l],k[4][l]);
			    	 
			    	 if(mode)
			    		 line_points[0][l] =(-1)*eigenvalvect[0][0];
			    	 else
			    		 line_points[0][l] = eigenvalvect[0][0];
			    	 n1 =eigenvalvect[1][0];
			    	 n2 = eigenvalvect[1][1];
			    	 line_points[1][l] = n1;
			    	 line_points[2][l] = n2;
			    	 
			         a = k[2][l]*n1*n1+2.0*k[3][l]*n1*n2+k[4][l]*n2*n2;
			         b = k[0][l]*n1+k[1][l]*n2;
			         t = (-1)*b/a;
			 		 line_points[5][l] = r+t*n1;
					 line_points[6][l] = c+t*n2;
			         
			    }
		 }

		return line_points;
		
	}
	
	
	
	 /** Translate row and column coordinates of an image into an index into its
   one-dimensional array. */
	public long LINCOOR(long row, long col,long width){ 
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

}
