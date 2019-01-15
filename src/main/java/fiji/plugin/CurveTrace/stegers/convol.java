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

public class convol {

	/** Derivative in row direction */
	public static final int DERIV_R  = 1;  
	/** Derivative in column direction */
	public static final int DERIV_C  = 2; 
	/** Second derivative in row direction */
	public static final int DERIV_RR = 3; 
	/** Second derivative in row and column direction */
	public static final int DERIV_RC = 4;  
	/** Second derivative in column direction */
	public static final int DERIV_CC = 5;  
	
	public static double SQRT_2_PI_INV =  0.398942280401432677939946059935;

/* Compute the integral of the Gaussian, i.e., the normal distribution. */

	public static double SQRTPI = 1.772453850905516027;

	public static double UPPERLIMIT=20.0; 

	public static double SQRT2 = 1.41421356237309504880;
	public static double P10 = 242.66795523053175;
	public static double P11 = 21.979261618294152;
	public static double P12 = 6.9963834886191355;
	public static double P13 = -0.035609843701815385;
	public static double Q10 = 215.05887586986120;
	public static double Q11 = 91.164905404514901;
	public static double Q12 = 15.082797630407787;
	public static double Q13 = 1.0;

	public static double P20 = 300.4592610201616005;
	public static double P21 = 451.9189537118729422;
	public static double P22 = 339.3208167343436870;
	public static double P23 = 152.9892850469404039;
	public static double P24 = 43.16222722205673530;
	public static double P25 = 7.211758250883093659;
	public static double P26 = 0.5641955174789739711;
	public static double P27 = -0.0000001368648573827167067;
	public static double Q20 = 300.4592609569832933;
	public static double Q21 = 790.9509253278980272;
	public static double Q22 = 931.3540948506096211;
	public static double Q23 = 638.9802644656311665;
	public static double Q24 = 277.5854447439876434;
	public static double Q25 = 77.00015293522947295;
	public static double Q26 = 12.78272731962942351;
	public static double Q27 = 1.0;

	public static double P30 = -0.00299610707703542174;
	public static double P31 = -0.0494730910623250734;
	public static double P32 = -0.226956593539686930;
	public static double P33 = -0.278661308609647788;
	public static double P34 = -0.0223192459734184686;
	public static double Q30 = 0.0106209230528467918;
	public static double Q31 = 0.191308926107829841;
	public static double Q32 = 1.05167510706793207;
	public static double Q33 = 1.98733201817135256;
	public static double Q34 = 1.0;

	/** Size for Gaussian mask */
	public static double MAX_SIZE_MASK_0=3.09023230616781; 
	/** Size for 1st derivative mask */
	public static double MAX_SIZE_MASK_1=3.46087178201605;    
	/** Size for 2nd derivative mask */
	public static double MAX_SIZE_MASK_2=3.82922419517181;    
	
	public static double MASK_SIZE(double MAX, double sigma) 
	{
		return Math.ceil(MAX*sigma);
	}

	
	
	public long width;
	public long height;
	
	public convol(long cwidth, long cheight)
	{
		width = cwidth;
		height = cheight; 
		
	}
	
	public static double normal(double x)
	{
	  int    sn;
	  double R1, R2, y, y2, y3, y4, y5, y6, y7;
	  double erf, erfc, z, z2, z3, z4;
	  double phi;

	  if (x < -UPPERLIMIT) return 0.0;
	  if (x > UPPERLIMIT) return 1.0;

	  y = x / SQRT2;
	  if (y < 0) {
	    y = -y;
	    sn = -1;
	  } else
	    sn = 1;

	  y2 = y * y;
	  y4 = y2 * y2;
	  y6 = y4 * y2;

	  if (y < 0.46875) {
	    R1 = P10 + P11 * y2 + P12 * y4 + P13 * y6;
	    R2 = Q10 + Q11 * y2 + Q12 * y4 + Q13 * y6;
	    erf = y * R1 / R2;
	    if (sn == 1)
	      phi = 0.5 + 0.5*erf;
	    else 
	      phi = 0.5 - 0.5*erf;
	  } else if (y < 4.0) {
	    y3 = y2 * y;
	    y5 = y4 * y;
	    y7 = y6 * y;
	    R1 = P20 + P21 * y + P22 * y2 + P23 * y3 + 
	      P24 * y4 + P25 * y5 + P26 * y6 + P27 * y7;
	    R2 = Q20 + Q21 * y + Q22 * y2 + Q23 * y3 + 
	      Q24 * y4 + Q25 * y5 + Q26 * y6 + Q27 * y7;
	    erfc = Math.exp(-y2) * R1 / R2;
	    if (sn == 1)
	      phi = 1.0 - 0.5*erfc;
	    else
	      phi = 0.5*erfc;
	  } else {
	    z = y4;
	    z2 = z * z;
	    z3 = z2 * z;
	    z4 = z2 * z2;
	    R1 = P30 + P31 * z + P32 * z2 + P33 * z3 + P34 * z4;
	    R2 = Q30 + Q31 * z + Q32 * z2 + Q33 * z3 + Q34 * z4;
	    erfc = (Math.exp(-y2)/y) * (1.0 / SQRTPI + R1 / (R2 * y2));
	    if (sn == 1)
	      phi = 1.0 - 0.5*erfc;
	    else 
	      phi = 0.5*erfc;
	  } 

	  return phi;
	}

	
	
	/** Integral of the Gaussian function */
	public static double phi0(double x,double sigma)
	 
	{
	  return normal(x/sigma);
	}
	
	/** The Gaussian function */
	public static double phi1(double x, double sigma) 
	{
	  double t;

	  t = x/sigma;
	  return SQRT_2_PI_INV/sigma*Math.exp(-0.5*t*t);
	  
	}
	
	/** First derivative of the Gaussian function */
	public static double phi2(double x, double sigma)	
	{
	  double t;

	  t = x/sigma;
	  return -x*SQRT_2_PI_INV/Math.pow(sigma,3.0)*Math.exp(-0.5*t*t);
	}
	
	/** Functions to compute the one-dimensional convolution masks of the 0th, 1st,
	   and 2nd derivative of the Gaussian kernel for a certain smoothing level
	   given by sigma.  The mask is allocated by the function and given as the
	   return value.  The caller must ensure that this memory is freed.  The
	   output is intended to be used as an array with range [-num:num].  Therefore,
	   the caller should add num to the return value.  Examples for the calling
	   sequence can be found in convolve_gauss.  Examples for the usage of the
	   masks are given in convolve_rows_gauss and convolve_cols_gauss. */

	/* Mask sizes in convol.c */


	/** Gaussian smoothing mask */
	public static double[] compute_gauss_mask_0(double sigma)
	{
	  int   i, n;
	  double limit;
	  double[] h;

	  limit = MASK_SIZE(MAX_SIZE_MASK_0,sigma); /* Error < 0.001 on each side */
	  n = (int) limit;
	  h = new double[2*n+1];
	  
	  for (i=-n+1;i<=n-1;i++)
	    h[i+n] = phi0(-i+0.5,sigma) - phi0(-i-0.5,sigma);
	  
	  h[0] = 1.0 - phi0(n-0.5,sigma);
	  h[2*n] = phi0(-n+0.5,sigma);

	  return h;
	}
	
	/** First derivative of Gaussian smoothing mask */
	public static double[] compute_gauss_mask_1(double sigma)

	{
	  int i, n;
	  double limit;
	  double []  h;

	  limit = MASK_SIZE(MAX_SIZE_MASK_1,sigma); /* Error < 0.001 on each side */
	  n = (int)limit;
	  h = new double[2*n+1];
	  	  
	  for (i=-n+1;i<=n-1;i++)
	    h[i+n] = phi1(-i+0.5,sigma) - phi1(-i-0.5,sigma);
	  
	  h[0] = -phi1(n-0.5,sigma);
	  h[2*n] = phi1(-n+0.5,sigma);

	  return h;
	}

	/** Second derivative of Gaussian smoothing mask */
	public static double[] compute_gauss_mask_2(double sigma)
	{
	  int i, n;
	  double limit;
	  double[] h;

	  limit = MASK_SIZE(MAX_SIZE_MASK_2,sigma); /* Error < 0.001 on each side */
	  n = (int)limit;
	  h = new double[2*n+1];
	  
	  for (i=-n+1;i<=n-1;i++)
	    h[i+n] = phi2(-i+0.5,sigma) - phi2(-i-0.5,sigma);
	  h[0] = -phi2(n-0.5,sigma);
	  h[2*n] = phi2(-n+0.5,sigma);
	  
	  return h;
	}
	

/** Convolve an image with the derivatives of a Gaussian smoothing kernel.
   Since all of the masks are separable, this is done in two steps in the
   function convolve_gauss.  Firstly, the rows of the image are convolved by
   an appropriate one-dimensional mask in convolve_rows_gauss, yielding an
   intermediate float-image h.  Then the columns of this image are convolved
   by another appropriate mask in convolve_cols_gauss to yield the final
   result k.  At the border of the image the gray values are mirrored. */
	
	
	/** Convolve the rows of an image with the derivatives of a Gaussian. **/
	public double[] convolve_rows_gauss(double[] image, double[] mask)
	{
	  long      j, r, c, l;
	  double    sum;
	  double [] h;
	  
	  
	  int n;
	  
	  h = new double[(int) (width*height)];
	  n=(int) ((mask.length-1)*0.5);
	  /* Inner region */
	  for (r=n; r<height-n; r++) {
	    for (c=0; c<width; c++) {
	      l = LINCOOR(r,c,width);
	      sum = 0.0;
	      for (j=-n;j<=n;j++)
	        sum += (double)(image[(int) (l+j*width)])*mask[(int) j+n];
	      h[(int) l] = sum;
	    }
	  }
	  /* Border regions */
	  for (r=0; r<n; r++) {
	    for (c=0; c<width; c++) {
	      l = LINCOOR(r,c,width);
	      sum = 0.0;
	      for (j=-n;j<=n;j++)
	        sum += (double)(image[(int) LINCOOR(BR(r+j),c,width)])*mask[(int) j+n];
	      h[(int) l] = sum;
	    }
	  }
	  for (r=height-n; r<height; r++) {
	    for (c=0; c<width; c++) {
	      l = LINCOOR(r,c,width);
	      sum = 0.0;
	      for (j=-n;j<=n;j++)
	        sum += (double)(image[(int) LINCOOR(BR(r+j),c,width)])*mask[(int) j+n];
	      h[(int) l] = sum;
	    }
	  }
	  return h;
	}
	

/** Convolve the columns of an image with the derivatives of a Gaussian. */
	public  double[] convolve_cols_gauss(double[] h,double[] mask) 
	{
	  long      j, r, c, l;
	  double    sum;
	  double[] k;
	  int n;
	  
	  n=(int) ((mask.length-1)*0.5);
	  
	  k=new double[(int) (width*height)];
	
	  // Inner region 
	  for (r=0; r<height; r++) {
	    for (c=n; c<width-n; c++) {
	      l = LINCOOR(r,c,width);
	      sum = 0.0;
	      for (j=-n;j<=n;j++)
	        sum += h[(int) (l+j)]*mask[(int) j+n];
	      k[(int) l] = (float)sum;
	    }
	  }
	  // Border regions 
	  for (r=0; r<height; r++) {
	    for (c=0; c<n; c++) {
	      l = LINCOOR(r,c,width);
	      sum = 0.0;
	      for (j=-n;j<=n;j++)
	        sum += h[(int) LINCOOR(r,BC(c+j),width)]*mask[(int) j+n];
	      k[(int) l] = (float)sum;
	    }
	  }
	  for (r=0; r<height; r++) {
	    for (c=width-n; c<width; c++) {
	      l = LINCOOR(r,c,width);
	      sum = 0.0;
	      for (j=-n;j<=n;j++)
	        sum += h[(int) LINCOOR(r,BC(c+j),width)]*mask[(int) j+n];
	      k[(int) l] = (float)sum;
	    }
	  }
	  
	  return k;
	}


/** Convolve an image with a derivative of the Gaussian. */
public double [] convolve_gauss(double[] image, double sigma, int deriv_type)
{
  double []  hr;
  double []  hc;
  
  
  double [] h;
  double [] k;

  switch (deriv_type) {
    case DERIV_R:
      hr = compute_gauss_mask_1(sigma);
      hc = compute_gauss_mask_0(sigma);
      break;
    case DERIV_C:
      hr = compute_gauss_mask_0(sigma);
      hc = compute_gauss_mask_1(sigma);
      break;
    case DERIV_RR:
      hr = compute_gauss_mask_2(sigma);
      hc = compute_gauss_mask_0(sigma);
      break;
    case DERIV_RC:
      hr = compute_gauss_mask_1(sigma);
      hc = compute_gauss_mask_1(sigma);
      break;
    case DERIV_CC:
      hr = compute_gauss_mask_0(sigma);
      hc = compute_gauss_mask_2(sigma);
      break;
    default: //just a stub
      hr = compute_gauss_mask_0(sigma);
      hc = compute_gauss_mask_0(sigma);
      break;
  }

  //maskr = hr + nr;
  //maskc = hc + nc;

  h = convolve_rows_gauss(image,hr);
  k = convolve_cols_gauss(h,hc);

  return k;
}
	
public double[][] get_all_derivatives(double[] image, double sigma) 
{
	double [][] k;
	
	k = new double[5][(int) (width*height)];
	
	k[0] = convolve_gauss(image,sigma,DERIV_R);
	k[1] = convolve_gauss(image,sigma,DERIV_C);
	k[2] = convolve_gauss(image,sigma,DERIV_RR);
	k[3] = convolve_gauss(image,sigma,DERIV_RC);
	k[4] = convolve_gauss(image,sigma,DERIV_CC);
	
	return k;
	
}
	
	
	
	
	
	 /** Translate row and column coordinates of an image into an index into its
    one-dimensional array. */
	public  long LINCOOR(long row, long col,long width){ 
		return (long)((row)*(width)+(col));				
	}
	/** Mirror the row coordinate at the borders of the image; height must be a
	   defined variable in the calling function containing the image height. */
	public  long BR(long row) {
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
	public  long BC(long col) { 
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
