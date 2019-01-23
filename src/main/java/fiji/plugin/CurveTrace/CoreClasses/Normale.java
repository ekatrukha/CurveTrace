package fiji.plugin.CurveTrace.CoreClasses;

/**class for normal angle at the point along the curve.
 * contains absolute angle of normal vector
 * and operations with it (averaging, distance, etc).
 * It is a circular variable (from -pi/2 to pi/2),
 * so it should be handled carefully.
 * 
 * Adopted from:
 * Circular Values Math and Statistics with C++11
 * Lior Kogan **/
public class Normale {
	
	/**left range of normale, range being [L,H) **/
	public static final double L = Math.PI*(-0.5);
	/**right range of normale, range being [L,H)**/
	public static final double H = Math.PI*(0.5);
	/**zero value**/
	public static final double Z = 0.0;
	/**range**/
	public static final double R = Math.PI;
	/**half range**/
	public static final double R_2 = Math.PI*0.5;
	
	public double val;
	
	public Normale(double in)
	{
		val = Normale.Wrap(in);
	}
	
	public double toDouble()
	{
		return val;
	}
	public float toFloat()
	{
		return (float)val;
	}
	

	/**Returns a Normale whose value is (this + sec).**/
	public Normale add(Normale sec)
	{
		return new Normale(val+sec.val- Normale.Z);
		
	}

	/**calculates average between two angles.
	 * in general, it should be two values: returned and
	 * (n2.val+0.5*Sdist(n2,n1), but they differ only when
	 * Sdist (n1,n2)==Sdist (n2,n1) = -R_2;
	 * This is a special case and we here ignore it 
	 * (should not be so crucial for lines
	 * **/
	public static Normale average(Normale n1, Normale n2)
	{
		return new Normale(n1.val+0.5*Sdist(n1,n2));
	}
	/**weighted
	 * **/
	public static Normale averageWeighted(Normale n1, Normale n2, double w)
	{
		return new Normale(n1.val+(1.0-w)*Sdist(n1,n2));
	}
	/**Carefully made Floating-point modulo function: 
	 * reminder = mod(dividend, divisor)
	 * taking into account limited accuracy of floating point**/
	static double Mod(double x, double y)
	{
		if (0.0 == y)
	        return x;

	    double m= x - y * Math.floor(x/y);
	    if (y > 0)              // modulo range: [0..y)
	    {
	        if (m>=y)           // Mod(-1e-16             , 360.    ): m= 360.
	            return 0;

	        if (m<0 )
	        {
	            if (y+m == y)
	                return 0  ; // just in case...
	            else
	                return y+m; // Mod(106.81415022205296 , _TWO_PI ): m= -1.421e-14 
	        }
	    }
	    else                    // modulo range: (y..0]
	    {
	        if (m<=y)           // Mod(1e-16              , -360.   ): m= -360.
	            return 0;

	        if (m>0 )
	        {
	            if (y+m == y)
	                return 0  ; // just in case...
	            else
	                return y+m; // Mod(-106.81415022205296, -_TWO_PI): m= 1.421e-14 
	        }
	    }
	    
	    return m;
		
	}
	/** whether the value is in range **/
	static boolean isInRange(double r)
	{
	    return (r>=Normale.L && r<Normale.H);
	}
	
	/** 'wraps' circular-value to [L,H) range **/
	static double Wrap(double r)
	{
	    // the next lines are for optimization and improved accuracy only
	    if (r>=Normale.L)
	    {
	             if (r< Normale.H        ) 
	            	 return r        ;
	        else if (r< Normale.H+Normale.R) 
	        	return r-Normale.R;
	    }
	    else
	             if (r>=Normale.L-Normale.R) 
	            	 return r+Normale.R;
	 
	    // general case
	    return Normale.Mod(r - Normale.L,Normale.R) + Normale.L;
	}
	/** the length of the directed walk from c1 to c2, with the lowest absolute-value length
	 return value is in [-R/2, R/2) **/
	public static double Sdist(Normale c1, Normale c2)
	{
	    double d= c2.val-c1.val;
	    if (d <  -Normale.R_2) 
	    	return d + Normale.R;
	    if (d >=  Normale.R_2) 
	    	return d - Normale.R;
	    return d;
	}
	/** the length of the increasing walk from c1 to c2 with the lowest length
	 return value is in [0, R) **/
	static double Pdist(Normale c1, Normale c2)
	{
	    return c2.val>=c1.val ? c2.val-c1.val : Normale.R-c1.val+c2.val;
	}
}
