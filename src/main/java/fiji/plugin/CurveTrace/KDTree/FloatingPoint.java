package fiji.plugin.CurveTrace.KDTree;

/**
 * Provides a standard API for evaluating double numbers when dealing with
 * floating point rounding error.
 * <p>
 * Specifically, this class enables one to compare two double numbers to 
 * see if they are within an epsilon of each other. 
 * 
 * @author George Heineman
 * @version 1.0, 6/15/08
 * @since 1.0
 */
public class FloatingPoint {
	
	/** Numbers within this amount are considered to be the same. */
	public static final double epsilon = 1E-9;
	
	/** Deny any attempt to instantiate this class. */
	private FloatingPoint(){};
	
	/**
	 * Return -1 if d1 < d2, 0 if d1 == d2, or +1 if d1 > d2.
	 * 
	 * @param d1    first number against which to compare
	 * @param d2    second number being compared
	 */
	public static int compare (double d1, double d2) {
		if (FloatingPoint.lesser(d1, d2)) { return -1; }
		if (FloatingPoint.same(d1, d2)) { return 0; }
		
		return +1;
	}
	
	/**
	 * When value won't work, because numbers are potentially infinite, then
	 * use this one.
	 * <p>
	 * Standarard means for comparing double when dealing with the special
	 * quantities, NaN and infinite numbers. Also properly ensures that
	 * numbers "close to zero" (within an epsilon) are to be treated as zero
	 * for this computation.
	 * 
	 * @param d1   first number being compared    
	 * @param d2   second number being compared
	 */
	public static boolean same (double d1, double d2) {
		// NaN numbers cannot be compared with '==' and must be treated separately.
		if (Double.isNaN(d1)) {
			return Double.isNaN(d2);
		}
		
		// this covers Infinite and NaN cases.
		if (d1 == d2) return true;
		
		// Infinity values can be compared with '==' as above
		if (Double.isInfinite(d1)) {
			return false;
		}

	
		// try normal value
		return value (d1-d2) == 0;
		
	}
	
	/**
	 * See if the value is close enough to actually be considered 0.0 and
	 * return 0.0 if need be.
	 * <p>
	 * Otherwise the value is returned.
	 * 
	 * @param x   value being considered
	 */
	public static double value(double x) {
		if ((x >= 0) && (x <= epsilon)) {
			return 0.0;
		}
		
		if ((x < 0) && (-x <= epsilon)) {
			return 0.0;
		}
		
		return x;
	}

	/** Given closeness-to-epsilon, is x >= y? */
	public static boolean greaterEquals(double x, double y) {
		return value(x-y) >= 0;
	}
	
	/** Given closeness-to-epsilon, is x > y? */
	public static boolean greater(double x, double y) {
		return value(x-y) > 0;
	}
	
	/** Given closeness-to-epsilon, is x <= y? */
	public static boolean lesserEquals(double x, double y) {
		return value(x-y) <= 0;
	}
	
	/** Given closeness-to-epsilon, is x < y? */
	public static boolean lesser(double x, double y) {
		return value(x-y) < 0;
	}
}
