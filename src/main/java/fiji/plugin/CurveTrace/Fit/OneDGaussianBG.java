package fiji.plugin.CurveTrace.Fit;

import jaolho.data.lma.LMAFunction;

public class OneDGaussianBG extends LMAFunction 
{
	@Override
	public double getY(double x, double[] a) {
		return Math.abs(a[3])+a[0] * Math.exp(-0.5* ( (x-a[1])*(x-a[1])/(a[2]*a[2]) ) ) ;
	}
	@Override
	public double getPartialDerivate(double x, double[] a, int parameterIndex) {
		switch (parameterIndex) {
			case 0: return Math.exp(-0.5* ( (x-a[1])*(x-a[1])/(a[2]*a[2]) ) );
			case 1: return a[0] * Math.exp(-0.5* ( (x-a[1])*(x-a[1])/(a[2]*a[2]) ) ) * ((x-a[1])/(a[2]*a[2]));
			case 2: return a[0] * Math.exp(-0.5* ( (x-a[1])*(x-a[1])/(a[2]*a[2]) ) ) * ((x-a[1])*(x-a[1])/(a[2]*a[2]*a[2]));
			case 3: return 1.0;
		}
		throw new RuntimeException("No such parameter index: " + parameterIndex);
	}
}
