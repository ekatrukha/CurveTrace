package fiji.plugin.CurveTrace.Fit;

import jaolho.data.lma.LMAFunction;

public class OneDGaussianBGMeanFixed extends LMAFunction {

	
	public double getY(double x, double[] a) {
		return Math.abs(a[2])+a[0] * Math.exp(-0.5* ( (x)*(x)/(a[1]*a[1]) ) ) ;
	}

	@Override
	public double getPartialDerivate(double x, double[] a, int parameterIndex) {
		switch (parameterIndex) {
			case 0: return Math.exp(-0.5* ( (x)*(x)/(a[1]*a[1]) ) );
			//case 1: return a[0] * Math.exp(-0.5* ( (x-dMean)*(x-dMean)/(a[2]*a[2]) ) ) * ((x-dMean)/(a[2]*a[2]));
			case 1: return a[0] * Math.exp(-0.5* ( (x)*(x)/(a[1]*a[1]) ) ) * ((x)*(x)/(a[1]*a[1]*a[1]));
			case 2: return 1.0;
		}
		throw new RuntimeException("No such parameter index: " + parameterIndex);
	}

}
