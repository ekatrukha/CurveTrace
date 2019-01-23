package fiji.plugin.CurveTrace.CoreClasses;

import fiji.plugin.CurveTrace.KDTree.IPoint;

/**
 *	Point class, contains coordinates and features of contour points
 */

public class Point implements Cloneable, IPoint {

	/**X and Y coordinate of a point **/
	public float[] coords;
	/**angle of normal (measured from the Y axis) in radians **/
	public Normale angle;
	/** width of curve at this point (SD of Gaussian if fitted) **/
	public float width;
	/** Gaussian amplitude at the point (if fitted) **/
	public float amp;
	/** Background value from fitted Gaussian (if fitted) **/
	public float bg;
	/** R2 of Gaussian fitting**/
	public float R2;
	/**Integrated intensity calculated from fitting**/
	public float integr_int;
	/** response from Steger's detection fitting **/
	public float response;

	
	//default constructor
	public Point()
	{
			coords = new float[2];
	}
	//constructor #1
	public Point(float [] coords_in, float angle_in, float width_in, float amp_in, float bg_in, float R2_in, float integr_int_in, float response_in)
	{
		coords = new float[2];
		coords[0]=coords_in[0];
		coords[1]=coords_in[1];
		angle= new Normale(angle_in);		
		width=width_in;
		amp=amp_in;
		bg=bg_in;
		R2=R2_in;
		integr_int= integr_int_in;
		response = response_in;
	}
	//constructor #2
	public Point(float x_in, float y_in, float angle_in, float response_in)
	{
		coords = new float[2];
		coords[0]=x_in;
		coords[1]=y_in;
		angle = new Normale(angle_in);
		response = response_in;
	}
	public Point clone()
	{
		Point pointret = new Point(coords, (float)angle.val, width, amp, bg, R2,integr_int,response);
		return pointret;
	}
	
	/** calculates average (in coordinates) point between
	 * current and provided **/
	public Point averagePoint(Point p2)
	{
		Point ave = new Point();
		//float ressum;
		
		//simple average
		ave.coords[0]=(float) (0.5*(coords[0]+p2.coords[0]));
		ave.coords[1]=(float) (0.5*(coords[1]+p2.coords[1]));
		
		//response weighted average
		//ressum = response + p2.response;
		//ave.coords[0]=(float) ((coords[0]*response+p2.coords[0]*p2.response)/ressum);
		//ave.coords[1]=(float) ((coords[1]*response+p2.coords[1]*p2.response)/ressum);

		// average angle
		// it should take into account circular manner
		// plus PI symmetry in the normale
		/*
		angle=(float) (angle%Math.PI);
		p2.angle=(float) (p2.angle%Math.PI);
		if(angle>0.5*Math.PI)
			angle=(float) (angle-Math.PI);
		if(angle<-0.5*Math.PI)
			angle=(float) (angle+Math.PI);
		if(p2.angle>0.5*Math.PI)
			p2.angle=(float) (p2.angle-Math.PI);	
		if(p2.angle<-0.5*Math.PI)
			p2.angle=(float) (p2.angle+Math.PI);	
		if(Math.abs(angle-p2.angle)>0.5*Math.PI)
			ave.angle = (float) (0.5*(angle+p2.angle+Math.PI));
		else
			ave.angle = (float) (0.5*(angle+p2.angle));
			*/
		ave.angle=Normale.average(this.angle,p2.angle);
		
		// for now just calculate average of everything
		
		ave.amp = (float) (0.5*(amp+p2.amp));
		
		ave.bg = (float) (0.5*(bg+p2.bg));
		ave.R2 = (float) (0.5*(R2+p2.R2));
		ave.integr_int = (float) (0.5*(integr_int+p2.integr_int));
		ave.response = (float) (0.5*(response+p2.response));
		return ave;		
	}
	/** calculates weighted average (in coordinates) point between
	 * current and provided, w is weight of current point, p2 gets 1-w **/
	public Point averageWeightedPoint(Point p2, float w)
	{
		Point ave = new Point();
		float w2=(float) (1.-w);
		
		//simple average
		ave.coords[0]=(float) ((w*coords[0]+w2*p2.coords[0]));
		ave.coords[1]=(float) ((w*coords[1]+w2*p2.coords[1]));
		
		//response weighted average
		//ressum = response + p2.response;
		//ave.coords[0]=(float) ((coords[0]*response+p2.coords[0]*p2.response)/ressum);
		//ave.coords[1]=(float) ((coords[1]*response+p2.coords[1]*p2.response)/ressum);

		// average angle
		// it should take into account circular manner
		// plus PI symmetry in the normale
		/*
		angle=(float) (angle%Math.PI);
		p2.angle=(float) (p2.angle%Math.PI);
		if(angle>0.5*Math.PI)
			angle=(float) (angle-Math.PI);
		if(p2.angle>0.5*Math.PI)
			p2.angle=(float) (p2.angle-Math.PI);		
		ave.angle = (float) ((w*angle+w2*p2.angle));
		*/
		ave.angle=Normale.averageWeighted(this.angle,p2.angle,w);
		
		// for now just calculate average of everything
		
		ave.amp = (float) ((w*amp+w2*p2.amp));
		
		ave.bg = (float) ((w*bg+w2*p2.bg));
		ave.R2 = (float) ((w*R2+w2*p2.R2));
		ave.integr_int = (float) ((w*integr_int+w2*p2.integr_int));
		ave.response = (float) ((w*response+w2*p2.response));
		return ave;		
		
	}
	@Override
	public double getX() {
		// TODO Auto-generated method stub
		return coords[0];
	}
	@Override
	public double getY() {
		// TODO Auto-generated method stub
		return coords[1];
	}
}
