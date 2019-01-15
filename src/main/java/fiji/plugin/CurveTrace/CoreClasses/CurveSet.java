package fiji.plugin.CurveTrace.CoreClasses;

import java.awt.Color;
import java.awt.Font;
import java.util.ArrayList;

import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.TextRoi;
/**
 *	CurveSet class, contains array of Curves
 */
public class CurveSet extends ArrayList<Curve> 
{
	/**frame number in stack**/	
	public int nFrame;
	
	public boolean bSelfLinked;
	
	public CurveSet()
	{
		super();
	}
	/** function adds curve stored inside object to Overlay **/
	public Overlay addCurvesToOverlay(Color color, Overlay inoverlay)
	{
		Overlay image_overlay = inoverlay.duplicate();
		int j,k,pointsN;
		
		PolygonRoi polyline_p;
		TextRoi textRoiN;
		Curve curve;
		Font font;
		font = new Font("SansSerif",Font.PLAIN,6);
		
		//cycling through curveset
		for (j=0;j<this.size();j++)
		{
			 curve=this.get(j);
			 pointsN=curve.size();
			 float[] px = new float[pointsN];
			 float[] py = new float[pointsN];
			 
			 
			 for (k=0;k<pointsN;k++)
			 {
				 px[k]=curve.get(k).coords[0];
				 py[k]=curve.get(k).coords[1];
			 }
			 textRoiN = new TextRoi(px[0],py[1],Integer.toString(j+1),font);
			 textRoiN.setStrokeColor(color);
			 textRoiN.setPosition(nFrame);
			 image_overlay.add(textRoiN);
			 
			 polyline_p = new PolygonRoi(px, py, Roi.POLYLINE);
			 polyline_p.setPosition(nFrame);
			 polyline_p.setStrokeColor(color);
			 polyline_p.setStrokeWidth(0.0);
			 image_overlay.add(polyline_p);
		}
		return image_overlay;
	}
	public void addCurveSet(CurveSet addedcurveset)	
	{
		int i;
		for (i=0;i<addedcurveset.size();i++)
		{
			this.add(addedcurveset.get(i).duplicate());
		}
		
	}
}
