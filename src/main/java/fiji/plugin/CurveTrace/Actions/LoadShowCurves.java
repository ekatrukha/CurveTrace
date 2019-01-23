package fiji.plugin.CurveTrace.Actions;

import fiji.plugin.CurveTrace.CoreClasses.CurveStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;



public class LoadShowCurves implements PlugIn
{
	
	public boolean bLSIgnoreFrame;
	public int nLS_OverlayLineType;
	
	@Override
	public void run(String arg0) 
	{
				
		/** main object storing all curves **/
		CurveStack curvestack = new CurveStack();
		
		if(!curvestack.loadCurvesFromRT())
			return;
		if(!loadShowCurvesDialog())
			return;
		
		if(!curvestack.showCurvesOverlay(nLS_OverlayLineType,bLSIgnoreFrame))
			return;
		
		
	}
	public boolean loadShowCurvesDialog()
	{
		GenericDialog loadshowD = new GenericDialog("Load curves from Results Table and add to overlay");
		String [] OverlayLineType = new String [] {
				"Nothing", "Only lines", "Lines of average width", "Lines with width border", "Lines with normal vectors"};
		
		loadshowD.addChoice("Add to overlay:", OverlayLineType, Prefs.get("CurveTrace.LS_OverlayLineType", "Only lines"));
		loadshowD.addCheckbox("Ignore frame number?", Prefs.get("CurveTrace.bLSIgnoreFrame", false));
		loadshowD.setResizable(false);
		loadshowD.showDialog();						
		if (loadshowD.wasCanceled())
            return false;
		nLS_OverlayLineType = loadshowD.getNextChoiceIndex();
		Prefs.set("CurveTrace.LS_OverlayLineType", OverlayLineType[nLS_OverlayLineType]);
		bLSIgnoreFrame = loadshowD.getNextBoolean();
		Prefs.set("CurveTrace.bLSIgnoreFrame", bLSIgnoreFrame);
		return true;
	}
	

}
