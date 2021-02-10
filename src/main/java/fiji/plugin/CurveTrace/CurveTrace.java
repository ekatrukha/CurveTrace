package fiji.plugin.CurveTrace;

import java.awt.Color;

import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.MultiDirectionalSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;

import fiji.plugin.CurveTrace.stegers.Steger_tracer;
import fiji.plugin.CurveTrace.stegers.contour;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.Arrow;
import ij.gui.GenericDialog;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;

public class CurveTrace implements PlugIn{

	/** main image window */
	public ImagePlus imp; 	
	/** currently active processor */
	public ImageProcessor ip; 
	/** current slice in stack */
	int nCurrentSlice;
	/** current slice (frame) */
	int nFrame = 0;
	/** total number of frames in image */
	int nStackSize;
	
	/* curve tracer */
	Steger_tracer ctracer;
	
	/** Total number of contours found */	
	int nContNumTotal = 0;
	
	/** link to roi manager */
	RoiManager roi_manager;
	/** overlay of main image */
	Overlay image_overlay; 
	
	/** Results table */
	ResultsTable ptable = ResultsTable.getResultsTable();
	/** whether results table had been cleaned*/
	boolean bClearResults = false;
	
	/** Type of lines added to RoiManager and Overlay*/
	public int nRoiLinesType = 1;
	/** Color of lines added to RoiManager and Overlay*/
	public int nRoiLinesColor = 0;
	
	
	public boolean DEBUG_MODE = false;
	public boolean DEBUG_show_subpixel = false;
	public boolean DEBUG_show_tracing = false;
	public boolean DEBUG_show_extensions = false;
	public boolean DEBUG_show_eigenvector = false;
	public boolean DEBUG_show_eigenvalue = false;
	public boolean DEBUG_show_junctions = false;
	
	@Override
	public void run(String arg) {
		int i;
		
		double [] bounds = new double [2];
		// TODO Auto-generated method stub
		// get active image
		imp = IJ.getImage();
		if(null == imp)
		{
			IJ.noImage();	    
			return;
		}		
		else if (imp.getType() != ImagePlus.GRAY8 && imp.getType() != ImagePlus.GRAY16 && imp.getType() != ImagePlus.GRAY32) 
		{
		    IJ.error("8, 16 or 32-bit grayscale image required");
		    return;
		}
		
		
		
		nCurrentSlice = imp.getCurrentSlice();
		ctracer = new Steger_tracer();
		
		// get active imageprocessor
		ctracer.ip = imp.getProcessor();
		
		nStackSize = imp.getStackSize();
		//width and height of image
		ctracer.width = ctracer.ip.getWidth();
		ctracer.height = ctracer.ip.getHeight();
		
		//show parameters dialog
		if(!show_parameters_dialog())
			return;
		ctracer.compute_line_points();
		
		if(!ctracer.get_threshold_levels())
			return ;
		
		
		
		IJ.showStatus("Finding lines...");
		IJ.showProgress(0, nStackSize);
		
		image_overlay = new Overlay();
		for(nFrame=0; nFrame<nStackSize; nFrame++)
		{
			imp.setSliceWithoutUpdate(nFrame+1);
			ctracer.ip = imp.getProcessor();
			ctracer.compute_line_points();
			bounds[0]=ctracer.low;
			bounds[1]=ctracer.high;
			ctracer.getismax(bounds); //move it to compute contours function

			
			
			if(DEBUG_show_subpixel)
				add_subpixelxy();
			
			ctracer.compute_contours();
			
			if(ctracer.compute_width)
			{
				ctracer.compute_line_width();
			}
				

			if(DEBUG_show_eigenvector)
				add_eigenvectors();
			if(DEBUG_show_junctions)
				add_junctions();
			
			
			/*
			//AUTOMATIC
			//STUB DEBUG
			ctracer.prepareOrigMSE();
			
			//IJ.log("MSE window_params: "+Double.toString(ctracer.value(bounds)));
			IJ.log("MSE manual_params: "+Double.toString(ctracer.value(new double[] { bounds[0], bounds[1]})));
			

			//ctracer.calcMap();
			//new ImagePlus("MSE_map", ctracer.ipMSE.duplicate()).show();
			//new ImagePlus("CurveN_map", ctracer.ipCurveN.duplicate()).show();
			bounds=optimizetest();
			ptable.reset();
			ptable.setPrecision(16);
			Double [] val;
			for (i=0;i<ctracer.ptable.size();i++)
			{
				val=ctracer.ptable.get(i);
				ptable.incrementCounter();
				ptable.addValue("EvalN", val[0]);
				ptable.addValue("MSE", val[1]);
				ptable.addValue("LB", val[2]);
				ptable.addValue("UB", val[3]);
				ptable.addValue("CurveN", val[4]);
			}*/
			
			//show_contours();
			nContNumTotal+= ctracer.cont.size();
			//clear results table if we found anything
			if (nContNumTotal>0 && !this.bClearResults)
			{
				ptable = ResultsTable.getResultsTable();
				ptable.setPrecision(5);
				ptable.reset(); // erase results table
				this.bClearResults = true;
			}
			//add lines characteristics to results table
			if (nContNumTotal>0)
				add_results_table();
			
			
	
			show_contours();
			
			IJ.showProgress(nFrame, nStackSize);
		}//nFrame cycle
	
		IJ.showProgress(nStackSize, nStackSize);
		IJ.showStatus("Done.");
		//ptable.show("Results");
		ptable.show("Results");
		imp.setPosition(nCurrentSlice);
		imp.draw();
		imp.setActivated();
		
	}
	
	public double [] optimizetest()
	{
		
		double minSal=0;
		double maxSal=0; 
		int i;
		double valMin = Double.MAX_VALUE;
		double valCurr;
		int indMin=0;

		
		double lowerTh=0;
		double fractionUp=0.1;
		

		double [] vals;
		//find the best estimate of a threshold
		
		for(i=0;i<4;i++)
		{
			minSal = ctracer.threshold_est_min[i];
			maxSal = ctracer.threshold_est_min[i] + fractionUp*(ctracer.threshold_est_max - ctracer.threshold_est_min[i]);
			valCurr=ctracer.value(new double[] { minSal, maxSal});
			if(valCurr<valMin)
			{
				valMin=valCurr;
				lowerTh=ctracer.threshold_est_min[i];
				indMin=i;
			}
		}

		IJ.log("MinInd="+Integer.toString(indMin)+" lowerth="+Double.toString(lowerTh));
		minSal=lowerTh;
		maxSal = lowerTh + fractionUp*(ctracer.threshold_est_max - lowerTh);
		
		double stepSal = 0.05*(ctracer.threshold_est_max - lowerTh);

		SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
		//PowellOptimizer optimizer = new PowellOptimizer(1e-10, 1e-30);
		//BOBYQAOptimizer optimizer = new BOBYQAOptimizer(5,(maxSal-minSal),1e-8);
		
		
		IJ.log(String.format("MSE before:%f",ctracer.value(new double[] { minSal, maxSal})));
		IJ.log(String.format("Ini min:%f, max:%f, step:%f",minSal,maxSal,stepSal));
		PointValuePair optimum=optimizer.optimize(new MaxEval(100),
				new MaxIter(100),
                new ObjectiveFunction(ctracer),
                GoalType.MINIMIZE,
                new InitialGuess(new double[] { minSal, maxSal }),
                new MultiDirectionalSimplex(new double[] { stepSal, (-1)*stepSal }));
                //new SimpleBounds(new double[] { ctracer.threshold_tmp_ip_min, ctracer.threshold_tmp_ip_min}, new double[] { ctracer.threshold_tmp_ip_max, ctracer.threshold_tmp_ip_max }));
                //new NelderMeadSimplex(new double[] { stepSal, (-1.0)*stepSal }));
		IJ.log("converged in "+Integer.toString(optimizer.getIterations())+" iterations");
		IJ.log("and "+Integer.toString(optimizer.getEvaluations())+" function evaluations");
		vals=optimum.getPoint();
		if (vals[0]>vals[1])
		{	stepSal=vals[1];
			vals[1]=vals[0];
			vals[0]=stepSal;
		}
		if(vals[0]<ctracer.threshold_tmp_ip_min)
			vals[0]=ctracer.threshold_tmp_ip_min;
		if(vals[1]>ctracer.threshold_tmp_ip_max)
			vals[1]=ctracer.threshold_tmp_ip_max;
		IJ.log(String.format("MSE after:%f",optimum.getValue()));
		IJ.log(String.format("Fin min:%f, max:%f",vals[0],vals[1]));
		
		return vals;
	}

	public boolean show_parameters_dialog()
	{
				
		int nDetectionType;

		String [] DetectionType = new String [] {
				"White lines on dark background", "Dark lines on white background"};
		String [] RoiType = new String [] {
				"Nothing", "Only lines", "Lines of average width", "Lines with width border"};
		String [] RoiColor = new String [] {
				"Rainbow", "Black", "White", "Blue", "Cyan","Green","Magenta","Orange","Pink","Red","Yellow"};
		/*
		String [] sLabelsCheckbox = new String [] {
				"subpixel x,y (red circles)","eigenvector (green arrows)", 
				"tracing lines (yellow circles)","eigenvalues",
				"extensions (green circles)","junctions"};
		boolean [] sLabelsDefault = new boolean [] {
				Prefs.get("stegers_source.show_subpixel", DEBUG_show_subpixel), Prefs.get("stegers_source.show_eigenvector", DEBUG_show_eigenvector), 
				Prefs.get("stegers_source.show_tracing", DEBUG_show_tracing), Prefs.get("stegers_source.show_eigenvalue", DEBUG_show_eigenvalue),
				Prefs.get("stegers_source.show_extensions", DEBUG_show_extensions),Prefs.get("stegers_source.show_junctions", DEBUG_show_junctions), 
				};
		*/

		
		GenericDialog paramD = new GenericDialog("Parameters");
		paramD.addChoice("Detection type:", DetectionType, Prefs.get("stegers_source.mode_light", "White lines on dark background"));
		//paramD.setInsets(10, 50, 0); 
		paramD.addNumericField("Line width (SD of profile)", Prefs.get("stegers_source.sigma", ctracer.SIGMA), 2, 4,"pixels");
		//paramD.setInsets(10, 100, 0);
		paramD.addNumericField("Maximum angle difference", Prefs.get("stegers_source.max_angle_difference", ctracer.MAX_ANGLE_DIFFERENCE*180/Math.PI), 1,4,"degrees");
		paramD.setInsets(10, 50, 0); 
		paramD.addCheckbox("Extend line ends", Prefs.get("stegers_source.extend_lines", ctracer.extend_lines));
		paramD.addNumericField("Maximum line extension", Prefs.get("stegers_source.max_line_extension", ctracer.MAX_LINE_EXTENSION/ctracer.SIGMA), 1,4,"*line width, pixels");
		paramD.setInsets(10, 50, 0);
		paramD.addCheckbox("Split lines at junctions", Prefs.get("stegers_source.split_lines", ctracer.split_lines));
		paramD.addNumericField("Minimum number of points in line", Prefs.get("stegers_source.min_nodes", ctracer.nMinNumberOfNodes), 0,3,"");
		paramD.setInsets(10, 50, 0);
		paramD.addCheckbox("Correct line position", Prefs.get("stegers_source.correct_pos", ctracer.correct_pos));
		paramD.setInsets(10, 50, 0);
		paramD.addCheckbox("Compute line width", Prefs.get("stegers_source.compute_width", ctracer.compute_width));
		paramD.addNumericField("Maximum width search", Prefs.get("stegers_source.max_line_width", ctracer.MAX_LINE_WIDTH/ctracer.SIGMA), 1,4,"*line width, pixels");
		paramD.setInsets(10, 0, 0);
		//MAX_LINE_WIDTH = (2.5*SIGMA);
		paramD.addChoice("Add to overlay and RoiManager:", RoiType, Prefs.get("stegers_source.roitype", "Only lines"));
		paramD.addChoice("Color of added lines:", RoiColor, Prefs.get("stegers_source.roicolor", "Rainbow"));
		
		//paramD.addMessage("~~~~~~~~~~~~ Learning/Debug  ~~~~~~~~~~~~");
		//paramD.addCheckboxGroup(3,3,sLabelsCheckbox,sLabelsDefault);
		

		paramD.setResizable(false);
		paramD.showDialog();
		if (paramD.wasCanceled())
            return false;
		
		nDetectionType = paramD.getNextChoiceIndex();
		Prefs.set("stegers_source.mode_light", DetectionType[nDetectionType]);
		if(nDetectionType ==0)
			ctracer.MODE_LIGHT = true;
		else
			ctracer.MODE_LIGHT = false;
		
		ctracer.SIGMA = paramD.getNextNumber();
		Prefs.set("stegers_source.sigma", ctracer.SIGMA);
		
		ctracer.MAX_ANGLE_DIFFERENCE = paramD.getNextNumber()*Math.PI/180;
		Prefs.set("stegers_source.max_angle_difference", ctracer.MAX_ANGLE_DIFFERENCE*180/Math.PI);		
				
		ctracer.extend_lines = paramD.getNextBoolean();
		Prefs.set("stegers_source.extend_lines", ctracer.extend_lines);
		
		ctracer.MAX_LINE_EXTENSION = paramD.getNextNumber()*ctracer.SIGMA;
		Prefs.set("stegers_source.max_line_extension", ctracer.MAX_LINE_EXTENSION/ctracer.SIGMA);	
		
		ctracer.split_lines = paramD.getNextBoolean();
		Prefs.set("stegers_source.split_lines", ctracer.split_lines);
		
		ctracer.nMinNumberOfNodes = (int) paramD.getNextNumber();
		Prefs.set("stegers_source.min_nodes", ctracer.nMinNumberOfNodes);
		
		ctracer.correct_pos = paramD.getNextBoolean();
		Prefs.set("stegers_source.correct_pos", ctracer.correct_pos);

		ctracer.compute_width = paramD.getNextBoolean();
		Prefs.set("stegers_source.compute_width", ctracer.compute_width);

		ctracer.MAX_LINE_WIDTH = paramD.getNextNumber()*ctracer.SIGMA;
		Prefs.set("stegers_source.max_line_width", ctracer.MAX_LINE_WIDTH/ctracer.SIGMA);		
		
		nRoiLinesType = paramD.getNextChoiceIndex();
		Prefs.set("stegers_source.roitype", RoiType[nRoiLinesType]);

		nRoiLinesColor = paramD.getNextChoiceIndex();
		Prefs.set("stegers_source.roicolor", RoiColor[nRoiLinesColor]);
		
		/*
		DEBUG_show_subpixel = paramD.getNextBoolean();
		Prefs.set("stegers_source.show_subpixel", DEBUG_show_subpixel);

		DEBUG_show_eigenvector = paramD.getNextBoolean();
		Prefs.set("stegers_source.show_eigenvector", DEBUG_show_eigenvector);
		
		DEBUG_show_tracing = paramD.getNextBoolean();
		Prefs.set("stegers_source.show_tracing", DEBUG_show_tracing);
		
		DEBUG_show_eigenvalue = paramD.getNextBoolean();
		Prefs.set("stegers_source.show_eigenvalue", DEBUG_show_eigenvalue);
		
		DEBUG_show_extensions = paramD.getNextBoolean();
		Prefs.set("stegers_source.show_extensions", DEBUG_show_extensions);
		
		DEBUG_show_junctions = paramD.getNextBoolean();
		Prefs.set("stegers_source.show_junctions", DEBUG_show_junctions);

		*/
		return true;
	
	}
	
	/**
	 * Function adding all found contours to overlay and roi manager
	 * */
	public void show_contours()
	{
			
		  PolygonRoi polyline_p;
		  contour tmp_cont;
		  long num_pnt;
		  String Roi_name;
		  int i,j;
		  
		  int color_n = 8;
		  double nx,ny;
		  double line_avg_width;
		
		  Color[] colors;
		  Color[] colors_bord;
		  Color single_main;
		  Color single_border;
		  colors = new Color[8];		  
		  colors_bord = new Color[8];
		  colors[0] = Color.BLUE;
		  colors[1] = Color.CYAN;
		  colors[2] = Color.GREEN;
		  colors[3] = Color.MAGENTA;
		  colors[4] = Color.ORANGE;
		  colors[5] = Color.PINK;
		  colors[6] = Color.RED;
		  colors[7] = Color.YELLOW;
		  
		  colors_bord[0] = new Color(0,0,139);//dark blue
		  colors_bord[1] = new Color(0,139,139); //dark cyan
		  colors_bord[2] = new Color(0,100,0); //dark green
		  colors_bord[3] = new Color(139,0,139); //dark magenta
		  colors_bord[4] = new Color(255,140,0); //dark orange
		  colors_bord[5] = new Color(199,21,133); //mediumvioletred
		  colors_bord[6] = new Color(139,0,0); //dark red
		  colors_bord[7] = new Color(153,153,0); //dark yellow
		  
		  
		  Roi_name  = "contour_f_";
		  Roi_name = Roi_name.concat(Integer.toString(nFrame+1));
		  Roi_name = Roi_name.concat("_c_");
		  if(nRoiLinesType>0)
		  {
			  //black
			  if(nRoiLinesColor==1)
			  {
				  for(i=0;i<8;i++)
				  {
					  colors[i]=Color.BLACK;
					  colors_bord[i]=Color.GRAY;
				  }
			  }
			  //black
			  if(nRoiLinesColor==2)
			  {
				  for(i=0;i<8;i++)
				  {
					  colors[i]=Color.WHITE;
					  colors_bord[i]=Color.GRAY;
				  }
			  }
			  if(nRoiLinesColor>=3)
			  {
				  single_main = colors[nRoiLinesColor-3];
				  single_border = colors_bord[nRoiLinesColor-3];
				  for(i=0;i<8;i++)
				  {
					  colors[i]=single_main;
					  colors_bord[i]=single_border;
				  }
			  }		  
			  
	  		  roi_manager = RoiManager.getInstance();
			  if(roi_manager == null) roi_manager = new RoiManager();
			  
			  
			  //show contour and add it to Roi manager		
			  for (i=0; i<ctracer.cont.size(); i++) 
			  {
				  
				  tmp_cont = ctracer.cont.get(i);  
				  num_pnt = tmp_cont.num;
				  float[] pxs = new float[(int) num_pnt];
				  float[] pxy = new float[(int) num_pnt];
				  for(j=0;j<num_pnt;j++)
				  {
					  pxs[j]=tmp_cont.row.get(j);
					  pxy[j]=tmp_cont.col.get(j);
				  }
					  			  
				  polyline_p = new PolygonRoi(pxy, pxs, Roi.POLYLINE);
				  polyline_p.setPosition(nFrame+1);
				  polyline_p.setStrokeColor(colors[i%color_n]);
		
				  polyline_p.setStrokeWidth(0.0);
				  polyline_p.setName(Roi_name.concat(Integer.toString(i+1)));
				  
				  if(nRoiLinesType == 1 || nRoiLinesType == 3)
					  	roi_manager.addRoi(polyline_p);
				  
				  if(ctracer.compute_width && nRoiLinesType ==2)
				  {
					  line_avg_width = 0;
					  for(j=0;j<num_pnt;j++)
					  {
						  line_avg_width +=tmp_cont.width_l.get(j)+tmp_cont.width_r.get(j);
					  }				 
					  line_avg_width =line_avg_width /num_pnt;
					 //polyline_p.fitSplineForStraightening();
					  polyline_p.setStrokeWidth(line_avg_width);
					  polyline_p.updateWideLine((float) line_avg_width);
				
	
					  //imp.setRoi(polyline_p);
					  //imp.getRoi();
					  
					  //roi_manager.runCommand("", colors[i%color_n].toString(), line_avg_width);
					  roi_manager.addRoi(polyline_p);
				  }
				  
				  
				  
				  image_overlay.add(polyline_p);
				  //polyline_p.setStrokeWidth(0.0);
				  
				  //add borders along width of line
				  if(ctracer.compute_width && nRoiLinesType ==3)
				  {
					  for(j=0;j<num_pnt;j++)
					  {
						nx = Math.cos(tmp_cont.angle.get(j));
				      	ny = Math.sin(tmp_cont.angle.get(j));		
					  	pxs[j]=(float) (tmp_cont.row.get(j) + nx*tmp_cont.width_l.get(j));
					  	pxy[j]=(float) (tmp_cont.col.get(j) + ny*tmp_cont.width_l.get(j));
					  }
					  polyline_p = new PolygonRoi(pxy, pxs, Roi.POLYLINE);
					  polyline_p.setStrokeColor(colors_bord[i%color_n]);
					  polyline_p.setPosition(nFrame+1);
					  polyline_p.setStrokeWidth(0.1);
					  image_overlay.add(polyline_p);				  
	
					  for(j=0;j<num_pnt;j++)
					  {
						nx = Math.cos(tmp_cont.angle.get(j));
				      	ny = Math.sin(tmp_cont.angle.get(j));
					  	pxs[j]=(float) (tmp_cont.row.get(j) - nx*tmp_cont.width_r.get(j));
					  	pxy[j]=(float) (tmp_cont.col.get(j) - ny*tmp_cont.width_r.get(j));
					  }
					  polyline_p = new PolygonRoi(pxy, pxs, Roi.POLYLINE);
					  polyline_p.setStrokeColor(colors_bord[i%color_n]);		
					  polyline_p.setStrokeWidth(0.1);
					  polyline_p.setPosition(nFrame+1);
					  image_overlay.add(polyline_p);				  
					  
				  }
					  			  			  
			  }//*/
			  imp.setOverlay(image_overlay);
			  imp.updateAndRepaintWindow();
			  imp.show();		  
			  
			  roi_manager.setVisible(true);
			  roi_manager.toFront();
		  }
	}
	
	/**
	 * Function adding all found contours to results table
	 * */
	public void add_results_table()
	{
		
		int i,j, num_pnt;
		contour tmp_cont;
		String Roi_name ="contour_f_";
		String Cont_name;
		Roi_name = Roi_name.concat(Integer.toString(nFrame+1));
		Roi_name = Roi_name.concat("_c_");

		for (i=0; i<ctracer.cont.size(); i++) 
		  {
			  Cont_name = Roi_name.concat(Integer.toString(i+1));
			  tmp_cont = ctracer.cont.get(i);  
			  num_pnt = (int) tmp_cont.num;
			  for(j=0;j<num_pnt;j++)
			  {
				  ptable.incrementCounter();
				  ptable.addLabel(Cont_name);
				  ptable.addValue("Frame Number", nFrame+1);
				  ptable.addValue("Contour Number", i+1);
				  ptable.addValue("Point Number", j+1);
				  ptable.addValue("X_(px)",tmp_cont.col.get(j));
				  ptable.addValue("Y_(px)",tmp_cont.row.get(j));
				  ptable.addValue("Angle_of_normal_(radians)",tmp_cont.angle.get(j));
				  ptable.addValue("Response",tmp_cont.response.get(j));
				  if(ctracer.compute_width)
				  {
					  ptable.addValue("Width_left_(px)",tmp_cont.width_l.get(j));
					  ptable.addValue("Width_right_(px)",tmp_cont.width_r.get(j));
					  if (ctracer.correct_pos)
					  {
						  ptable.addValue("Assymetry",tmp_cont.asymmetry.get(j));
						  //does not make much sense
						  ptable.addValue("Contrast",tmp_cont.contrast.get(j));
					  }
				  }
				  
			  }
		  }

	}
	
	
	/** Function adding to Roi subpixel position of line points */
	private void add_subpixelxy() 
	{
		
		
		OvalRoi resolved_p;

		int l;
			
		
		 for (int zzpx=0; zzpx<ctracer.width; zzpx++) 
		    {
				for (int zzpy=0; zzpy<ctracer.height; zzpy++) 
				{
					l = (int) Steger_tracer.LINCOOR(zzpy,zzpx,ctracer.width);
			    	
					//* show super-resolved lines positions
			    	 			    	 
					 resolved_p = new OvalRoi(ctracer.posy[l]-0.25,ctracer.posx[l]-0.25, 0.5, 0.5);
					 resolved_p.setPosition(nFrame+1);
			    	 resolved_p.setStrokeColor(Color.RED);
			    	 image_overlay.add(resolved_p);	
				}
		    }
		
		
	}
	/** Function adding to Roi subpixel position of line points */
	private void add_junctions() 
	{
		
		Roi junc_roi;
		
		for(int i=0; i<ctracer.junc.size();i++)
		{
			junc_roi = new Roi(ctracer.junc.get(i).y-0.5,ctracer.junc.get(i).x-0.5,1,1);
			junc_roi.setPosition(nFrame+1);
			junc_roi.setStrokeColor(Color.MAGENTA);
	    	image_overlay.add(junc_roi);	
			
		}
		
	}	
	/** Function adding to Roi transverse eigenvectors at each pixel */
	private void add_eigenvectors() 
	{
		
		
		Arrow eigenvector_arr;

		int l;
			
		float[] dxx =new float[2];
		float[] dyy =new float[2];
		 for (int zzpx=0; zzpx<ctracer.width; zzpx++) 
		    {
				for (int zzpy=0; zzpy<ctracer.height; zzpy++) 
				{
					l = (int) Steger_tracer.LINCOOR(zzpy,zzpx,ctracer.width);
			    	
			    	//* Show vectors directions
			    	 dxx[0] = (float) (zzpx+0.5*ctracer.normx[l]+0.5);
			    	 dxx[1] = (float) (zzpx-0.5*ctracer.normx[l]+0.5);
			    	 dyy[0] = (float) (zzpy-0.5*ctracer.normy[l]+0.5);
			    	 dyy[1] = (float) (zzpy+0.5*ctracer.normy[l]+0.5);
			    	 
			    	 eigenvector_arr = new Arrow(dxx[0], dyy[0], dxx[1], dyy[1]);		    	
					 eigenvector_arr.setStrokeColor(Color.GREEN);
					 eigenvector_arr.setPosition(nFrame+1);
					 eigenvector_arr.setStrokeWidth(0.05);
					 eigenvector_arr.setHeadSize(0.5);
					 image_overlay.add(eigenvector_arr);
				}
		    }
		
		
	}
}
