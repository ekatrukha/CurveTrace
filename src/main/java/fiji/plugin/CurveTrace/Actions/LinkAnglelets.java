package fiji.plugin.CurveTrace.Actions;

import java.awt.Color;
import java.util.ArrayList;


import fiji.plugin.CurveTrace.CoreClasses.Chain;
import fiji.plugin.CurveTrace.CoreClasses.Curve;
import fiji.plugin.CurveTrace.CoreClasses.CurveSet;
import fiji.plugin.CurveTrace.CoreClasses.CurveStack;
import fiji.plugin.CurveTrace.CoreClasses.CurvesOverlap;
import fiji.plugin.CurveTrace.CoreClasses.Point;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.ContrastEnhancer;
import ij.plugin.PlugIn;
import ij.plugin.ZProjector;
import ij.process.ImageProcessor;

public class LinkAnglelets implements PlugIn {
	
	/** main image window */
	ImagePlus imp; 	
	
	/** main object storing all curves **/
	CurveStack curvestack;
	
	/** linking result as set **/
	CurveSet finalcurveset;
	
	/** points below this distance a considered overlapping **/
	double dPointDistance;
	/** threshold difference in normale angle (more is not considered an overlap), radians **/
	double dAngleDistance;
	/**threshold of fraction of curve length in the FULL overlap,
	 * if for one of curves overlap is higher, curves will be merged **/
	double dFullOverlapF;
	/**relative distance from the edge of overlap
	 * till curve end in relative curve length to consider it EXTENSION ovelap **/
	double dExtensionOverlapF;
	
	/** minimum number of points to consider overlap **/
	int nOverLapThreshold;
	
	/** connectivity matrix, specifying which curves should be linked**/
	int [][] dConnectivityMatrix;
	
	/** array storing all current overlaps between two curvesets**/
	ArrayList<CurvesOverlap> alloverlaps;
	
	/** array containing all chains of curves for linking **/
	ArrayList<Chain> allchains;
	/**whether to show intermediate linking results **/
	boolean bLinkDebug;
	
	/**whether to generate cumulative SUM z-projection**/
	boolean bCumulProj;
	
	/** overlay of main image */
	Overlay image_overlay; 
	
	@Override
	public void run(String arg) {
		
		int i;
		

		
		CurveStack outputcurvestack = new CurveStack();

		finalcurveset = new CurveSet();

		/** input curves **/
		curvestack = new CurveStack();
		

		//load curves from results table
		if(!curvestack.loadCurvesFromRT())
			return;
		
		if(curvestack.size()<=1)
		{
		    IJ.error("Curve detection over multiple frames is required! Aborting.");
		    return;
		}
		
		//ask for fitting parameters
		if(!showLinkingDialog())
			return;
		
		
		finalcurveset =curvestack.get(0);
		
		image_overlay = new Overlay();
		if(bLinkDebug)
		{			
			imp = IJ.getImage();
			if(null == imp)
			{
				IJ.noImage();	    
				return;
			}	
			if(imp.getStackSize()!=curvestack.nFramesTotal)
			{
				IJ.error("Frame number of image is not the same as curve tracings!");
				return;
			}
			
			finalcurveset.nFrame=1;
			image_overlay=finalcurveset.addCurvesToOverlay(Color.GREEN, image_overlay);
			if(bCumulProj)
				constructZprojects();
		}
		IJ.showStatus("Linking curves...");
		IJ.showProgress(0, curvestack.size());
		for (i=1;i<curvestack.size();i++)
		{

			finalcurveset = linkTwoCurveSets(finalcurveset,curvestack.get(i));
			/*finalcurveset.bSelfLinked = false;
			while (!finalcurveset.bSelfLinked)
			{
				finalcurveset=selfLinkCurveSet(finalcurveset);
			}*/
			
			IJ.showProgress(i, curvestack.size());
			if(bLinkDebug)
			{
				curvestack.get(i).nFrame=i;
				image_overlay=curvestack.get(i).addCurvesToOverlay(Color.RED, image_overlay);
				finalcurveset.nFrame=i+1;
				image_overlay=finalcurveset.addCurvesToOverlay(Color.GREEN, image_overlay);
			}
			
		}
		IJ.showStatus("Linking curves...Done.");
		IJ.showProgress(curvestack.size(), curvestack.size());
		finalcurveset.nFrame=1;
		outputcurvestack.add(finalcurveset);
		outputcurvestack.nFramesTotal=1;
		if(!bLinkDebug)			
			outputcurvestack.showCurvesOverlay(1, true);
		else
		{
			imp.setOverlay(image_overlay);
			imp.updateAndRepaintWindow();
			imp.show();
		}
		//outputcurvestack.toResultsTable();
		
	}
	

	
	/** function linking two curvesets to one **/
	public CurveSet linkTwoCurveSets(CurveSet curveset1, CurveSet curveset2)
	{
		CurveSet mergedcs = new CurveSet();
		 
		int [] nCurves;
		/** marks if curve was processed by joining **/
		boolean [] curveschecked1,curveschecked2;
		int i,j;
		int nOverlapsCount=0;
		
		CurvesOverlap overlap;
		nCurves = new int [2];
		nCurves[0]=curveset1.size();
		nCurves[1]=curveset2.size();
		curveschecked1=new boolean[nCurves[0]];
		curveschecked2=new boolean[nCurves[1]];
		dConnectivityMatrix = new int [nCurves[0]][nCurves[1]];
		alloverlaps= new ArrayList<CurvesOverlap>();
		
		//checking for ovelaps between curves in two sets and 
		//filling connectivity matrix
		for (i=0;i<nCurves[0];i++)
		{
		
			for (j=0;j<nCurves[1];j++)
			{
				overlap = calculateOverlap(curveset1.get(i),curveset2.get(j));
				overlap.curves_ind[0]=i;
				overlap.curves_ind[1]=j;
				
				//DEBUG
				if (overlap.bDirectConflict)
				{
					IJ.log("--------CHAIN TO CHAIN BEG ---------");
					IJ.log("Frame "+Integer.toString(curveset1.nFrame) + " and "+Integer.toString(curveset2.nFrame));
					IJ.log("Curve "+Integer.toString(i+1) + " and "+Integer.toString(j+1));
					IJ.log("Overlap type "+Integer.toString(overlap.nOvelapType) );
					String indices = "indx ";
					for (int h=0;h<overlap.dirChange.size();h++)
						indices = indices +Integer.toString(overlap.dirChange.get(h))+" ";
					IJ.log("Indices "+indices + Integer.toString(overlap.nOverlapSize));
					IJ.log("--------CHAIN TO CHAIN END ---------");
				}
				
				//consider only full or extension overlaps
				if (overlap.nOvelapType<2)
				{
					//add to connectivity matrix
					nOverlapsCount++;
					dConnectivityMatrix[i][j]=nOverlapsCount;
					alloverlaps.add(overlap);
					curveschecked1[i]=true;
					curveschecked2[j]=true;
					/*
					mergedcs.add(mergeCurves(curveset1.get(i),curveset2.get(j), overlap));
					curveschecked1[i]=true;
					curveschecked2[j]=true;*/
					
				}							
			}//nCurves2 cycle end
	
		
		}//nCurves1 cycle end
		
		//build chains of connectivity
		extractChains(nCurves);
		if(bLinkDebug)
		{
			addOverlapsToOverlay(curveset1.nFrame);
		}
		
		for (i=0;i<allchains.size();i++)
		{
			mergedcs.addCurveSet(mergeChain(curveset1,curveset2, allchains.get(i)));
		}
		
		//let's add curves that where not merged
		for (i=0;i<nCurves[0];i++)
		{
			if(!curveschecked1[i])
				mergedcs.add(curveset1.get(i).duplicate());
		}
		for (i=0;i<nCurves[1];i++)
		{
			if(!curveschecked2[i])
				mergedcs.add(curveset2.get(i).duplicate());
		}
		
		return mergedcs;
	}

	/** given two curves and overlap object between them returns
	 * new curve that is a merge between them**/
	public CurveSet mergeChain(CurveSet curveset1, CurveSet curveset2, Chain chain)
	{
		CurveSet simplemerge = new CurveSet ();
		CurveSet finalcurveset = new CurveSet ();
		Curve finalcurve = new Curve ();
		//Curve intermcurve = new Curve ();
		int i;
		CurvesOverlap overlapx;
		
		//merge all pairs
		for (i=0;i<chain.size();i++)
		{
			simplemerge.add(mergeCurves(curveset1.get(chain.get(i)[0]), curveset2.get(chain.get(i)[1]), alloverlaps.get(chain.get(i)[2]-1)));
		}
		
		//first pair
		finalcurve =simplemerge.get(0);
		//finalcurve = mergeCurves(curveset1.get(chain.get(0)[0]), curveset2.get(chain.get(0)[1]), alloverlaps.get(chain.get(0)[2]-1));
		for (i=1;i<chain.size();i++)
		{
			overlapx=calculateOverlap(finalcurve,simplemerge.get(i));
			if(overlapx.bDirectConflict)
			{
				//IJ.log("Conflict points "+Integer.toString(overlapx.nDirectConfCount) );
				//IJ.log("Overlap size "+Integer.toString(overlapx.nOverlapSize) );
				IJ.log("--------MERGE STAGE BEG------");
				IJ.log("Overlap type "+Integer.toString(overlapx.nOvelapType) );
				String indices = "indx ";
				for (int h=0;h<overlapx.dirChange.size();h++)
					indices = indices +Integer.toString(overlapx.dirChange.get(h))+" ";
				IJ.log("Indices "+indices + Integer.toString(overlapx.nOverlapSize));
				IJ.log("--------MERGE STAGE END------");
			}
			if (overlapx.nOvelapType>=2)
			{
				finalcurveset.add(finalcurve);
				finalcurve = simplemerge.get(i);
				
			}
			else
			{
				finalcurve = mergeCurves(finalcurve,simplemerge.get(i), overlapx);
			}
			
			/*
			overlapx = calculateOverlap(finalcurve,curveset1.get(chain.get(i)[0]));
			finalcurve = mergeCurves(finalcurve,curveset1.get(chain.get(i)[0]), overlapx);
			overlapx = calculateOverlap(finalcurve,curveset2.get(chain.get(i)[1]));
			finalcurve = mergeCurves(finalcurve,curveset2.get(chain.get(i)[1]), overlapx);
*/
		/*	
			intermcurve = mergeCurves(curveset1.get(chain.get(i)[0]), curveset2.get(chain.get(i)[1]), alloverlaps.get(chain.get(i)[2]-1));
			overlapx = calculateOverlap(finalcurve,intermcurve);
			if (overlapx.nOvelapType>=2)
			{
				IJ.log("frames "+Integer.toString(curveset1.nFrame)+" & "+Integer.toString(curveset2.nFrame));
				IJ.log("curve pair "+Integer.toString(chain.get(i)[0]+1)+"   "+Integer.toString(chain.get(i)[1]+1));
				IJ.log("overlap type "+Integer.toString(overlapx.nOvelapType));
			}
			finalcurve = mergeCurves(finalcurve,intermcurve, overlapx);
			*/
		}
		finalcurveset.add(finalcurve);
		return finalcurveset;
	
	}
	
	/** function extracting chains (of connected curves)
	 * from connectivity matrix using Breadth-first search (BFS) algorithm **/
	public void extractChains(int [] nCurves)
	{
		int i,j,k, nPair;
		int [] pair;
		int nC1,nC2;
		
		Chain chain;
		int nChainCurrLength;
		allchains = new ArrayList<Chain>();
		
		for (i=0; i<nCurves[0]; i++)
			for (j=0; j<nCurves[1]; j++)
			{
				if(dConnectivityMatrix[i][j]>0)
				{
					chain = new Chain(i, j, dConnectivityMatrix[i][j]);
					dConnectivityMatrix[i][j]=0;
					while (!chain.isProcessed())
					{
						nChainCurrLength = chain.size();
						for (nPair=0;nPair<nChainCurrLength;nPair++)
						{
							pair=chain.get(nPair);
							if(pair[3]>0)
							{
								//going along columns
								nC1=pair[0];
								for(k=0;k<nCurves[1];k++)
								{
									if(dConnectivityMatrix[nC1][k]>0)
									{
										chain.addPair(nC1,k,dConnectivityMatrix[nC1][k]);
										dConnectivityMatrix[nC1][k]=0;										
									}
								}
								//going along rows
								nC2=pair[1];
								for(k=0;k<nCurves[0];k++)
								{
									if(dConnectivityMatrix[k][nC2]>0)
									{
										chain.addPair(k,nC2,dConnectivityMatrix[k][nC2]);
										dConnectivityMatrix[k][nC2]=0;										
									}
								}
								
								//mark this pair as processed
								pair[3]=0;
							}//if(pair[3]>0)
							
						}
					}//end while (!chain.isProcessed())
					
					allchains.add(chain);
					
				}//end of if(dConnectivityMatrix[i][j]>0)
				
			}//end of for (j=0; j<nCurves[1]; j++) 			
	}// end of for (i=0; i<nCurves[0]; i++)
	
	
	/** given two curves and overlap object between them returns
	 * new curve that is a merge between them**/
	public Curve mergeCurves(Curve curve1, Curve curve2, CurvesOverlap overlap)
	{
		Curve finalcurve = new Curve ();
		
		//for simplicity make array of curves
		Curve [] curves = new Curve [2];
		int [] nPoints = new int [2];
		int nIndBeg,nIndEnd;
		int i;
		float [] weights;
		
		curves[0]= curve1;
		curves[1]= curve2;
		
		nPoints[0]=curve1.size();
		nPoints[1]=curve2.size();
		
		//let's determine which curve we use as a beginning of merged curve
		// and which curve for the end (stored in nIndBeg, nIndEnd
		
		//same numbering direction
		if(overlap.nOverlapDirection==0)
		{
			if (overlap.indminmax[0][0]>overlap.indminmax[1][0])
				nIndBeg=0;
			else
				nIndBeg=1;
			//adding beginning from longer curve
/*			for (i=0;i<overlap.indminmax[nIndBeg][0];i++)
			{
				finalcurve.add(curves[nIndBeg].get(i).duplicate());
			}
	*/		
		}
		//numbering is opposite
		else
		{
			//beginning of curve 1 is longer
			if((overlap.indminmax[0][0]-1)>=nPoints[1]-overlap.indminmax[1][1])
			{
				nIndBeg=0;
				/*for (i=0;i<overlap.indminmax[0][0];i++)
				{
					finalcurve.add(curves[0].get(i).duplicate());
				}*/
				
			}
			//end of curve 2 is longer
			else
			{
				nIndBeg=1;
				//bOverlapReverse = true;
				/*for (i=nPoints[1]-1;i>overlap.indminmax[1][1];i--)
				{
					finalcurve.add(curves[1].get(i).duplicate());
				}
				*/
			}
			
		}//end of else, case of overlap.nOverlapDirection==1
		

		//adding the last part of curves
		//same numbering direction
		if(overlap.nOverlapDirection==0)
		{
			if((nPoints[0]-overlap.indminmax[0][1])>(nPoints[1]-overlap.indminmax[1][1]))
			{
				nIndEnd=0;
			}
			else
			{
				nIndEnd=1;
			}
	/*		
			for(i=overlap.indminmax[nIndEnd][1]+1;i<nPoints[nIndEnd];i++)
			{
				finalcurve.add(curves[nIndEnd].get(i).duplicate());
			}*/
		}
		else
		{
			//end of curve 1 is longer
			if((nPoints[0]-overlap.indminmax[0][1])>=overlap.indminmax[1][0]-1)
			{
				nIndEnd=0;
				/*for(i=overlap.indminmax[0][1]+1;i<nPoints[0];i++)
				{
					finalcurve.add(curves[0].get(i).duplicate());
				}*/
				
			}
			//beginning of curve 2 is longer
			else
			{
				nIndEnd=1;
				
				/*
				for(i=overlap.indminmax[1][0]-1;i>-1;i--)
				{
					finalcurve.add(curves[1].get(i).duplicate());
				}*/
				
			}
			
			
		}//end of else, case of overlap.nOverlapDirection==1
		
		weights = new float[overlap.nOverlapSize];
		//let's recalculate overlay
		//calculating weights
		float fOverlapSizeMOne=(float) ((float)overlap.nOverlapSize-1.0);
		for (i=0;i<overlap.nOverlapSize;i++)
		{
			if (nIndBeg==nIndEnd)
				if(nIndBeg==0)
					weights[i]=(float) 1.0;
				else
					weights[i]=(float) 0.0;
			else
			{
				if (nIndBeg==0)
					weights[i]=(fOverlapSizeMOne-(float) i)/fOverlapSizeMOne;
				else
					weights[i]=((float) i)/fOverlapSizeMOne;
			}
		}
		//calculate overlap position based on weighted average of points
		int [] nOvelapInd;// =new int [2];
		Point point;
		overlap.averaged_points = new ArrayList<Point>();
		for (i=0;i<overlap.nOverlapSize;i++)
		{
			nOvelapInd = overlap.indexes.get(i);
			point = curve1.get(nOvelapInd[0]).averageWeightedPoint(curve2.get(nOvelapInd[1]), weights[i]);
			//overlap.averaged_points.set(i,point.duplicate());
			overlap.averaged_points.add(point.clone());
			//overlap.averaged_points.get(i) = curve1.get(nOvelapInd[0]).averageWeightedPoint();
		}
		
		//now let's add curves
		
		///////////////////beginning
		if(overlap.nOverlapDirection==0)
		{
			//adding beginning from longer curve
			for (i=0;i<overlap.indminmax[nIndBeg][0];i++)
			{
				finalcurve.add(curves[nIndBeg].get(i).clone());
			}
		}
		else
		{
			if(nIndBeg==0)
				for (i=0;i<overlap.indminmax[0][0];i++)
				{
					finalcurve.add(curves[0].get(i).clone());
				}
			else
				for (i=nPoints[1]-1;i>overlap.indminmax[1][1];i--)
				{
					finalcurve.add(curves[1].get(i).clone());
				}
			
		}
		
		//////////////////////overlap
		for(i=0;i<overlap.nOverlapSize;i++)
		{
			finalcurve.add(overlap.averaged_points.get(i).clone());
		}
		
		//adding the last part of curves
				//same numbering direction
		if(overlap.nOverlapDirection==0)
		{
			for(i=overlap.indminmax[nIndEnd][1]+1;i<nPoints[nIndEnd];i++)
			{
				finalcurve.add(curves[nIndEnd].get(i).clone());
			}	
		}
		else
		{
			if(nIndEnd==0)
				for(i=overlap.indminmax[0][1]+1;i<nPoints[0];i++)
				{
					finalcurve.add(curves[0].get(i).clone());
				}
			else
				for(i=overlap.indminmax[1][0]-1;i>-1;i--)
				{
					finalcurve.add(curves[1].get(i).clone());
				}
				
		}
		
		return finalcurve;
		
	}
	
	
	
	/** function calculating overlap between two curves **/
	public CurvesOverlap calculateOverlap(Curve curve1, Curve curve2)
	{
		CurvesOverlap overlap = new CurvesOverlap();
		/** number of points in each curve**/
		int [] nPoints = new int[2];
		int i,j;
		double distance, dMinDistance;
		double dAnglDiff;
		int nMinIndex=0;
		Point point1, point2;
		int [] indexes;
		double dWPointDistance;
		
		
		
		nPoints[0] = curve1.size();
		nPoints[1] = curve2.size();
		
		for (i=0;i<nPoints[0];i++)
		{
			point1=curve1.get(i);
			dMinDistance = Double.MAX_VALUE;
			for (j=0;j<nPoints[1];j++)
			{
				point2=curve2.get(j);
				// euclidian distance
				distance = Math.sqrt(Math.pow(point1.coords[0]-point2.coords[0], 2)+Math.pow(point1.coords[1]-point2.coords[1], 2));

				if (distance<dMinDistance)
				{
					/*
					//angle distance
					point1.angle=(float) (point1.angle%Math.PI);
					point2.angle=(float) (point2.angle%Math.PI);
					if(point1.angle>0.5*Math.PI)
						point1.angle=(float) (point1.angle-Math.PI);
					if(point2.angle>0.5*Math.PI)
						point2.angle=(float) (point2.angle-Math.PI);
					dAnglDiff =Math.abs(point1.angle-point2.angle);
					*/
					dAnglDiff=Math.cos(point1.angle-point2.angle);
					dAnglDiff=Math.acos(dAnglDiff);
					dAnglDiff=Math.min(dAnglDiff, Math.PI-dAnglDiff);
					if(dAnglDiff<dAngleDistance)
					{
						dMinDistance=distance;
						nMinIndex=j;
					}
				}
			}//end of nPoints1 cycle
			
			//flexible ends
			/*
			dWPointDistance=dPointDistance;
			if (i<5 || i>(nPoints[0]-5))
			{
				dWPointDistance*=2.0;
				
			}
			if (j<5 || i>(nPoints[1]-5))
			{
				dWPointDistance*=2.0;
				
			}
			*/
			//found overlapping points
			
			//if(dMinDistance<dWPointDistance)
			if(dMinDistance<dPointDistance)
			{
				point2=curve2.get(nMinIndex);
				indexes=new int[2];
				indexes[0]=i;
				indexes[1]=nMinIndex;
				overlap.indexes.add(indexes);
				overlap.averaged_points.add(point1.averagePoint(point2));				
			}
		}//end of nPoints2 cycle
		
		overlap.nOverlapSize=overlap.indexes.size();
		
		overlap.analyzeOverlap();
		
		if (overlap.nOverlapSize>=nOverLapThreshold)
		{

			//overlap analysis
			overlap.calculateOverlapDirection();
			overlap.calculateMinMaxIndexes();
			
			//figure out overlap type
			overlap.overlapFraction1=(double)overlap.nOverlapSize/(double)nPoints[0];
			overlap.overlapFraction2=(double)overlap.nOverlapSize/(double)nPoints[1];
			
			//FULL overlap
			if(Math.max(overlap.overlapFraction1, overlap.overlapFraction2)>=dFullOverlapF)
			{
				overlap.nOvelapType= CurvesOverlap.FULL;
				//IJ.log("Gap " + Integer.toString(nStepDiff) +"L "+Integer.toString(overlap.nOverlapSize) + " T"+Integer.toString(overlap.nOvelapType));
				//double ratio =(double)nStepDiff/(double)overlap.nOverlapSize; 
				//IJ.log("" + Integer.toString(nStepDiff) +"/"+Integer.toString(overlap.nOverlapSize) +'='+Double.toString(ratio) +" ("+Integer.toString(overlap.nOvelapType)+")");
				//big gap in overlap
				if (overlap.dMaxJumpFraction>0.5)
				//if (ratio>0.5 && !overlap.bDirectConflict)
				{
					overlap.nOvelapType= CurvesOverlap.NONE;
				}
				return overlap;
			}
			//smaller overlap
			//let's see if overlap is close to ends of curves
			int [] nCurvesEndInd=new int[2];
			nCurvesEndInd[0]=-1;
			nCurvesEndInd[1]=-1;
			for (j=0;j<2;j++)
			{
				//minimum index
				if((double)overlap.indminmax[j][0]/(double)nPoints[j]<=dExtensionOverlapF)
				{
					nCurvesEndInd[j]=0;
				}
				//maximum index
				if((double)overlap.indminmax[j][1]/(double)nPoints[j]>=(1.0-dExtensionOverlapF))
				{
					nCurvesEndInd[j]=1;
				}
			}
			//crossing overlap, skip
			if(nCurvesEndInd[0]<0 || nCurvesEndInd[1]<0)
			{
				overlap.nOvelapType=CurvesOverlap.CROSS;
				return overlap;
			}
			
			//see if it is branching or extension overlap
			if(overlap.nOverlapDirection==0 && nCurvesEndInd[0]==nCurvesEndInd[1])
			{
				overlap.nOvelapType=CurvesOverlap.BRANCHING;
				return overlap;
				
			}
			if(overlap.nOverlapDirection==1 && nCurvesEndInd[0]!=nCurvesEndInd[1])
			{
				overlap.nOvelapType=CurvesOverlap.BRANCHING;
				return overlap;
			}
			else
			{
				//ok, it is extension overlap
				overlap.nOvelapType=CurvesOverlap.EXTENSION;
				//IJ.log("Gap " + Integer.toString(nStepDiff) +"L "+Integer.toString(overlap.nOverlapSize) + " T"+Integer.toString(overlap.nOvelapType));
				//double ratio =(double)nStepDiff/(double)overlap.nOverlapSize; 
				//IJ.log("" + Integer.toString(nStepDiff) +"/"+Integer.toString(overlap.nOverlapSize) +'='+Double.toString(ratio) +" ("+Integer.toString(overlap.nOvelapType)+")");
				//big gap in overlap
				//if (ratio>0.5 && !overlap.bDirectConflict)
				
				//see if it is consistent (no big point indices jumps inside)
				if (overlap.dMaxJumpFraction>0.5)
				{
					overlap.nOvelapType= CurvesOverlap.NONE;
				}
				return overlap;				
			}

			
		}
		//no overlap
		else
		{
			overlap.nOvelapType= CurvesOverlap.NONE;
		}
		
		
		
		return overlap;
		
	}
	
	public void addOverlapsToOverlay(int nFrame) 
	{
		int i,j,k;
		Chain chain;
		CurvesOverlap overlap;
		Color[] colors;
		PolygonRoi polyline_p;
		int pointsN;
		
		colors = new Color[8];		  
		colors[0] = new Color(0,0,255,126);//Color.BLUE;
		colors[1] = new Color(0,255,255,126);//Color.CYAN;
		colors[2] = new Color(0,255,0,126);//Color.GREEN;
		colors[3] = new Color(255,0,255,126);//Color.MAGENTA;
		colors[4] = new Color(255,165,0,126);//Color.ORANGE;
		colors[5] = new Color(255,182,193,126);//Color.PINK;	255-182-193
		colors[6] = new Color(255,0,0,126);//Color.RED;
		colors[7] = new Color(255,255,0,126);//Color.YELLOW;
		
		
		for (i=0;i<allchains.size();i++)
		{
			chain = allchains.get(i);
			for (j=0;j<chain.size();j++)
			{
				overlap = alloverlaps.get(chain.get(j)[2]-1);
				pointsN= overlap.averaged_points.size();
				float[] px = new float[pointsN];
				float[] py = new float[pointsN];
				for (k=0;k<pointsN;k++)
				{
					 px[k]=overlap.averaged_points.get(k).coords[0];
					 py[k]=overlap.averaged_points.get(k).coords[1];
				}
				 polyline_p = new PolygonRoi(px, py, Roi.POLYLINE);
				 polyline_p.setPosition(nFrame);
				 //polyline_p.setStrokeColor(new Color(0,0,255,126));
				 polyline_p.setStrokeColor(colors[i%8]);
				 polyline_p.setStrokeWidth(dPointDistance);
				 //polyline_p.setFillColor(new Color(0,0,255,126));
				 image_overlay.add(polyline_p);
			}
		}
		
	}

	
	/**method constructs intermediate concatenating z-projects (SUMS) of stack,
	 * i.e. 1, 1-2, 1-2-3, etc and stores it as new imp class variable**/
	public void constructZprojects()
	{
		ZProjector zproj = new ZProjector(imp);
		ImageStack imstack = new ImageStack(imp.getWidth(), imp.getHeight());
		ImageProcessor ip;
		int nSlice;
		String imTitle;
		ContrastEnhancer contrastEnh = new ContrastEnhancer();
		
		imp.setSliceWithoutUpdate(1);
		ip=imp.getProcessor();
		//imstack.addSlice((FloatProcessor) ip.duplicate().convertToFloat());
		
		imstack.addSlice(ip.duplicate().convertToByte(true));
		zproj.setMethod(ZProjector.SUM_METHOD);
		zproj.setStartSlice(1);
		for (nSlice=2; nSlice<=imp.getStackSize(); nSlice++)
		{
			zproj.setStopSlice(nSlice);
			zproj.doProjection(false);
			ip = zproj.getProjection().getProcessor().duplicate();
			//to make sure intensity in the same range
			//ip.multiply(1.0/(double)nSlice);
			contrastEnh.stretchHistogram(ip, 0.35);
			//IJ.run(imp, "Enhance Contrast", "saturated=0.35");
			
			//imstack.addSlice(ip.convertToByte(true));
			ip=ip.convertToByte(true);
			imstack.addSlice(ip);
		}
		imTitle = imp.getTitle();
		imp=new ImagePlus();
		imTitle = imTitle+"_accumulated_proj";
		imp.setStack(imTitle, imstack);
		
		
	}
	
	/** Dialog with linking parameters **/
	public boolean showLinkingDialog()
	{
		
		GenericDialog linkingD = new GenericDialog("Linking parameters");
		
		linkingD.addNumericField("Max distance between points in overlap:", Prefs.get("CurveTrace.linkPointDistance", 4), 2, 4,"pixels");
		linkingD.addNumericField("Max normale angle difference between points in overlap:", Prefs.get("CurveTrace.linkAngleDistance", 4), 2, 4,"degrees");
		linkingD.addNumericField("Minimum number of points in overlap:", Prefs.get("CurveTrace.linkOverLapThreshold", 4), 0, 3,"points");
		linkingD.addNumericField("Fraction of overlap to full merge curves (0-1):", Prefs.get("CurveTrace.linkFullOverlapF", 0.8), 2, 4,"fraction");
		linkingD.addNumericField("Overlap distance to ends for extention overlap (0-1):", Prefs.get("CurveTrace.linkExtensionOverlapF", 0.05), 2, 4,"fraction");		
		linkingD.addCheckbox("Show intermediate linking information?", Prefs.get("CurveTrace.linkDebug", false));
		linkingD.addCheckbox("Generate cumulative projection?", Prefs.get("CurveTrace.linkCumulProj", false));
		linkingD.setResizable(false);
		linkingD.showDialog();	
		if (linkingD.wasCanceled())
            return false;
		
		dPointDistance = linkingD.getNextNumber();
		Prefs.set("CurveTrace.linkPointDistance", dPointDistance);
		dAngleDistance = linkingD.getNextNumber();
		Prefs.set("CurveTrace.linkAngleDistance", dAngleDistance);
		//convert to radians
		dAngleDistance=Math.PI*dAngleDistance/180.0;
		nOverLapThreshold = (int)linkingD.getNextNumber();
		Prefs.set("CurveTrace.linkOverLapThreshold", nOverLapThreshold);

		dFullOverlapF = linkingD.getNextNumber();
		Prefs.set("CurveTrace.linkFullOverlapF", dFullOverlapF);
		dExtensionOverlapF = linkingD.getNextNumber();
		Prefs.set("CurveTrace.linkExtensionOverlapF", dExtensionOverlapF);
		bLinkDebug = linkingD.getNextBoolean();
		Prefs.set("CurveTrace.linkDebug", bLinkDebug);
		bCumulProj= linkingD.getNextBoolean();
		Prefs.set("CurveTrace.linkCumulProj", bCumulProj);
		
		return true;
	}

}
