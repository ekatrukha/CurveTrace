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


import ij.ImageListener;
import ij.ImagePlus;

import ij.gui.GenericDialog;

import ij.measure.Measurements;
import ij.plugin.PlugIn;



import java.awt.Scrollbar;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

public class ThresholdSal extends GenericDialog implements PlugIn, Measurements,
Runnable, ActionListener, AdjustmentListener, ItemListener, ImageListener {
	

	
	public ThresholdSal() {
					
		super("Threshold saliency map");
	
	
	}
	@Override
	public void imageOpened(ImagePlus imp) {
		// TODO Auto-generated method stub
	
	}

	@Override
	public void imageClosed(ImagePlus imp) {
		// TODO Auto-generated method stub

	}

	@Override
	public void imageUpdated(ImagePlus imp) {
		// TODO Auto-generated method stub

	}

	@Override
	public void itemStateChanged(ItemEvent arg0) {
		// TODO Auto-generated method stub

	}

	@Override
	public void adjustmentValueChanged(AdjustmentEvent e) {
		
		int uv, lv;
		Scrollbar sbupper = (Scrollbar)slider.elementAt(0);
		Scrollbar sblower = (Scrollbar)slider.elementAt(1);
		Scrollbar updateanother = null;
		boolean need_adjustment = false;
		uv =sbupper.getValue();
		lv =sblower.getValue();
		Object source = e.getSource();
		
		//one slider should always be behind
		
		if(uv<lv)
		{ 
			need_adjustment = true;
			for (int i=0; i<slider.size(); i++) 
			{
				if (source==slider.elementAt(i)) 
				{				
						if(i==0)
						{
							sblower.setValue(uv);
							updateanother = sblower;
						}
						else
						{
							sbupper.setValue(lv);
							updateanother = sbupper;
						}
				}
			}
			
		}
		
		super.adjustmentValueChanged(e);
		
		if (need_adjustment)
		{
			e.setSource(updateanother);
			super.adjustmentValueChanged(e);
		}

		
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		super.actionPerformed(e);
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void run(String arg) {
		// TODO Auto-generated method stub
		
	}


}
