#!/usr/bin/python
import argparse
import importlib
import ROOT
from ROOT import TTree, TBranch
ROOT.gSystem.Load("/nfs/slac/g/ldmx/users/smidd/ldmx-sw/install/lib/libEvent.so")
import os
import math
import sys
#import numpy as np
from array import array
from optparse import OptionParser
sys.path.insert(0, '../')

class WabAnalysis:
    """Configuration class for WabAnalysis

    Attributes
    ----------
    Document your parameters here!
    

    def __init__(self, name) :
        super().__init__(name,'ldmx::WabAnalysis')

        #change MODULE to the module you are in
        from LDMX.MODULE import include
        include.library()

        #Put parameters and their defaults here, for example
        self.my_parameter = 5 #accessed in C++ with getParameter<int>("my_parameter")
      """  
    def __init__(self, fn, ofn, tag):

        self.fin = ROOT.TFile(fn);
        self.tin = self.fin.Get("LDMX_Events")
        self.tag = int(tag);

        self.fn_out = ofn;
        self.fout = ROOT.TFile("hist_"+self.fn_out,"RECREATE");

        self.fn_sklim = ofn;
        self.fsklim = ROOT.TFile("sklim_"+self.fn_sklim,"RECREATE");
        self.tsklim = self.tin.CloneTree(0);

        self.evHeader = ROOT.ldmx.EventHeader()
        self.simParticles = ROOT.std.map(1,'ldmx::SimParticle');
        self.hcalHits = ROOT.std.vector('ldmx::HcalHit')();
        self.hcalSimHits = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.tin.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader ));
        self.tin.SetBranchAddress("HcalSimHits_wab",  ROOT.AddressOf( self.hcalSimHits ));
        self.tin.SetBranchAddress("HcalRecHits_wab",  ROOT.AddressOf( self.hcalHits ));
		# histograms
        self.histograms = {};
       
        #self.tout = ROOT.TTree("tree2", "tree1");
        #self.x = [];
        #self.y = [];
        #self.branchx = self.tout.Branch("HcalRecHitsX",  self.x, "self.x/F");
        self.histograms["hit_x"]= ROOT.TH1F("rechit_x","; rechit x (mm); N", 400, -2000, 2000);
        self.histograms["hit_y"] = ROOT.TH1F("rechit_y","; rechit y (mm); N", 400, -2000, 2000);
        self.histograms["hit_z"] = ROOT.TH1F("rechit_z","; rechit z (mm); N", 400, 0, 4000);
        self.histograms["total_event_energy_rec"]= ROOT.TH1F("Total Event E RecHits","; Total Event Energy [GeV]; N", 400, 0,6);
        self.histograms["rechit_edep"] = ROOT.TH1F(" E Dep Rec Hits","; Rec Hit Energy Dep [GeV]; N", 400, 0,6);
        self.histograms["simhit_x"]= ROOT.TH1F("hit_x","; hit x (mm); N", 400, -800, 800);
        self.histograms["simhit_y"] = ROOT.TH1F("hit_y","; hit y (mm); N", 400, -800,800);
        self.histograms["simhit_z"] = ROOT.TH1F("hit_z","; hit z (mm); N", 400, 0, 800);
        self.histograms["simhit_edep"] = ROOT.TH1F(" E DepSimHits","; Sim Hit Energy Dep [GeV]; N", 400, 0,6);
        self.histograms["total_event_energy_sim"]= ROOT.TH1F("Total Event E SimHits","; Total Event Energy [GeV]; N", 400, 0,6);
        self.histograms["ehit_x"]= ROOT.TH1F("ehit_x","; ehit x (mm); N", 400, -2000, 2000);
        self.histograms["ehit_y"] = ROOT.TH1F("ehit_y","; ehit y (mm); N", 400, -2000, 2000);
        self.histograms["ehit_z"] = ROOT.TH1F("ehit_z","; ehit z (mm); N", 400, 0, 4000);
        self.histograms["esimhit_x"]= ROOT.TH1F("ehit_x","; ehit x (mm); N", 400, -800, 800);
        self.histograms["esimhit_y"] = ROOT.TH1F("ehit_y","; ehit y (mm); N", 400, -800,800);
        self.histograms["esimhit_z"] = ROOT.TH1F("ehit_z","; ehit z (mm); N", 400, 0, 800);
        self.histograms["simhit_theta"] = ROOT.TH1F("hit_theta","; hit theta (rad); N", 400, 0, math.pi);
        self.histograms["hit_strip"]   = ROOT.TH1F("hit_strip","; hit strip; N", 100, 0, 100);		
        self.histograms["hit_layer"]   = ROOT.TH1F("hit_layer","; hit_layer; N", 100, 0, 100);		
        self.histograms["hit_section"] = ROOT.TH1F("hit_section","; hit_section; N", 10, 0, 10);		

        self.histograms["hit_PEs"] = ROOT.TH1F("hit_PEs","; hit PEs; N", 80, 0, 4000);		    
        self.histograms["hit_MinPEs"] = ROOT.TH1F("hit_MinPEs","; hit minPEs; N", 50, 0, 2000);		    
        self.histograms["event_MaxPE"] = ROOT.TH1F("event_MaxPE","; event maxPEs; N", 50, 0, 5000);		    

        self.histograms["hit_PEs_back"] = ROOT.TH1F("hit_PEs_back","; hit PEs; N", 10, 0, 10);		    
        self.histograms["hit_MinPEs_back"] = ROOT.TH1F("hit_MinPEs_back","; hit minPEs; N", 10, 0, 10);		    
        self.histograms["event_MaxPE_back"] = ROOT.TH1F("event_MaxPE_back","; event maxPEs; N", 10, 0, 10);		    
        self.histograms["hit_PEs_side"] = ROOT.TH1F("hit_PEs_side","; hit PEs; N", 10, 0, 10);		    
        self.histograms["hit_MinPEs_side"] = ROOT.TH1F("hit_MinPEs_side","; hit minPEs; N", 10, 0, 10);		    
        self.histograms["event_MaxPE_side"] = ROOT.TH1F("event_MaxPE_side","; event maxPEs; N", 10, 0, 10);	

        self.histograms["hit_strip_v_z"]   = ROOT.TH2F("hit_strip_v_z","; hit strip; hit z (mm)", 100, 0, 100, 400, 0, 4000);		
        self.histograms["hit_strip_v_layer"]   = ROOT.TH2F("hit_strip_v_layer","; hit strip; hit layer", 100, 0, 100, 100, 0, 100);		

        self.histograms["hit_PEs_v_hit_MinPEs"]      = ROOT.TH2F("hit_PEs_v_hit_MinPEs","; hit PEs; hit MinPEs", 10, 0, 10, 10, 0, 10);		    		
        self.histograms["hit_PEs_v_hit_MinPEs_back"] = ROOT.TH2F("hit_PEs_v_hit_MinPEs_back","; hit PEs; hit MinPEs", 10, 0, 10, 10, 0, 10);		    
        self.histograms["hit_PEs_v_hit_MinPEs_side"] = ROOT.TH2F("hit_PEs_v_hit_MinPEs_side","; hit PEs; hit MinPEs", 10, 0, 10, 10, 0, 10);		    


        ## gen information
        # inclusive
        # not vetoed by hcal

        self.loop();
        self.writeOutHistos();

    def writeOutHistos(self):

        self.fout.cd();
        for key in self.histograms: 
            print(key);
            self.histograms[key].Write();
        
        #self.tout.Write();
        self.fout.Close();

        self.fsklim.cd();
        self.tsklim.Write();
        self.fsklim.Close();

        ####################################################
        ### looping 
        ####################################################
    def loop(self):

        nent = self.tin.GetEntriesFast();
        # print "nent = ", nent
        #np.empty((1), dtype="float32");
        #self.tree.Branch("x", self.x, "x/F");
        #self.tree.Branch("y", self.y, "y/F");
        for i in range(nent):
            if(nent/100 > 0 and i % (1 * nent/100) == 0):
	            sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done")
	            sys.stdout.flush()

            self.tin.GetEntry(i);

            # hcal veto info
            maxPEs = 0;
            maxPEs_back = 0;
            maxPEs_side = 0;
            energy_sum_sim = 0;
            energy_sum_rec = 0;
            for ih,hit in enumerate(self.hcalSimHits):
                [x, y, z] = hit.getPosition();			
                self.histograms["simhit_x"].Fill(x);
                self.histograms["simhit_y"].Fill(y);
                self.histograms["simhit_z"].Fill(z);
                self.histograms["simhit_edep"].Fill(hit.getEdep());
                energy_sum_sim += hit.getEdep();
                self.histograms["total_event_energy_sim"].Fill(energy_sum_sim);
                r = math.sqrt(x*x+y*y+z*z);
                theta = math.acos(z/r);
                self.histograms["simhit_theta"].Fill(theta);
            for ih,hit in enumerate(self.hcalHits):		
                if (hit.isNoise()==0 and hit.getX()!=0 and hit.getY()!=0 and hit.getZ()!=0):
                    print("hit PEs = ", hit.getPE(),", hit MinPEs = ",hit.getMinPE(),",\t id = ", hit.getID(), ", layer = ", hit.getLayer(), ", section = ", hit.getSection(), ", strip = ", hit.getStrip(), ", hit X = ", hit.getX(),", hit Y = ", hit.getY(),", hit Z = ", hit.getZ());
                    self.histograms["hit_x"].Fill(hit.getX());
                    self.histograms["hit_y"].Fill(hit.getY());
                    self.histograms["hit_z"].Fill(hit.getZ());
                    self.histograms["ehit_x"].Fill(hit.getX());
                    self.histograms["ehit_y"].Fill(hit.getY());
                    self.histograms["ehit_z"].Fill(hit.getZ());
                    self.histograms["hit_strip"].Fill(hit.getStrip());
                    self.histograms["hit_layer"].Fill(hit.getLayer());
                    self.histograms["hit_section"].Fill(hit.getSection());
                    self.histograms["hit_PEs"].Fill(hit.getPE());
                    self.histograms["hit_MinPEs"].Fill(hit.getMinPE());
                    energy_sum_rec += hit.getEnergy();
                    self.histograms["rechit_edep"].Fill(hit.getEnergy())
                    self.histograms["hit_PEs_v_hit_MinPEs"].Fill(hit.getPE(),hit.getMinPE());
                    if hit.getSection() == 0:
                        self.histograms["hit_PEs_back"].Fill(hit.getPE());
                        self.histograms["hit_MinPEs_back"].Fill(hit.getMinPE());
                        self.histograms["hit_PEs_v_hit_MinPEs_back"].Fill(hit.getPE(),hit.getMinPE());
                    else:					
                        self.histograms["hit_PEs_side"].Fill(hit.getPE());
                        self.histograms["hit_MinPEs_side"].Fill(hit.getMinPE());
                        self.histograms["hit_PEs_v_hit_MinPEs_side"].Fill(hit.getPE(),hit.getMinPE());

                        self.histograms["hit_strip_v_z"].Fill(hit.getStrip(),hit.getZ());
                        self.histograms["hit_strip_v_layer"].Fill(hit.getStrip(),hit.getLayer());


                        if hit.getMinPE() > 0 and hit.getPE() > maxPEs: maxPEs = hit.getPE();
                        if hit.getMinPE() > 0 and hit.getPE() > maxPEs_back and hit.getSection() == 0: maxPEs_back = hit.getPE();
                    if hit.getMinPE() > 0 and hit.getPE() > maxPEs_side and hit.getSection() > 0:  maxPEs_side = hit.getPE();
			    
                    self.histograms["event_MaxPE"].Fill(maxPEs);
                    self.histograms["event_MaxPE_back"].Fill(maxPEs_back);
                    self.histograms["event_MaxPE_side"].Fill(maxPEs_side);
                    self.histograms["total_event_energy_sim"].Fill(energy_sum_sim);
                    self.histograms["total_event_energy_rec"].Fill(energy_sum_rec);
                    if i > 1000: break;
			    
                #print( "\n");
			
			
def main(options,args) : 
	
	#print "Hello!"
	sc = WabAnalysis(options.ifile,options.ofile,options.tag) ;

def makeCanvas(hists, tags, norm=False, logy=True):

	colors = [1,2,4,6,7]

	leg = ROOT.TLegend(0.6,0.7,0.85,0.9);
	leg.SetBorderSize(0);
	leg.SetFillStyle(0);
	leg.SetTextSize(0.04);
	leg.SetTextFont(42);
	tmax = -999;
	tmin = 9999;
	for i in range(len(hists)):
		leg.AddEntry(hists[i],tags[i],'l');
		if norm: hists[i].Scale(1./hists[i].Integral());
		if hists[i].GetMaximum() > tmax: tmax = hists[i].GetMaximum();
		if hists[i].GetMinimum() < tmin: tmin = hists[i].GetMinimum();

	c = ROOT.TCanvas("c","c",1000,800);
	hists[0].SetMinimum(1e-7);
	hists[0].SetMaximum(tmax*1.2);
	for i in range(len(hists)):
		if i == 0: hists[i].Draw();
		else: 
			hists[i].SetLineColor( colors[i] );
			hists[i].Draw('sames');

	leg.Draw();
	ROOT.gPad.SetLogy();
	c.SaveAs("plots/"+hists[0].GetName()+".pdf")
	c.SaveAs("plots/"+hists[0].GetName()+".png")


if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

	parser.add_option('-i','--ifile', dest='ifile', default = 'file.root',help='directory with data', metavar='idir')
	parser.add_option('-o','--ofile', dest='ofile', default = 'ofile.root',help='directory to write plots', metavar='odir')

	parser.add_option('--tag', dest='tag', default = '1',help='file tag', metavar='tag')

	(options, args) = parser.parse_args()

	
	ROOT.gStyle.SetPadTopMargin(0.10)
	ROOT.gStyle.SetPadLeftMargin(0.16)
	ROOT.gStyle.SetPadRightMargin(0.10)
	ROOT.gStyle.SetPalette(1)
	ROOT.gStyle.SetPaintTextFormat("1.1f")
	ROOT.gStyle.SetOptFit(0000)
	ROOT.gROOT.SetBatch()

	# Get the Event library 
	ROOT.gSystem.Load("/nfs/slac/g/ldmx/users/smidd/ldmx-sw/install/lib/libEvent.so");	

	main(options,args);			
			
			
			
			
			
			
			
