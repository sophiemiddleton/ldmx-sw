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

lumin = 0.22

#signal
proton_signal = 1
xsec_signal = 22203.
n_signal = 15129
event_weight_signal = xsec_signal/n_signal
factor_signal = event_weight_signal * proton_signal * lumin

#wabs
proton_wabs = 74.
xsec_wab = 1795.
n_wab = 10000.
lumin_wab = 0.22
event_weight_wab = xsec_wab/n_wab
factor_wab = event_weight_wab * proton_wabs * lumin
class WabSigCompare:
   
    def __init__(self, fn1, fn2, ofn, tag):

        #wab:
        self.fin1 = ROOT.TFile(fn1);
        self.tin1 = self.fin1.Get("LDMX_Events")
        
        #sig:
        self.fin2 = ROOT.TFile(fn2);
        self.tin2 = self.fin2.Get("LDMX_Events")
        
        self.tag = int(tag);

        self.fn_out = ofn;
        self.fout = ROOT.TFile("hist_"+self.fn_out,"RECREATE");

        
        #wab:
        self.evHeader1 = ROOT.ldmx.EventHeader()
        self.simParticles1 = ROOT.std.map(int,'ldmx::SimParticle')();
        self.hcalHits1 = ROOT.std.vector('ldmx::HcalHit')();
        self.hcalSimHits1 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.tin1.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader1 ));
        self.tin1.SetBranchAddress("HcalSimHits_wab",  ROOT.AddressOf( self.hcalSimHits1 ));
        self.tin1.SetBranchAddress("HcalRecHits_wab",  ROOT.AddressOf( self.hcalHits1 ));
        #sig:
        self.evHeader2 = ROOT.ldmx.EventHeader()
        self.simParticles2 = ROOT.std.map(int,'ldmx::SimParticle')();
        self.hcalHits2 = ROOT.std.vector('ldmx::HcalHit')();
        self.hcalSimHits2 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.tin2.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader2 ));
        self.tin2.SetBranchAddress("HcalSimHits_wab",  ROOT.AddressOf( self.hcalSimHits2 ));
        self.tin2.SetBranchAddress("HcalRecHits_wab",  ROOT.AddressOf( self.hcalHits2 ));

        self.histograms = {};

        #SimHits WAB:
        self.histograms["simhit_N_wab"]= ROOT.TH1F("HCAL simhit_HitsPerEvent_wab",";HCAL SimHits per Event; Weighted N", 400, 0,1000);
        self.histograms["simhit_x_wab"]= ROOT.TH1F("HCAL simhit_x_wab",";HCAL simhit x [mm]; Weighted N", 400, -800, 800);
        self.histograms["simhit_y_wab"] = ROOT.TH1F("HCAL simhit_y_wab",";HCAL simhit y [mm]; Weighted N", 400, -800,800);
        self.histograms["simhit_z_wab"] = ROOT.TH1F("HCAL simhit_z_wab",";HCAL simhit z [mm]; Weighted N", 400, 0, 2000);
        self.histograms["simhit_edep_wab"] = ROOT.TH1F("HCAL Energy Dep per SimHit_wab",";HCAL Sim Hit Energy Dep [MeV]; Weighted N", 400, 0,2);
        self.histograms["simhit_totE_wab"]= ROOT.TH1F("HCAL Total Event E SimHits_wab","; HCAL Total Event SimHit Energy [MeV]; Weighted N", 400, 0,400);
        self.histograms["simhit_theta_wab"] = ROOT.TH1F("HCAL simhit_theta_wab",";HCAL simhit Polar Angle [rad]; Weighted N", 400, 0, math.pi);
        
        # Sim Hits Signal:
        
        self.histograms["simhit_N_sig"]= ROOT.TH1F("HCAL simhit_HitsPerEvent_sig",";HCAL SimHits per Event; Weighted N", 400, 0,1000);
        self.histograms["simhit_x_sig"]= ROOT.TH1F("HCAL simhit_x_sig",";HCAL simhit x [mm]; Weighted N", 400, -800, 800);
        self.histograms["simhit_y_sig"] = ROOT.TH1F("HCAL simhit_y_sig",";HCAL simhit y [mm]; Weighted N", 400, -800,800);
        self.histograms["simhit_z_sig"] = ROOT.TH1F("HCAL simhit_z_sig",";HCAL simhit z [mm]; Weighted N", 400, 0, 2000);
        self.histograms["simhit_edep_sig"] = ROOT.TH1F("HCAL Energy Dep per SimHit_sig",";HCAL Sim Hit Energy Dep [MeV]; Weighted N", 400, 0,2);
        self.histograms["simhit_totE_sig"]= ROOT.TH1F("HCAL Total Event E SimHits_sig","; HCAL Total Event SimHit Energy [MeV]; Weighted N", 400, 0,400);
        self.histograms["simhit_theta_sig"] = ROOT.TH1F("HCAL simhit_theta_sig",";HCAL simhit Polar Angle [rad]; Weighted N", 400, 0, math.pi);
        
        #RecoHits WABS:
        self.histograms["rechit_N_wab"]= ROOT.TH1F("HCAL rechit_HitsPerEvent_wab",";HCAL RecHits per Event; Weighted N", 400, 0,100);
        self.histograms["rechit_x_wab"]= ROOT.TH1F("HCAL rechit_x_wab",";HCAL rechit x [mm]; Weighted N", 400, -1500, 1500);
        self.histograms["rechit_y_wab"] = ROOT.TH1F("HCAL rechit_y_wab",";HCAL rechit y [mm]; Weighted N", 400, -1500, 1500);
        self.histograms["rechit_z_wab"] = ROOT.TH1F("HCAL rechit_z_wab",";HCAL rechit z [mm]; Weighted N", 400, 0, 2000);
        self.histograms["rechit_totE_wab"]= ROOT.TH1F("HCAL Total Event E RecHits_wab","; HCAL Total Event RecHit Energy [MeV]; Weighted N", 400, 0,400);
        self.histograms["rechit_edep_wab"] = ROOT.TH1F("HCAL Energy Dep per Rec Hits_wab","; HCAL Rec Hit Energy Dep [MeV]; Weighted N", 400, 0,10);
        self.histograms["rechit_theta_wab"] = ROOT.TH1F(" rechit_theta_wab",";HCAL Rec Hit Polar Angle [rad]; Weighted N", 400, 0,math.pi);
        
        #RecoHits Signal:
        self.histograms["rechit_N_sig"]= ROOT.TH1F("HCAL rechit_HitsPerEvent_sig","; HCAL RecHits per Event; Weighted N", 400, 0,100);
        self.histograms["rechit_x_sig"]= ROOT.TH1F("HCAL rechit_x_sig","; HCAL rechit x [mm]; Weighted N", 400, -1500, 1500);
        self.histograms["rechit_y_sig"] = ROOT.TH1F("HCAL rechit_y_sig","; HCAL rechit y [mm]; Weighted N", 400, -1500, 1500);
        self.histograms["rechit_z_sig"] = ROOT.TH1F("HCAL rechit_z_sig","; HCAL rechit z [mm]; Weighted N", 400, 0, 2000);
        self.histograms["rechit_totE_sig"]= ROOT.TH1F("HCAL Total Event E RecHits_sig","; HCAL Total Event RecHit Energy [MeV]; Weighted N", 400, 0,400);
        self.histograms["rechit_edep_sig"] = ROOT.TH1F("HCAL Energy Dep per Rec Hits_sig","; HCAL Rec Hit Energy Dep [MeV]; Weighted N", 400, 0,10);
        self.histograms["rechit_theta_sig"] = ROOT.TH1F("HCAL rechit_theta_sig","; HCAL Rec Hit Polar Angle [rad]; Weighted N", 400, 0,math.pi);

        self.loop();
        self.writeOutHistos();

    def writeOutHistos(self):

        self.fout.cd();
        for key in self.histograms: 
            print(key);
            self.histograms[key].Write();
        
    #WABs:
    def loop(self):
        nentWAB = self.tin1.GetEntriesFast();
        for i in range(nentWAB):
           
            self.tin1.GetEntry(i);
            energy_sum_sim = 0;
            energy_sum_rec = 0;
            simN = 0;
            recN = 0;
            for ih,hit in enumerate(self.hcalSimHits1):
                [x, y, z] = hit.getPosition();	
                simN +=1;
                self.histograms["simhit_x_wab"].Fill(x);
                self.histograms["simhit_y_wab"].Fill(y);
                self.histograms["simhit_z_wab"].Fill(z);
                self.histograms["simhit_edep_wab"].Fill(hit.getEdep());
                energy_sum_sim += hit.getEdep();
                r = math.sqrt(x*x+y*y+z*z);
                thetasim = math.acos(z/r);
                self.histograms["simhit_theta_wab"].Fill(thetasim);
            for ih,hit in enumerate(self.hcalHits1):		
                if (hit.isNoise()==0 and hit.getX()!=0 and hit.getY()!=0 and hit.getZ()!=0):
                    recN +=1;
                    self.histograms["rechit_x_wab"].Fill(hit.getX());
                    self.histograms["rechit_y_wab"].Fill(hit.getY());
                    self.histograms["rechit_z_wab"].Fill(hit.getZ());
                    d = math.sqrt(hit.getX()*hit.getX()+hit.getY()*hit.getY()+hit.getZ()*hit.getZ());
                    thetarec = math.acos(hit.getZ()/d);
                    self.histograms["rechit_theta_wab"].Fill(thetarec);
                    self.histograms["rechit_edep_wab"].Fill(hit.getEnergy());
                    energy_sum_rec += hit.getEnergy();
            self.histograms["simhit_totE_wab"].Fill(energy_sum_sim);
            self.histograms["rechit_totE_wab"].Fill(energy_sum_rec);
            self.histograms["rechit_N_wab"].Fill(recN);
            self.histograms["simhit_N_wab"].Fill(simN);
        nentSIG = self.tin2.GetEntriesFast();
        for i in range(nentSIG):
           
            self.tin2.GetEntry(i);
            energy_sum_sim = 0;
            energy_sum_rec = 0;
            simN = 0;
            recN = 0;
            for ih,hit in enumerate(self.hcalSimHits2):
                [x, y, z] = hit.getPosition();	
                simN +=1;
                self.histograms["simhit_x_sig"].Fill(x);
                self.histograms["simhit_y_sig"].Fill(y);
                self.histograms["simhit_z_sig"].Fill(z);
                self.histograms["simhit_edep_sig"].Fill(hit.getEdep());
                energy_sum_sim += hit.getEdep();
                r = math.sqrt(x*x+y*y+z*z);
                thetasim = math.acos(z/r);
                self.histograms["simhit_theta_sig"].Fill(thetasim);
            for ih,hit in enumerate(self.hcalHits2):		
                if (hit.isNoise()==0 and hit.getX()!=0 and hit.getY()!=0 and hit.getZ()!=0):
                    recN +=1;
                    self.histograms["rechit_x_sig"].Fill(hit.getX());
                    self.histograms["rechit_y_sig"].Fill(hit.getY());
                    self.histograms["rechit_z_sig"].Fill(hit.getZ());
                    d = math.sqrt(hit.getX()*hit.getX()+hit.getY()*hit.getY()+hit.getZ()*hit.getZ());
                    thetarec = math.acos(hit.getZ()/d);
                    self.histograms["rechit_theta_sig"].Fill(thetarec);
                    self.histograms["rechit_edep_sig"].Fill(hit.getEnergy());
                    energy_sum_rec += hit.getEnergy();
            self.histograms["simhit_totE_sig"].Fill(energy_sum_sim);
            self.histograms["rechit_totE_sig"].Fill(energy_sum_rec);
            self.histograms["rechit_N_sig"].Fill(recN);
            self.histograms["simhit_N_sig"].Fill(simN);		
    
            
def makeCanvas(hists, tags, title, weights, logy=False, norm=False):

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
        hists[i].Scale(1/weights[i])
        if norm: hists[i].Scale(1./hists[i].Integral());
        if hists[i].GetMaximum() > tmax: tmax = hists[i].GetMaximum();
        if hists[i].GetMinimum() < tmin: tmin = hists[i].GetMinimum();

    c = ROOT.TCanvas("c","c",1000,800);
    hists[0].SetMinimum(1e-7);
    hists[0].SetMaximum(tmax*1.2);
    for i in range(len(hists)):
        if i == 0: hists[i].Draw("HIST");
        else: 
            hists[i].SetLineColor( colors[i] );
            hists[i].Draw('HISTsames');

    leg.Draw();
    if logy:
	    ROOT.gPad.SetLogy();
    c.SaveAs("compare_plot_HCAL"+str(title)+".root")

			
def main(options,args) : 
    sc = WabSigCompare(options.ifile1, options.ifile2, options.ofile,options.tag) ;

    histoE1 = [];
    tags = ["WAB, FF3", "Dark Photon Sig."];
    weights = [factor_wab, factor_signal];
    histoE1.append(sc.histograms["simhit_totE_wab"]);
    histoE1.append(sc.histograms["simhit_totE_sig"]);
    makeCanvas(histoE1, tags, "SimHitTotE",weights, True);
    histoX1 = [];
    histoX1.append(sc.histograms["simhit_x_wab"]);
    histoX1.append(sc.histograms["simhit_x_sig"]);
    makeCanvas(histoX1, tags, "SimHitposX", weights);
    histoY1 = [];
    histoY1.append(sc.histograms["simhit_y_wab"]);
    histoY1.append(sc.histograms["simhit_y_sig"]);
    makeCanvas(histoY1, tags, "SimHitposY", weights);
    histoZ1 = [];
    histoZ1.append(sc.histograms["simhit_z_wab"]);
    histoZ1.append(sc.histograms["simhit_z_sig"]);
    makeCanvas(histoZ1, tags, "SimHitposZ", weights);
    histoE2 = [];
    tags = ["WAB, FF3", "Dark Photon Sig."];
    histoE2.append(sc.histograms["rechit_totE_wab"]);
    histoE2.append(sc.histograms["rechit_totE_sig"]);
    makeCanvas(histoE2, tags, "RecHitTotE", weights, True);
    histoX2 = [];
    histoX2.append(sc.histograms["rechit_x_wab"]);
    histoX2.append(sc.histograms["rechit_x_sig"]);
    makeCanvas(histoX2, tags, "RecHitposX", weights);
    histoY2 = [];
    histoY2.append(sc.histograms["rechit_y_wab"]);
    histoY2.append(sc.histograms["rechit_y_sig"]);
    makeCanvas(histoY2, tags, "RecHitposY", weights);
    histoZ2 = [];
    histoZ2.append(sc.histograms["rechit_z_wab"]);
    histoZ2.append(sc.histograms["rechit_z_sig"]);
    makeCanvas(histoZ2, tags, "RecHitposZ", weights);
    sc.fout.Close();
if __name__ == "__main__":
	
    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    # input data files (4)
    parser.add_option('-a','--ifile1', dest='ifile1', default = 'file1.root',help='directory with data1', metavar='idir1')
    parser.add_option('-e','--ifile2', dest='ifile2', default = 'file2.root',help='directory with data2', metavar='idir2')
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
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    # Get the Event library 
    ROOT.gSystem.Load("/nfs/slac/g/ldmx/users/smidd/ldmx-sw/install/lib/libEvent.so");	

    main(options,args);			
	    
			
			
			
			
			
			
