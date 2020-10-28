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
        self.simParticles1 = ROOT.std.map(1,'ldmx::SimParticle');
        self.ecalHits1 = ROOT.std.vector('ldmx::EcalHit')();
        self.ecalSimHits1 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.tin1.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader1 ));
        self.tin1.SetBranchAddress("EcalSimHits_wab",  ROOT.AddressOf( self.ecalSimHits1 ));
        self.tin1.SetBranchAddress("EcalRecHits_wab",  ROOT.AddressOf( self.ecalHits1 ));
        #sig:
        self.evHeader2 = ROOT.ldmx.EventHeader()
        self.simParticles2 = ROOT.std.map(1,'ldmx::SimParticle');
        self.ecalHits2 = ROOT.std.vector('ldmx::EcalHit')();
        self.ecalSimHits2 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.tin2.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader2 ));
        self.tin2.SetBranchAddress("EcalSimHits_wab",  ROOT.AddressOf( self.ecalSimHits2 ));
        self.tin2.SetBranchAddress("EcalRecHits_wab",  ROOT.AddressOf( self.ecalHits2 ));

        self.histograms = {};

        #SimHits WAB:
        self.histograms["simhit_N_wab"]= ROOT.TH1F("ECAL simhit_HitsPerEvent_wab","; ECAL SimHits per Event; Weighted N", 400, 0,1000);
        self.histograms["simhit_x_wab"]= ROOT.TH1F("ECAL simhit_x_wab","; ECAL simhit x [mm]; Weighted N", 400, -500, 500);
        self.histograms["simhit_y_wab"] = ROOT.TH1F("ECAL simhit_y_wab","; ECAL simhit y [mm]; Weighted N", 400, -500,500);
        self.histograms["simhit_z_wab"] = ROOT.TH1F("ECAL simhit_z_wab","; ECAL simhit z [mm]; Weighted N", 400, 0, 800);
        self.histograms["simhit_edep_wab"] = ROOT.TH1F("ECAL Energy Dep per SimHit_wab","; ECAL Sim Hit Energy Dep [MeV]; Weighted N", 400, 0,2);
        self.histograms["simhit_totE_wab"]= ROOT.TH1F("ECAL Total Event E SimHits_wab","; ECAL Total Event Energy [MeV]; Weighted N", 400, 0,4000);
        self.histograms["simhit_theta_wab"] = ROOT.TH1F("ECAL simhit_theta_wab","; ECAL simhit Polar Angle [rad]; Weighted N", 400, 0, math.pi);
        
        # Sim Hits Signal:
        
        self.histograms["simhit_N_sig"]= ROOT.TH1F("ECAL simhit_HitsPerEvent_sig","; ECAL SimHits per Event; Weighted N", 400, 0,1000);
        self.histograms["simhit_x_sig"]= ROOT.TH1F("ECAL simhit_x_sig","; ECAL simhit x [mm]; Weighted N", 400, -500, 500);
        self.histograms["simhit_y_sig"] = ROOT.TH1F("ECAL simhit_y_sig","; ECAL simhit y [mm]; Weighted N", 400, -500,500);
        self.histograms["simhit_z_sig"] = ROOT.TH1F("ECAL simhit_z_sig","; ECAL simhit z [mm]; Weighted N", 400, 0, 500);
        self.histograms["simhit_edep_sig"] = ROOT.TH1F("ECAL Energy Dep per SimHit_sig","; ECAL Sim Hit Energy Dep [MeV]; Weighted N", 400, 0,2);
        self.histograms["simhit_totE_sig"]= ROOT.TH1F("ECAL Total Event E SimHits_sig","; ECAL Total Event Energy [MeV]; Weighted N", 400, 0,4000);
        self.histograms["simhit_theta_sig"] = ROOT.TH1F("ECAL simhit_theta_sig","; ECAL simhit Polar Angle [rad]; Weighted N", 400, 0, math.pi);
        
        #RecoHits WABS:
        self.histograms["rechit_N_wab"]= ROOT.TH1F("ECAL rechit_HitsPerEvent_wab","; ECAL RecHits per Event; Weighted N", 400, 0,100);
        self.histograms["rechit_x_wab"]= ROOT.TH1F("ECAL rechit_x_wab","; ECAL rechit x [mm]; Weighted N", 400, -500, 500);
        self.histograms["rechit_y_wab"] = ROOT.TH1F("ECAL rechit_y_wab","; ECAL rechit y [mm]; Weighted N", 400, -500, 500);
        self.histograms["rechit_z_wab"] = ROOT.TH1F("ECAL rechit_z_wab","; ECAL rechit z [mm]; Weighted N", 400, 0, 2000);
        self.histograms["rechit_totE_wab"]= ROOT.TH1F("ECAL Total Event E RecHits_wab","; ECAL Total Event Energy [MeV]; Weighted N", 400, 0,4000);
        self.histograms["rechit_edep_wab"] = ROOT.TH1F("ECAL Energy Dep per Rec Hits_wab","; ECAL RecHit Energy Dep [MeV]; Weighted N", 400, 0,10);
        self.histograms["rechit_theta_wab"] = ROOT.TH1F("ECAL rechit_theta_wab","; ECAL Rec Hit Polar Angle [rad]; Weighted N", 400, 0,math.pi);
        
        #RecoHits Signal:
        self.histograms["rechit_N_sig"]= ROOT.TH1F("ECAL rechit_HitsPerEvent_sig","; ECAL RecHits per Event; Weighted N", 400, 0,100);
        self.histograms["rechit_x_sig"]= ROOT.TH1F("ECAL rechit_x_sig","; ECAL rechit x [mm]; Weighted N", 400, -500, 500);
        self.histograms["rechit_y_sig"] = ROOT.TH1F("ECAL rechit_y_sig","; ECAL rechit y [mm]; Weighted N", 400, -500, 500);
        self.histograms["rechit_z_sig"] = ROOT.TH1F("ECAL rechit_z_sig","; ECAL rechit z [mm]; Weighted N", 400, 0, 2000);
        self.histograms["rechit_totE_sig"]= ROOT.TH1F("ECAL Total Event E RecHits_sig","; ECAL Total Event Energy [MeV]; Weighted N", 400, 0,4000);
        self.histograms["rechit_edep_sig"] = ROOT.TH1F("ECAL Energy Dep per Rec Hits_sig","; ECAL Rec Hit Energy Dep [MeV]; Weighted N", 400, 0,10);
        self.histograms["rechit_theta_sig"] = ROOT.TH1F("ECAL rechit_theta_sig","; ECAL Rec Hit Polar Angle [rad]; Weighted N", 400, 0,math.pi);

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
            for ih,hit in enumerate(self.ecalSimHits1):
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
            for ih,hit in enumerate(self.ecalHits1):		
                if (hit.isNoise()==0 and hit.getXPos()!=0 and hit.getYPos()!=0 and hit.getZPos()!=0):
                    recN +=1;
                    self.histograms["rechit_x_wab"].Fill(hit.getXPos());
                    self.histograms["rechit_y_wab"].Fill(hit.getYPos());
                    self.histograms["rechit_z_wab"].Fill(hit.getZPos());
                    d = math.sqrt(hit.getXPos()*hit.getXPos()+hit.getYPos()*hit.getYPos()+hit.getZPos()*hit.getZPos());
                    thetarec = math.acos(hit.getZPos()/d);
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
            for ih,hit in enumerate(self.ecalSimHits2):
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
            for ih,hit in enumerate(self.ecalHits2):		
                if (hit.isNoise()==0 and hit.getXPos()!=0 and hit.getYPos()!=0 and hit.getZPos()!=0):
                    recN +=1;
                    self.histograms["rechit_x_sig"].Fill(hit.getXPos());
                    self.histograms["rechit_y_sig"].Fill(hit.getYPos());
                    self.histograms["rechit_z_sig"].Fill(hit.getZPos());
                    d = math.sqrt(hit.getXPos()*hit.getXPos()+hit.getYPos()*hit.getYPos()+hit.getZPos()*hit.getZPos());
                    thetarec = math.acos(hit.getZPos()/d);
                    self.histograms["rechit_theta_sig"].Fill(thetarec);
                    self.histograms["rechit_edep_sig"].Fill(hit.getEnergy());
                    energy_sum_rec += hit.getEnergy();
            self.histograms["simhit_totE_sig"].Fill(energy_sum_sim);
            self.histograms["rechit_totE_sig"].Fill(energy_sum_rec);
            self.histograms["rechit_N_sig"].Fill(recN);
            self.histograms["simhit_N_sig"].Fill(simN);		
    
            
def makeCanvas(hists, tags, title, weights, logy=False, norm=False):

    colors = [2,4]

    leg = ROOT.TLegend(0.6,0.7,0.85,0.9);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.SetTextFont(42);
    tmax = -999;
    tmin = 9999;

    for i in range(len(hists)):
        leg.AddEntry(hists[i],tags[i],'l');
        hists[i].Scale(1/weights[i]);
        if norm: hists[i].Scale(1./hists[i].Integral());
        if hists[i].GetMaximum() > tmax: tmax = hists[i].GetMaximum();
        if hists[i].GetMinimum() < tmin: tmin = hists[i].GetMinimum();

    c = ROOT.TCanvas("c","c",1000,800);
    hists[0].SetMinimum(1e-7);
    hists[0].SetMaximum(tmax*1.2);
    for i in range(len(hists)):
        if i == 0: 
            hists[i].SetLineColor( colors[i] );
            hists[i].Draw('HIST');
        else: 
            hists[i].SetLineColor( colors[i] );
            hists[i].Draw('HISTsames');
    
    leg.Draw();
    if logy:
	    ROOT.gPad.SetLogy();
    c.SaveAs("compare_plot_ECAL_"+str(title)+".root")

			
def main(options,args) : 
    sc = WabSigCompare(options.ifile1, options.ifile2, options.ofile,options.tag) ;

    histoE = [];
    tags = [ "Dark Photon Signal","WAB, FF3"];
    weights = [factor_signal, factor_wab ];
    histoE.append(sc.histograms["simhit_totE_sig"]);
    histoE.append(sc.histograms["simhit_totE_wab"]);
    makeCanvas(histoE, tags, "SimHitTotE", weights, True);
    histoX = [];
    histoX.append(sc.histograms["simhit_x_sig"]);
    histoX.append(sc.histograms["simhit_x_wab"]);
    makeCanvas(histoX, tags, "SimHitposX", weights);
    histoY = [];
    histoY.append(sc.histograms["simhit_y_sig"]);
    histoY.append(sc.histograms["simhit_y_wab"]);
    makeCanvas(histoY, tags, "SimHitposY", weights);
    histoZ = [];
    histoZ.append(sc.histograms["simhit_z_sig"]);
    histoZ.append(sc.histograms["simhit_z_wab"]);
    makeCanvas(histoZ, tags, "SimHitposZ", weights);
    histoE = [];
    tags = ["Dark Photon Sig.", "WAB, FF3"];
    histoE.append(sc.histograms["rechit_totE_sig"]);
    histoE.append(sc.histograms["rechit_totE_wab"]);
    makeCanvas(histoE, tags, "RecHitTotE", weights, True);
    histoX = [];
    histoX.append(sc.histograms["rechit_x_sig"]);
    histoX.append(sc.histograms["rechit_x_wab"]);
    makeCanvas(histoX, tags, "RecHitposX", weights);
    histoY = [];
    histoY.append(sc.histograms["rechit_y_sig"]);
    histoY.append(sc.histograms["rechit_y_wab"]);
    makeCanvas(histoY, tags, "RecHitposY", weights);
    histoZ = [];
    histoZ.append(sc.histograms["rechit_z_sig"]);
    histoZ.append(sc.histograms["rechit_z_wab"]);
    makeCanvas(histoZ, tags, "RecHitposZ", weights);
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



			
			
			
			
