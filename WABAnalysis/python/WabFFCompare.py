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
proton_wabs = 74.
n_wab = 10000.

#FF1
xsec_FF1 = 2811000
event_weight_FF1 = xsec_FF1/n_wab
factor_wab_FF1 = event_weight_FF1 * proton_wabs * lumin

#FF2
xsec_FF2 = 2.394
event_weight_FF2 = xsec_FF2/n_wab
factor_wab_FF2 = event_weight_FF2 * proton_wabs * lumin

#FF3
xsec_FF3 = 1795.
event_weight_FF3 = xsec_FF3/n_wab
factor_wab_FF3 = event_weight_FF3 * proton_wabs * lumin

#FF4
xsec_FF4 = 15.83
event_weight_FF4 = xsec_FF4/n_wab
factor_wab_FF4 = event_weight_FF4 * proton_wabs * lumin

class WabFFCompare:
   
    def __init__(self, fn1, fn2, fn3, fn4, ofn, tag):

        #FF1:
        self.fin1 = ROOT.TFile(fn1);
        self.tin1 = self.fin1.Get("LDMX_Events")
        #FF2:
        self.fin2 = ROOT.TFile(fn2);
        self.tin2 = self.fin2.Get("LDMX_Events")
        #FF3:
        self.fin3 = ROOT.TFile(fn3);
        self.tin3 = self.fin3.Get("LDMX_Events")        
        #FF4:
        self.fin4 = ROOT.TFile(fn4);
        self.tin4 = self.fin4.Get("LDMX_Events")   
        
        self.tag = int(tag);

        self.fn_out = ofn;
        self.fout = ROOT.TFile("hist_"+self.fn_out,"RECREATE");
 
        #FF1:
        self.evHeader1 = ROOT.ldmx.EventHeader()
        self.simParticles1 = ROOT.std.map(1,'ldmx::SimParticle');
        self.hcalHits1 = ROOT.std.vector('ldmx::HcalHit')();
        self.hcalSimHits1 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.tin1.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader1 ));
        self.tin1.SetBranchAddress("HcalSimHits_wab",  ROOT.AddressOf( self.hcalSimHits1 ));
        self.tin1.SetBranchAddress("HcalRecHits_wab",  ROOT.AddressOf( self.hcalHits1 ));
        #FF2:
        self.evHeader2 = ROOT.ldmx.EventHeader()
        self.simParticles2 = ROOT.std.map(1,'ldmx::SimParticle');
        self.hcalHits2 = ROOT.std.vector('ldmx::HcalHit')();
        self.hcalSimHits2 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.tin2.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader2 ));
        self.tin2.SetBranchAddress("HcalSimHits_wab",  ROOT.AddressOf( self.hcalSimHits2 ));
        self.tin2.SetBranchAddress("HcalRecHits_wab",  ROOT.AddressOf( self.hcalHits2 ));
        #FF3:
        self.evHeader3 = ROOT.ldmx.EventHeader()
        self.simParticles3 = ROOT.std.map(1,'ldmx::SimParticle');
        self.hcalHits3 = ROOT.std.vector('ldmx::HcalHit')();
        self.hcalSimHits3 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.tin3.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader3 ));
        self.tin3.SetBranchAddress("HcalSimHits_wab",  ROOT.AddressOf( self.hcalSimHits3 ));
        self.tin3.SetBranchAddress("HcalRecHits_wab",  ROOT.AddressOf( self.hcalHits3 ));
        #FF4:
        self.evHeader4 = ROOT.ldmx.EventHeader()
        self.simParticles4 = ROOT.std.map(1,'ldmx::SimParticle');
        self.hcalHits4 = ROOT.std.vector('ldmx::HcalHit')();
        self.hcalSimHits4 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.tin4.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader4 ));
        self.tin4.SetBranchAddress("HcalSimHits_wab",  ROOT.AddressOf( self.hcalSimHits4 ));
        self.tin4.SetBranchAddress("HcalRecHits_wab",  ROOT.AddressOf( self.hcalHits4 ));

        self.histogramsFF1 = {};
        self.histogramsFF2 = {};
        self.histogramsFF3 = {};
        self.histogramsFF4 = {};
        #SimHits:
        self.histogramsFF1["simhit_N"]= ROOT.TH1F("simhit_HitsPerEvent","; SimHits per Event; N", 400, 0,1000);
        self.histogramsFF1["simhit_x"]= ROOT.TH1F("simhit_x","; simhit x [mm]; N", 400, -800, 800);
        self.histogramsFF1["simhit_y"] = ROOT.TH1F("simhit_y","; simhit y [mm]; N", 400, -800,800);
        self.histogramsFF1["simhit_z"] = ROOT.TH1F("simhit_z","; simhit z [mm]; N", 400, 0, 800);
        self.histogramsFF1["simhit_edep"] = ROOT.TH1F("Energy Dep per SimHit","; Sim Hit Energy Dep [MeV]; N", 400, 0,2);
        self.histogramsFF1["simhit_totE"]= ROOT.TH1F("Total Event E SimHits","; Total Event Energy [MeV]; N", 400, 0,400);
        self.histogramsFF1["simhit_theta"] = ROOT.TH1F("simhit_theta","; simhit Polar Angle [rad]; N", 400, 0, math.pi);
        
        self.histogramsFF2["simhit_N"]= ROOT.TH1F("simhit_HitsPerEvent","; SimHits per Event; N", 400, 0,1000);
        self.histogramsFF2["simhit_x"]= ROOT.TH1F("simhit_x","; simhit x [mm]; N", 400, -800, 800);
        self.histogramsFF2["simhit_y"] = ROOT.TH1F("simhit_y","; simhit y [mm]; N", 400, -800,800);
        self.histogramsFF2["simhit_z"] = ROOT.TH1F("simhit_z","; simhit z [mm]; N", 400, 0, 800);
        self.histogramsFF2["simhit_edep"] = ROOT.TH1F("Energy Dep per SimHit","; Sim Hit Energy Dep [MeV]; N", 400, 0,2);
        self.histogramsFF2["simhit_totE"]= ROOT.TH1F("Total Event E SimHits","; Total Event Energy [MeV]; N", 400, 0,400);
        self.histogramsFF2["simhit_theta"] = ROOT.TH1F("simhit_theta","; simhit Polar Angle [rad]; N", 400, 0, math.pi);
        
        self.histogramsFF3["simhit_N"]= ROOT.TH1F("simhit_HitsPerEvent","; SimHits per Event; N", 400, 0,1000);
        self.histogramsFF3["simhit_x"]= ROOT.TH1F("simhit_x","; simhit x [mm]; N", 400, -800, 800);
        self.histogramsFF3["simhit_y"] = ROOT.TH1F("simhit_y","; simhit y [mm]; N", 400, -800,800);
        self.histogramsFF3["simhit_z"] = ROOT.TH1F("simhit_z","; simhit z [mm]; N", 400, 0, 800);
        self.histogramsFF3["simhit_edep"] = ROOT.TH1F("Energy Dep per SimHit","; Sim Hit Energy Dep [MeV]; N", 400, 0,2);
        self.histogramsFF3["simhit_totE"]= ROOT.TH1F("Total Event E SimHits","; Total Event Energy [MeV]; N", 400, 0,400);
        self.histogramsFF3["simhit_theta"] = ROOT.TH1F("simhit_theta","; simhit Polar Angle [rad]; N", 400, 0, math.pi);

        self.histogramsFF4["simhit_N"]= ROOT.TH1F("simhit_HitsPerEvent","; SimHits per Event; N", 400, 0,1000);
        self.histogramsFF4["simhit_x"]= ROOT.TH1F("simhit_x","; simhit x [mm]; N", 400, -800, 800);
        self.histogramsFF4["simhit_y"] = ROOT.TH1F("simhit_y","; simhit y [mm]; N", 400, -800,800);
        self.histogramsFF4["simhit_z"] = ROOT.TH1F("simhit_z","; simhit z [mm]; N", 400, 0, 800);
        self.histogramsFF4["simhit_edep"] = ROOT.TH1F("Energy Dep per SimHit","; Sim Hit Energy Dep [MeV]; N", 400, 0,2);
        self.histogramsFF4["simhit_totE"]= ROOT.TH1F("Total Event E SimHits","; Total Event Energy [MeV]; N", 400, 0,400);
        self.histogramsFF4["simhit_theta"] = ROOT.TH1F("simhit_theta","; simhit Polar Angle [rad]; N", 400, 0, math.pi);
        
        
        #RecoHits
        self.histogramsFF1["rechit_N"]= ROOT.TH1F("rechit_HitsPerEvent","; RecHits per Event; N", 400, 0,100);
        self.histogramsFF1["rechit_x"]= ROOT.TH1F("rechit_x","; rechit x [mm]; N", 400, -2000, 2000);
        self.histogramsFF1["rechit_y"] = ROOT.TH1F("rechit_y","; rechit y [mm]; N", 400, -2000, 2000);
        self.histogramsFF1["rechit_z"] = ROOT.TH1F("rechit_z","; rechit z [mm]; N", 400, 0, 4000);
        self.histogramsFF1["rechit_totE"]= ROOT.TH1F("Total Event E RecHits","; Total Event Energy [MeV]; N", 400, 0,600);
        self.histogramsFF1["rechit_edep"] = ROOT.TH1F(" Energy Dep per Rec Hits","; Rec Hit Energy Dep [MeV]; N", 400, 0,10);
        self.histogramsFF1["rechit_theta"] = ROOT.TH1F(" rechit_theta","; Rec Hit Polar Angle [rad]; N", 400, 0,math.pi);
        
        self.histogramsFF2["rechit_N"]= ROOT.TH1F("rechit_HitsPerEvent","; RecHits per Event; N", 400, 0,100);
        self.histogramsFF2["rechit_x"]= ROOT.TH1F("rechit_x","; rechit x [mm]; N", 400, -2000, 2000);
        self.histogramsFF2["rechit_y"] = ROOT.TH1F("rechit_y","; rechit y [mm]; N", 400, -2000, 2000);
        self.histogramsFF2["rechit_z"] = ROOT.TH1F("rechit_z","; rechit z [mm]; N", 400, 0, 4000);
        self.histogramsFF2["rechit_totE"]= ROOT.TH1F("Total Event E RecHits","; Total Event Energy [MeV]; N", 400, 0,600);
        self.histogramsFF2["rechit_edep"] = ROOT.TH1F(" Energy Dep per Rec Hits","; Rec Hit Energy Dep [MeV]; N", 400, 0,10);
        self.histogramsFF2["rechit_theta"] = ROOT.TH1F(" rechit_theta","; Rec Hit Polar Angle [rad]; N", 400, 0,math.pi);
        
        self.histogramsFF3["rechit_N"]= ROOT.TH1F("rechit_HitsPerEvent","; RecHits per Event; N", 400, 0,100);
        self.histogramsFF3["rechit_x"]= ROOT.TH1F("rechit_x","; rechit x [mm]; N", 400, -2000, 2000);
        self.histogramsFF3["rechit_y"] = ROOT.TH1F("rechit_y","; rechit y [mm]; N", 400, -2000, 2000);
        self.histogramsFF3["rechit_z"] = ROOT.TH1F("rechit_z","; rechit z [mm]; N", 400, 0, 4000);
        self.histogramsFF3["rechit_totE"]= ROOT.TH1F("Total Event E RecHits","; Total Event Energy [MeV]; N", 400, 0,600);
        self.histogramsFF3["rechit_edep"] = ROOT.TH1F(" Energy Dep per Rec Hits","; Rec Hit Energy Dep [MeV]; N", 400, 0,10);
        self.histogramsFF3["rechit_theta"] = ROOT.TH1F(" rechit_theta","; Rec Hit Polar Angle [rad]; N", 400, 0,math.pi);

        self.histogramsFF4["rechit_N"]= ROOT.TH1F("rechit_HitsPerEvent","; RecHits per Event; N", 400, 0,100);
        self.histogramsFF4["rechit_x"]= ROOT.TH1F("rechit_x","; rechit x [mm]; N", 400, -2000, 2000);
        self.histogramsFF4["rechit_y"] = ROOT.TH1F("rechit_y","; rechit y [mm]; N", 400, -2000, 2000);
        self.histogramsFF4["rechit_z"] = ROOT.TH1F("rechit_z","; rechit z [mm]; N", 400, 0, 4000);
        self.histogramsFF4["rechit_totE"]= ROOT.TH1F("Total Event E RecHits","; Total Event Energy [MeV]; N", 400, 0,600);
        self.histogramsFF4["rechit_edep"] = ROOT.TH1F(" Energy Dep per Rec Hits","; Rec Hit Energy Dep [MeV]; N", 400, 0,10);
        self.histogramsFF4["rechit_theta"] = ROOT.TH1F(" rechit_theta","; Rec Hit Polar Angle [rad]; N", 400, 0,math.pi);

        self.loop();
        self.writeOutHistos();

    def writeOutHistos(self):

        self.fout.cd();
        for key in self.histogramsFF1: 
            print(key);
            self.histogramsFF1[key].Write();
        for key in self.histogramsFF2: 
            print(key);
            self.histogramsFF2[key].Write();
        for key in self.histogramsFF3: 
            print(key);
            self.histogramsFF3[key].Write();
        for key in self.histogramsFF4: 
            print(key);
            self.histogramsFF4[key].Write();
            
    def loop(self):
        nentFF1 = self.tin1.GetEntriesFast();
        print(nentFF1)
        for i in range(nentFF1):
            self.tin1.GetEntry(i);
            energy_sum_sim = 0;
            energy_sum_rec = 0;
            simN = 0;
            recN = 0;
            for ih,hit in enumerate(self.hcalSimHits1):
                [x, y, z] = hit.getPosition();	
                simN +=1;
                self.histogramsFF1["simhit_x"].Fill(x);
                self.histogramsFF1["simhit_y"].Fill(y);
                self.histogramsFF1["simhit_z"].Fill(z);
                self.histogramsFF1["simhit_edep"].Fill(hit.getEdep());
                energy_sum_sim += hit.getEdep();
                r = math.sqrt(x*x+y*y+z*z);
                thetasim = math.acos(z/r);
                self.histogramsFF1["simhit_theta"].Fill(thetasim);
            for ih,hit in enumerate(self.hcalHits1):		
                if (hit.isNoise()==0 and hit.getX()!=0 and hit.getY()!=0 and hit.getZ()!=0):
                    recN +=1;
                    self.histogramsFF1["rechit_x"].Fill(hit.getX());
                    self.histogramsFF1["rechit_y"].Fill(hit.getY());
                    self.histogramsFF1["rechit_z"].Fill(hit.getZ());
                    d = math.sqrt(hit.getX()*hit.getX()+hit.getY()*hit.getY()+hit.getZ()*hit.getZ());
                    thetarec = math.acos(hit.getZ()/d);
                    self.histogramsFF1["rechit_theta"].Fill(thetarec);
                    self.histogramsFF1["rechit_edep"].Fill(hit.getEnergy());
                    energy_sum_rec += hit.getEnergy();
            self.histogramsFF1["simhit_totE"].Fill(energy_sum_sim);
            self.histogramsFF1["rechit_totE"].Fill(energy_sum_rec);
            self.histogramsFF1["rechit_N"].Fill(recN);
            self.histogramsFF1["simhit_N"].Fill(simN);
            	
        nentFF2 = self.tin2.GetEntriesFast();
        print(nentFF2)
        for i in range(nentFF2):
            self.tin2.GetEntry(i);
            energy_sum_sim = 0;
            energy_sum_rec = 0;
            simN = 0;
            recN = 0;
            for ih,hit in enumerate(self.hcalSimHits2):
                [x, y, z] = hit.getPosition();	
                simN +=1;
                self.histogramsFF2["simhit_x"].Fill(x);
                self.histogramsFF2["simhit_y"].Fill(y);
                self.histogramsFF2["simhit_z"].Fill(z);
                self.histogramsFF2["simhit_edep"].Fill(hit.getEdep());
                energy_sum_sim += hit.getEdep();
                r = math.sqrt(x*x+y*y+z*z);
                thetasim = math.acos(z/r);
                self.histogramsFF2["simhit_theta"].Fill(thetasim);
            for ih,hit in enumerate(self.hcalHits2):		
                if (hit.isNoise()==0 and hit.getX()!=0 and hit.getY()!=0 and hit.getZ()!=0):
                    recN +=1;
                    self.histogramsFF2["rechit_x"].Fill(hit.getX());
                    self.histogramsFF2["rechit_y"].Fill(hit.getY());
                    self.histogramsFF2["rechit_z"].Fill(hit.getZ());
                    d = math.sqrt(hit.getX()*hit.getX()+hit.getY()*hit.getY()+hit.getZ()*hit.getZ());
                    thetarec = math.acos(hit.getZ()/d);
                    self.histogramsFF2["rechit_theta"].Fill(thetarec);
                    self.histogramsFF2["rechit_edep"].Fill(hit.getEnergy());
                    energy_sum_rec += hit.getEnergy();
            self.histogramsFF2["simhit_totE"].Fill(energy_sum_sim);
            self.histogramsFF2["rechit_totE"].Fill(energy_sum_rec);
            self.histogramsFF2["rechit_N"].Fill(recN);
            self.histogramsFF2["simhit_N"].Fill(simN);	
            
        nentFF3 = self.tin3.GetEntriesFast();
        print(nentFF3)
        for i in range(nentFF3):
            self.tin3.GetEntry(i);
            energy_sum_sim = 0;
            energy_sum_rec = 0;
            simN = 0;
            recN = 0;
            for ih,hit in enumerate(self.hcalSimHits3):
                [x, y, z] = hit.getPosition();	
                simN +=1;
                self.histogramsFF3["simhit_x"].Fill(x);
                self.histogramsFF3["simhit_y"].Fill(y);
                self.histogramsFF3["simhit_z"].Fill(z);
                self.histogramsFF3["simhit_edep"].Fill(hit.getEdep());
                energy_sum_sim += hit.getEdep();
                r = math.sqrt(x*x+y*y+z*z);
                thetasim = math.acos(z/r);
                self.histogramsFF3["simhit_theta"].Fill(thetasim);
            for ih,hit in enumerate(self.hcalHits3):		
                if (hit.isNoise()==0 and hit.getX()!=0 and hit.getY()!=0 and hit.getZ()!=0):
                    recN +=1;
                    self.histogramsFF3["rechit_x"].Fill(hit.getX());
                    self.histogramsFF3["rechit_y"].Fill(hit.getY());
                    self.histogramsFF3["rechit_z"].Fill(hit.getZ());
                    d = math.sqrt(hit.getX()*hit.getX()+hit.getY()*hit.getY()+hit.getZ()*hit.getZ());
                    thetarec = math.acos(hit.getZ()/d);
                    self.histogramsFF3["rechit_theta"].Fill(thetarec);
                    self.histogramsFF3["rechit_edep"].Fill(hit.getEnergy());
                    energy_sum_rec += hit.getEnergy();
            self.histogramsFF3["simhit_totE"].Fill(energy_sum_sim);
            self.histogramsFF3["rechit_totE"].Fill(energy_sum_rec);
            self.histogramsFF3["rechit_N"].Fill(recN);
        self.histogramsFF3["simhit_N"].Fill(simN);
        
        nentFF4 = self.tin4.GetEntriesFast();
        print(nentFF4)
        for i in range(nentFF4):
            self.tin4.GetEntry(i);
            energy_sum_sim = 0;
            energy_sum_rec = 0;
            simN = 0;
            recN = 0;
            for ih,hit in enumerate(self.hcalSimHits4):
                [x, y, z] = hit.getPosition();	
                simN +=1;
                self.histogramsFF4["simhit_x"].Fill(x);
                self.histogramsFF4["simhit_y"].Fill(y);
                self.histogramsFF4["simhit_z"].Fill(z);
                self.histogramsFF4["simhit_edep"].Fill(hit.getEdep());
                energy_sum_sim += hit.getEdep();
                r = math.sqrt(x*x+y*y+z*z);
                thetasim = math.acos(z/r);
                self.histogramsFF3["simhit_theta"].Fill(thetasim);
            for ih,hit in enumerate(self.hcalHits4):		
                if (hit.isNoise()==0 and hit.getX()!=0 and hit.getY()!=0 and hit.getZ()!=0):
                    recN +=1;
                    self.histogramsFF4["rechit_x"].Fill(hit.getX());
                    self.histogramsFF4["rechit_y"].Fill(hit.getY());
                    self.histogramsFF4["rechit_z"].Fill(hit.getZ());
                    d = math.sqrt(hit.getX()*hit.getX()+hit.getY()*hit.getY()+hit.getZ()*hit.getZ());
                    thetarec = math.acos(hit.getZ()/d);
                    self.histogramsFF4["rechit_theta"].Fill(thetarec);
                    self.histogramsFF4["rechit_edep"].Fill(hit.getEnergy());
                    energy_sum_rec += hit.getEnergy();
            self.histogramsFF4["simhit_totE"].Fill(energy_sum_sim);
            self.histogramsFF4["rechit_totE"].Fill(energy_sum_rec);
            self.histogramsFF4["rechit_N"].Fill(recN);
        self.histogramsFF4["simhit_N"].Fill(simN);	
        
def makeCanvas(hists, tags, title, weights, logy=False, norm=False):

    colors = [2,3,6,7,8]

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
    c.SaveAs("compare_plot_HCALFF"+str(title)+".root")
			
def main(options,args) : 
    sc = WabFFCompare(options.ifile1, options.ifile2, options.ifile3, options.ifile4, options.ofile,options.tag) ;
    histoE = [];
    tags = ["FF1", "FF3", "FF4", "FF2"];
    weights = [factor_wab_FF1, factor_wab_FF3, factor_wab_FF4, factor_wab_FF2];
    histoE.append(sc.histogramsFF1["rechit_totE"]);
    histoE.append(sc.histogramsFF3["rechit_totE"]);
    histoE.append(sc.histogramsFF4["rechit_totE"]);
    histoE.append(sc.histogramsFF2["rechit_totE"]);
    makeCanvas(histoE, tags, "RecHitsTotE", weights, True);
    
    histoX = [];
    histoX.append(sc.histogramsFF1["rechit_x"]);
    histoX.append(sc.histogramsFF3["rechit_x"]);
    histoX.append(sc.histogramsFF4["rechit_x"]);
    histoX.append(sc.histogramsFF2["rechit_x"]);
    makeCanvas(histoX, tags, "RecHitsTotx", weights, True);
    
    histoY = [];
    histoY.append(sc.histogramsFF1["rechit_y"]);
    histoY.append(sc.histogramsFF3["rechit_y"]);
    histoY.append(sc.histogramsFF4["rechit_y"]);
    histoY.append(sc.histogramsFF2["rechit_y"]);
    makeCanvas(histoY, tags, "RecHitsToty", weights, True);
    
    histoZ = [];
    histoZ.append(sc.histogramsFF1["rechit_z"]);
    histoZ.append(sc.histogramsFF3["rechit_z"]);
    histoZ.append(sc.histogramsFF4["rechit_z"]);
    histoZ.append(sc.histogramsFF2["rechit_z"]);
    makeCanvas(histoZ, tags, "RecHitsTotz", weights, True);
    
    sc.fout.Close();
if __name__ == "__main__":
	
    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    # input data files (4)
    parser.add_option('-a','--ifile1', dest='ifile1', default = 'file1.root',help='directory with data1', metavar='idir1')
    parser.add_option('-e','--ifile2', dest='ifile2', default = 'file2.root',help='directory with data2', metavar='idir2')
    parser.add_option('-c','--ifile3', dest='ifile3', default = 'file3.root',help='directory with data3', metavar='idir3')
    parser.add_option('-f','--ifile4', dest='ifile4', default = 'file4.root',help='directory with data4', metavar='idir4')
    # 1 output data file
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
		    
		    
		    
			
			
			
			
