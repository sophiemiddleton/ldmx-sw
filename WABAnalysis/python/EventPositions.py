#!/usr/bin/python
import argparse
import importlib
import ROOT
from ROOT import TTree, TBranch
ROOT.gSystem.Load("/nfs/slac/g/ldmx/users/smidd/NewLDMX/ldmx-sw/install/lib/libEvent.so")
import os
import math
import sys
#import numpy as np
import matplotlib.pyplot as plt
from array import array
from optparse import OptionParser
sys.path.insert(0, '../')

lumin = 0.22

#signal
proton_signal = 74.
xsec_signal = 3678.1344
n_signal = 10000
event_weight_signal = xsec_signal/n_signal
factor_signal = event_weight_signal * proton_signal * lumin

#wabs
proton_wabs = 74.
xsec_wab = 1795.
n_wab = 10000.
lumin_wab = 0.22
event_weight_wab = xsec_wab/n_wab
factor_wab = event_weight_wab * proton_wabs * lumin

class Event():
    def __init__(self):
        self.x = []
        self.y = []
        self.z = []
        self.e = []
        self.pid = []
        
class WabEvent:
   
    def __init__(self, fn1, ofn, tag):

        self.fin1 = ROOT.TFile(fn1);
        self.tin1 = self.fin1.Get("LDMX_Events")

        self.tag = int(tag);

        self.fn_out = ofn;
        self.fout = ROOT.TFile("hist_"+self.fn_out,"RECREATE");
        self.EventList = []
        self.RecHitList = []
        self.HCALPlaneHitListElectron = []
        self.HCALPlaneHitListPhoton = []
        self.ECALPlaneHitListElectron = []
        self.ECALPlaneHitListPhoton = []
        self.TargetPlaneHitListElectron = []
        self.TargetPlaneHitListPhoton = []
        #wab:
        self.evHeader1 = ROOT.ldmx.EventHeader()
        self.simParticles1 = ROOT.std.map(int,'ldmx::SimParticle')();
        self.hcalHits1 = ROOT.std.vector('ldmx::HcalHit')();
        self.hcalSimHits1 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.hcalScoringPlaneHits1 = ROOT.std.vector('ldmx::SimTrackerHit')();
        self.ecalScoringPlaneHits1 = ROOT.std.vector('ldmx::SimTrackerHit')();
        self.targetScoringPlaneHits1 = ROOT.std.vector('ldmx::SimTrackerHit')();
        self.tin1.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader1 ));
        self.tin1.SetBranchAddress("HcalSimHits_v12",  ROOT.AddressOf( self.hcalSimHits1 ));
        self.tin1.SetBranchAddress("HcalRecHits_v12",  ROOT.AddressOf( self.hcalHits1 ));
        self.tin1.SetBranchAddress("HcalScoringPlaneHits_v12",  ROOT.AddressOf(self.hcalScoringPlaneHits1));
        self.tin1.SetBranchAddress("EcalScoringPlaneHits_v12",  ROOT.AddressOf(self.ecalScoringPlaneHits1));
        self.tin1.SetBranchAddress("TargetScoringPlaneHits_v12",  ROOT.AddressOf(self.targetScoringPlaneHits1));
        self.tin1.SetBranchAddress("SimParticles_v12",  ROOT.AddressOf( self.simParticles1 ));
        self.RecHits_x = []
        self.RecHits_y = []
        self.RecHits_z = []
        
        self.SimHits_x = []
        self.SimHits_y = []
        self.SimHits_z = []
        
        self.Elec_HCALSPHits_x = []
        self.Elec_HCALSPHits_y = []
        self.Elec_HCALSPHits_z = []
        
        self.Elec_ECALSPHits_x = []
        self.Elec_ECALSPHits_y = []
        self.Elec_ECALSPHits_z = []
        
        self.Phot_HCALSPHits_x = []
        self.Phot_HCALSPHits_y = []
        self.Phot_HCALSPHits_z = []
        
        self.Phot_ECALSPHits_x = []
        self.Phot_ECALSPHits_y = []
        self.Phot_ECALSPHits_z = []
        
        
        self.loop();
        self.writeOutHistos();

    def writeOutHistos(self):

        self.fout.cd();
        
    def loop(self):
        partList = []
        nentWAB = self.tin1.GetEntriesFast();

      
        """
        for i in range(nentWAB):
            Electron = Event()
            Photon = Event()
            self.tin1.GetEntry(i);
            for ih,hit in enumerate(self.targetScoringPlaneHits1):
                j = hit.getID() & 0xFFF;
                if j ==1 and hit.getTrackID()==1:
                    #print("Track ID", i, hit.getTrackID(), hit.getPdgID());
                    
                    Electron.x.append(hit.getPosition()[0])
                    Electron.y.append(hit.getPosition()[1])
                    Electron.z.append(hit.getPosition()[2])
                    Electron.e.append(hit.getEnergy())
                    Electron.pid.append(hit.getPdgID())
            
                if j ==1 and hit.getTrackID()==2:
                    Photon.x.append(hit.getPosition()[0])
                    Photon.y.append(hit.getPosition()[1])
                    Photon.z.append(hit.getPosition()[2])
                    Photon.e.append(hit.getEnergy())
                    Photon.pid.append(hit.getPdgID())
            self.TargetPlaneHitListElectron.append(Electron) 
            self.TargetPlaneHitListPhoton.append(Photon)  
        """
        for i in range(10):
           
            self.tin1.GetEntry(i);
            for ih,hit in enumerate(self.hcalScoringPlaneHits1):
                j = hit.getID() & 0xFFF;
                if hit.getTrackID()==1:
                    self.Elec_HCALSPHits_x.append(hit.getPosition()[0])
                    self.Elec_HCALSPHits_y.append(hit.getPosition()[1])
                    self.Elec_HCALSPHits_z.append(hit.getPosition()[2])
                    
                if hit.getTrackID()==2:
                    self.Phot_HCALSPHits_x.append(hit.getPosition()[0])
                    self.Phot_HCALSPHits_y.append(hit.getPosition()[1])
                    self.Phot_HCALSPHits_z.append(hit.getPosition()[2])
                    
            
        for i in range(10):
           
            self.tin1.GetEntry(i);
            for ih,hit in enumerate(self.ecalScoringPlaneHits1):
                if hit.getTrackID()==1:
                    self.Elec_ECALSPHits_x.append(hit.getPosition()[0])
                    self.Elec_ECALSPHits_y.append(hit.getPosition()[1])
                    self.Elec_ECALSPHits_z.append(hit.getPosition()[2])
                    
                if hit.getTrackID()==2:
                    self.Phot_ECALSPHits_x.append(hit.getPosition()[0])
                    self.Phot_ECALSPHits_y.append(hit.getPosition()[1])
                    self.Phot_ECALSPHits_z.append(hit.getPosition()[2])
        
            
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
    c.SaveAs("compare_plot_HCALWas"+str(title)+".root")

def make_single(hist,  title, logy=False, norm=False):
    c = ROOT.TCanvas("c","c",1000,800);
    hist.GetXaxis().SetTitle("x [mm]")
    hist.GetYaxis().SetTitle("y [mm]")
    hist.GetZaxis().SetTitle("z [mm]")
    hist.Draw("COLZ")
    c.SaveAs("event"+str(title)+".root")
    
		
def main(options,args) : 
    sc = WabEvent(options.ifile1, options.ofile,options.tag) ;
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.RecHits_x,
            bins=50, 
            range=(-800,800), 
            label="HCAL Reco Hits (x)")
   
    ax.set_xlabel('HCAL Rec Hit x')  
    plt.savefig("HCALRecHitEvent_x.png")
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.RecHits_y,
            bins=50, 
            range=(-800,800), 
            label="HCAL Reco Hits (y)")
   
    ax.set_xlabel('HCAL Rec Hit y')  
    plt.savefig("HCALRecHitEvent_y.png")
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.RecHits_z,
            bins=50, 
            range=(0,2000), 
            label="HCAL Reco Hits (z)")
    
    ax.set_xlabel('HCAL Rec Hit z')  
    plt.savefig("HCALRecHitEvent_z.png")
    
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.SimHits_x,
            bins=50, 
            range=(-800,800), 
            label="HCAL Sim Hits (x)")
   
    ax.set_xlabel('HCAL Sim Hit x')  
    plt.savefig("HCALSimHitEvent_x.png")
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.SimHits_y,
            bins=50, 
            range=(-800,800), 
            label="HCAL Sim Hits (y)")
   
    ax.set_xlabel('HCAL Sim Hit y')  
    plt.savefig("HCALSimHitEvent_y.png")
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.SimHits_z,
            bins=50, 
            range=(0,2000), 
            label="HCAL Sim Hits (z)")
   
    ax.set_xlabel('HCAL Sim Hit z')  
    plt.savefig("HCALSimHitEvent_z.png")
    
    
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.Phot_HCALSPHits_x,
            bins=50, 
            range=(-800,800), 
            label="HCAL SP WAB Photons (x)")
   
    ax.set_xlabel('HCAL SP WAB Photons  x')  
    plt.savefig("HCALSPWABPhotons_x.png")
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.Phot_HCALSPHits_y,
            bins=50, 
            range=(-800,800), 
            label="HCAL SP WAB Photons  (y)")
   
    ax.set_xlabel('HCAL SP WAB Photons  y')  
    plt.savefig("HCALSPWABPhotons_y.png")
    
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.Phot_ECALSPHits_x,
            bins=50, 
            range=(-800,800), 
            label="ECAL SP WAB Photons (x)")
   
    ax.set_xlabel('ECAL SP WAB Photons  x')  
    plt.savefig("ECALSPWABPhotons_x.png")
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.Phot_ECALSPHits_y,
            bins=50, 
            range=(-800,800), 
            label="ECAL SP WAB Photons  (y)")
   
    ax.set_xlabel('ECAL SP WAB Photons  y')  
    plt.savefig("ECALSPWABPhotons_y.png")

    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.Elec_HCALSPHits_x,
            bins=50, 
            range=(-800,800), 
            label="HCAL SP WAB Photons (x)")
   
    ax.set_xlabel('HCAL SP WAB Electrons  x')  
    plt.savefig("HCALSPWABElectrons_x.png")
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.Elec_HCALSPHits_y,
            bins=50, 
            range=(-800,800), 
            label="HCAL SP WAB Electrons  (y)")
   
    ax.set_xlabel('HCAL SP WAB Electrons  y')  
    plt.savefig("HCALSPWABElectrons_y.png")
    
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.Elec_ECALSPHits_x,
            bins=50, 
            range=(-800,800), 
            label="ECAL SP WAB Electrons (x)")
   
    ax.set_xlabel('ECAL SP WAB Electrons  x')  
    plt.savefig("ECALSPWABElectrons_x.png")
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.Elec_ECALSPHits_y,
            bins=50, 
            range=(-800,800), 
            label="ECAL SP WAB Electrons  (y)")
   
    ax.set_xlabel('ECAL SP WAB Electrons  y')  
    plt.savefig("ECALSPWABElectrons_y.png")
    """
    for i, event in enumerate(sc.EventList):
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        p = ax.scatter3D(event.x, event.z, event.y, c=event.e, cmap='inferno')
        ax.set_xlabel('HCAL Sim Hit x')
        ax.set_ylabel('HCAL Sim Hit z')
        ax.set_zlabel('HCAL Sim Hit y')    
        fig.suptitle('HCAL Sim Hits Event '+str(i))  
        fig.colorbar(p)
        plt.savefig("SimHitEvent"+str(i)+".png")
        
    for i, event in enumerate(sc.HCALPlaneHitListElectron):
        fig = plt.figure()
        ax = plt.axes()
        p = ax.scatter(event.x, event.y, c=event.e, cmap='inferno')
        ax.set_xlabel('HCAL Plane Hit x')
        ax.set_ylabel('HCAL Plane Hit y')
        #ax.set_zlabel('HCAL Plane Hit y')  
        fig.suptitle('HCAL Primary Electron at Scoring Plane Event '+str(i))    
        fig.colorbar(p)
        plt.savefig("HCALPlaneHitEvent_Electron"+str(i)+".png")
        
    for i, event in enumerate(sc.HCALPlaneHitListPhoton):
        fig = plt.figure()
        ax = plt.axes()
        p = ax.scatter(event.x, event.y, c=event.e, cmap='inferno')
        ax.set_xlabel('HCAL Plane Hit x')
        ax.set_ylabel('HCAL Plane Hit y')
        fig.suptitle('HCAL Primary Photon at Scoring Plane Event '+str(i))
        fig.colorbar(p)
        plt.savefig("HCALPlaneHitEvent_Photon"+str(i)+".png")
        
    for i, event in enumerate(sc.ECALPlaneHitListElectron):
        fig = plt.figure()
        ax = plt.axes()
        p = ax.scatter(event.x, event.y, c=event.e, cmap='inferno')
        ax.set_xlabel('ECAL Plane Hit x')
        ax.set_ylabel('ECAL Plane Hit y')
        fig.suptitle('ECAL Primary Electron at Scoring Plane Event '+str(i))  
        fig.colorbar(p)
        plt.savefig("ECALPlaneHitEvent_Electron"+str(i)+".png")
        
    for i, event in enumerate(sc.ECALPlaneHitListPhoton):
        fig = plt.figure()
        ax = plt.axes()
        p = ax.scatter(event.x, event.y, c=event.e, cmap='inferno')
        ax.set_xlabel('ECAL Plane Hit x')
        ax.set_ylabel('ECAL Plane Hit y')
        fig.suptitle('ECAL Primary Photon at Scoring Plane Event '+str(i)) 
        fig.colorbar(p)
        plt.savefig("ECALPlaneHitEvent_Photon"+str(i)+".png")
        
        
    for i, event in enumerate(sc.TargetPlaneHitListElectron):
        fig = plt.figure()
        ax = plt.axes()
        p = ax.scatter(event.x, event.y)
        ax.set_xlabel('Target Plane 1 Hit x')
        ax.set_ylabel('Target Plane 1 Hit y')
        fig.suptitle('Target Primary Electron at Scoring Plane Event '+str(i))    
        fig.colorbar(p)
        plt.savefig("TargetPlaneHitEvent_Electrons"+str(i)+".png")
        
    for i, event in enumerate(sc.TargetPlaneHitListPhoton):
        fig = plt.figure()
        ax = plt.axes()
        p = ax.scatter(event.x, event.y)
        ax.set_xlabel('Target Plane 1 Hit x')
        ax.set_ylabel('Target Plane 1 Hit y')
        fig.suptitle('Target Primary Photon at Scoring Plane Event '+str(i))
        fig.colorbar(p)
        plt.savefig("TargetPlaneHitEvent_Photons"+str(i)+".png")
    """
    sc.fout.Close();
    
if __name__ == "__main__":
	
    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    # input data files (4)
    parser.add_option('-a','--ifile1', dest='ifile1', default = 'file1.root',help='directory with data1', metavar='idir1')
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
    ROOT.gSystem.Load("/nfs/slac/g/ldmx/users/smidd/NewLDMX/ldmx-sw/install/lib/libEvent.so");	

    main(options,args);			
	    
			
			
			
			
			
			
