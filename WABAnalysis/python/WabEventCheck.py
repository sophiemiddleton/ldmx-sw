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
        self.events = [];
        self.histograms = {};
        for i in range(10):
            #SimHits WABS:
            self.histograms["simhit_event_"+str(i)]= ROOT.TH3F("HCAL  SimHits Event_"+str(i),";HCAL SimHits Event"+str(i), 400, -800, 800, 400,0,800, 200, -800,800);
            
            #RecoHits WABS:
            self.histograms["rechit_event_"+str(i)]= ROOT.TH3F("HCAL  RecHits Event_"+str(i),";HCAL RecHits Event"+str(i), 400, -800, 800, 400,0,800, 200, -800,800);
            

        self.loop();
        self.writeOutHistos();

    def writeOutHistos(self):

        self.fout.cd();
        for key in self.histograms: 
            print(key);
            self.histograms[key].Write();
            #for ih, particle in enumerate(self.simParticles1):
                #idcode = particle.first     
    def loop(self):
        partList = []
       
        for i in range(10):
            GenData = Event();
            self.tin1.GetEntry(i);
            for ih, particle in enumerate(self.simParticles1):
                idcode = particle.first 
                if  particle.second.getProcessType() == 0:
                    partList = [idcode];
                if particle.second.getProcessType() ==0 and (particle.second.getPdgID()!=22 or particle.second.getPdgID() !=11):
                    print(i, ih, idcode, particle.second.getPdgID() )
        for i in range(10):
            GenData = Event()
            self.tin1.GetEntry(i);
            for ih,hit in enumerate(self.hcalSimHits1):
                [x, y, z] = hit.getPosition();
                GenData.x.append(x)
                GenData.y.append(y)
                GenData.z.append(z)
                GenData.e.append(hit.getEdep())
                self.histograms["simhit_event_"+str(i)].Fill(x,z,y,hit.getEdep());
            self.EventList.append(GenData)
            
        for i in range(10):
            GenData = Event()
            self.tin1.GetEntry(i);
            for ih,hit in enumerate(self.hcalHits1):
                 if (hit.isNoise()==0 and hit.getXPos()!=0 and hit.getYPos()!=0 and hit.getZPos()!=0):
                    GenData.x.append(hit.getXPos())
                    GenData.y.append(hit.getYPos())
                    GenData.z.append(hit.getZPos())
                    GenData.e.append(hit.getEnergy())
                    self.histograms["rechit_event_"+str(i)].Fill(hit.getXPos(), hit.getZPos(), hit.getYPos(), hit.getEnergy());
            self.RecHitList.append(GenData)
        for i in range(10):
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
            
        for i in range(10):
            Electron = Event()
            Photon = Event()
            self.tin1.GetEntry(i);
            for ih,hit in enumerate(self.hcalScoringPlaneHits1):
                j = hit.getID() & 0xFFF;
                if j==45:
                    print(j)
                if j==42 and  j==45 and  hit.getTrackID()==1:
                    Electron.x.append(hit.getPosition()[0])
                    Electron.y.append(hit.getPosition()[1])
                    Electron.z.append(hit.getPosition()[2])
                    Electron.e.append(hit.getEnergy())
                    Electron.pid.append(hit.getPdgID())
                    # Z= 240.00050354003906

                if j==42 and  j==45 and hit.getTrackID()==2:
                    Photon.x.append(hit.getPosition()[0])
                    Photon.y.append(hit.getPosition()[1])
                    Photon.z.append(hit.getPosition()[2])
                    Photon.e.append(hit.getEnergy())
                    Photon.pid.append(hit.getPdgID())
            self.HCALPlaneHitListElectron.append(Electron) 
            self.HCALPlaneHitListPhoton.append(Photon)  
            
        for i in range(10):
            Electron = Event()
            Photon = Event()
            self.tin1.GetEntry(i);
            for ih,hit in enumerate(self.ecalScoringPlaneHits1):
                j = hit.getID() & 0xFFF;
                if j==31 and hit.getTrackID()==1:
                    Electron.x.append(hit.getPosition()[0])
                    Electron.y.append(hit.getPosition()[1])
                    Electron.z.append(hit.getPosition()[2])
                    Electron.e.append(hit.getEnergy())
                    Electron.pid.append(hit.getPdgID())
            
                if j==31 and hit.getTrackID()==2:
                    Photon.x.append(hit.getPosition()[0])
                    Photon.y.append(hit.getPosition()[1])
                    Photon.z.append(hit.getPosition()[2])
                    Photon.e.append(hit.getEnergy())
                    Photon.pid.append(hit.getPdgID())
            self.ECALPlaneHitListElectron.append(Electron) 
            self.ECALPlaneHitListPhoton.append(Photon)
        
        
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
    for i, event in enumerate(sc.RecHitList):
        
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        primary_electron = ax.scatter3D(sc.HCALPlaneHitListElectron[i].x, sc.HCALPlaneHitListElectron[i].z, sc.HCALPlaneHitListElectron[i].y)
        primary_photon = ax.scatter3D(sc.HCALPlaneHitListPhoton[i].x, sc.HCALPlaneHitListPhoton[i].z, sc.HCALPlaneHitListPhoton[i].y)
        p = ax.scatter3D(event.x, event.z, event.y, c=event.e, cmap='inferno')
        
        #print("point at",event.x, event.z, event.y)
        #print(sc.HCALPlaneHitListElectron[i].x, sc.HCALPlaneHitListElectron[i].z, sc.HCALPlaneHitListElectron[i].y)
        ax.set_xlabel('HCAL Rec Hit x')
        ax.set_ylabel('HCAL Rec Hit z')
        ax.set_zlabel('HCAL Rec Hit y') 
        fig.suptitle('HCAL Reco Hits Event ' +str(i))   
        fig.colorbar(p)
        plt.xlim([-800,800])
        plt.ylim([-800,800])
        plt.savefig("RecHitEvent"+str(i)+".png")
        
    for i, event in enumerate(sc.EventList):
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        p = ax.scatter3D(event.x, event.z, event.y, c=event.e, cmap='inferno')
        ax.set_xlabel('HCAL Sim Hit x')
        ax.set_ylabel('HCAL Sim Hit z')
        ax.set_zlabel('HCAL Sim Hit y')    
        fig.suptitle('HCAL Sim Hits Event '+str(i))  
        fig.colorbar(p)
        plt.xlim([-800,800])
        plt.ylim([-800,800])
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
        
    for i, event in enumerate(sc.ECALPlaneHitListPhoton):
        fig = plt.figure()
        ax = plt.axes()
        n, bins, patches = ax.hist(event.x,
            bins=50, 
            range=(-800,800), 
            label="photon")
        ax.set_xlabel('ECAL Plane Hit x')
        #print(i, "photon", event.x)
        for k, elec in enumerate(sc.ECALPlaneHitListElectron):
            if k == i:
                n, bins, patches = ax.hist(elec.x,
                    bins=50, 
                    range=(-800,800), 
                    label="electron")
                #print(i,"electron", elec.x)
        ax.legend()
        fig.suptitle('ECAL Primary Photon at Scoring Plane Event '+str(i))
        plt.savefig("1D_ECALPlaneHitEvent_x"+str(i)+".png")
        
    for i, event in enumerate(sc.ECALPlaneHitListPhoton):
        fig = plt.figure()
        ax = plt.axes()
        n, bins, patches = ax.hist(event.y,
            bins=50, 
            range=(-800,800), 
            label="photon")
        ax.set_xlabel('ECAL Plane Hit y')
        for k, elec in enumerate(sc.ECALPlaneHitListElectron):
            if k == i:
                n, bins, patches = ax.hist(elec.y,
                    bins=50, 
                    range=(-800,800), 
                    label="electron")
        ax.legend()
        fig.suptitle('ECAL Primary Photon at Scoring Plane Event '+str(i))
        plt.savefig("1D_ECALPlaneHitEvent_y"+str(i)+".png")
        
    for i, event in enumerate(sc.ECALPlaneHitListPhoton):
        fig = plt.figure()
        ax = plt.axes()
        n, bins, patches = ax.hist(event.x,
            bins=50, 
            range=(-800,800), 
            label="photon")
        ax.set_xlabel('HCAL Plane Hit x')
        for k, elec in enumerate(sc.ECALPlaneHitListElectron):
            if k == i:
                n, bins, patches = ax.hist(elec.x,
                    bins=50, 
                    range=(-800,800), 
                    label="electron")
        ax.legend()
        fig.suptitle('ECAL Primary Photon at Scoring Plane Event '+str(i))
        plt.savefig("1D_ECALPlaneHitEvent_x"+str(i)+".png")
        
    for i, event in enumerate(sc.ECALPlaneHitListPhoton):
        fig = plt.figure()
        ax = plt.axes()
        n, bins, patches = ax.hist(event.y,
            bins=50, 
            range=(-800,800), 
            label="photon")
        ax.set_xlabel('ECAL Plane Hit y')
        for k, elec in enumerate(sc.ECALPlaneHitListElectron):
            if k == i:
                n, bins, patches = ax.hist(elec.y,
                    bins=50, 
                    range=(-800,800), 
                    label="electron")
        ax.legend()
        fig.suptitle('ECAL Primary Photon at Scoring Plane Event '+str(i))
        plt.savefig("1D_ECALPlaneHitEvent_y"+str(i)+".png")
        
          
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
	    
			
			
			
			
			
			
