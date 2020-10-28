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


lumin = 0.22
proton_wabs = 74.
n_wab = 10000.

#FF1
xsec_FF1 = 2811000
event_weight_FF1 = xsec_FF1/n_wab
factor_wab_FF1 = event_weight_FF1 * 1 * lumin

#FF2
xsec_FF2 = 2.394
event_weight_FF2 = xsec_FF2/n_wab
factor_wab_FF2 = event_weight_FF2 * 1 * lumin

#FF3
xsec_FF3 = 1795.
event_weight_FF3 = xsec_FF3/n_wab
factor_wab_FF3 = event_weight_FF3 * proton_wabs * lumin

#FF4
xsec_FF4 = 15.83
event_weight_FF4 = xsec_FF4/n_wab
factor_wab_FF4 = event_weight_FF4 * proton_wabs * lumin

class Event():
    def __init__(self):
        self.x = []
        self.y = []
        self.z = []
        self.e = []
        self.pid = []
        
class WabEvent:
   
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
        self.hcalHits1 = ROOT.std.vector('ldmx::HcalHit')();
        self.hcalSimHits1 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.hcalScoringPlaneHits1 = ROOT.std.vector('ldmx::SimTrackerHit')();
        self.ecalScoringPlaneHits1 = ROOT.std.vector('ldmx::SimTrackerHit')();

        self.tin1.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader1 ));
        self.tin1.SetBranchAddress("HcalSimHits_v12",  ROOT.AddressOf( self.hcalSimHits1 ));
        self.tin1.SetBranchAddress("HcalRecHits_v12",  ROOT.AddressOf( self.hcalHits1 ));
        self.tin1.SetBranchAddress("HcalScoringPlaneHits_v12",  ROOT.AddressOf(self.hcalScoringPlaneHits1));
        self.tin1.SetBranchAddress("EcalScoringPlaneHits_v12",  ROOT.AddressOf(self.ecalScoringPlaneHits1));
        
         #FF2:
        self.evHeader2 = ROOT.ldmx.EventHeader()
        self.hcalHits2 = ROOT.std.vector('ldmx::HcalHit')();
        self.hcalSimHits2 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.hcalScoringPlaneHits2 = ROOT.std.vector('ldmx::SimTrackerHit')();
        self.ecalScoringPlaneHits2 = ROOT.std.vector('ldmx::SimTrackerHit')();

        self.tin2.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader2));
        self.tin2.SetBranchAddress("HcalSimHits_v12",  ROOT.AddressOf( self.hcalSimHits2 ));
        self.tin2.SetBranchAddress("HcalRecHits_v12",  ROOT.AddressOf( self.hcalHits2));
        self.tin2.SetBranchAddress("HcalScoringPlaneHits_v12",  ROOT.AddressOf(self.hcalScoringPlaneHits2));
        self.tin2.SetBranchAddress("EcalScoringPlaneHits_v12",  ROOT.AddressOf(self.ecalScoringPlaneHits2));
        
        
         #FF3:
        self.evHeader3 = ROOT.ldmx.EventHeader()
        self.hcalHits3 = ROOT.std.vector('ldmx::HcalHit')();
        self.hcalSimHits3 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.hcalScoringPlaneHits3 = ROOT.std.vector('ldmx::SimTrackerHit')();
        self.ecalScoringPlaneHits3 = ROOT.std.vector('ldmx::SimTrackerHit')();

        self.tin3.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader3 ));
        self.tin3.SetBranchAddress("HcalSimHits_v12",  ROOT.AddressOf( self.hcalSimHits3 ));
        self.tin3.SetBranchAddress("HcalRecHits_v12",  ROOT.AddressOf( self.hcalHits3 ));
        self.tin3.SetBranchAddress("HcalScoringPlaneHits_v12",  ROOT.AddressOf(self.hcalScoringPlaneHits3));
        self.tin3.SetBranchAddress("EcalScoringPlaneHits_v12",  ROOT.AddressOf(self.ecalScoringPlaneHits3));
        
         #FF4:
        self.evHeader4 = ROOT.ldmx.EventHeader()
        self.hcalHits4 = ROOT.std.vector('ldmx::HcalHit')();
        self.hcalSimHits4 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.hcalScoringPlaneHits4 = ROOT.std.vector('ldmx::SimTrackerHit')();
        self.ecalScoringPlaneHits4 = ROOT.std.vector('ldmx::SimTrackerHit')();

        self.tin4.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader4 ));
        self.tin4.SetBranchAddress("HcalSimHits_v12",  ROOT.AddressOf( self.hcalSimHits4 ));
        self.tin4.SetBranchAddress("HcalRecHits_v12",  ROOT.AddressOf( self.hcalHits4 ));
        self.tin4.SetBranchAddress("HcalScoringPlaneHits_v12",  ROOT.AddressOf(self.hcalScoringPlaneHits4));
        self.tin4.SetBranchAddress("EcalScoringPlaneHits_v12",  ROOT.AddressOf(self.ecalScoringPlaneHits4));

      
        self.RecHits_x = []
        self.RecHits_y = []
        self.RecHits_z = []
        
        self.SimHits_x = []
        self.SimHits_y = []
        self.SimHits_z = []
        # FF1:
        self.Elec_HCALSPHits_x_FF1 = []
        self.Elec_HCALSPHits_y_FF1  = []
        self.Elec_HCAL_weights_FF1  = []
        self.Elec_ECALSPHits_x_FF1  = []
        self.Elec_ECALSPHits_y_FF1  = []
        self.Elec_ECAL_weights_FF1  = []
        
        self.Phot_HCALSPHits_x_FF1  = []
        self.Phot_HCALSPHits_y_FF1  = []
        self.Phot_HCAL_weights_FF1  = []
        self.Phot_ECALSPHits_x_FF1  = []
        self.Phot_ECALSPHits_y_FF1  = []
        self.Phot_ECAL_weights_FF1  = []
        
         # FF2:
        self.Elec_HCALSPHits_x_FF2 = []
        self.Elec_HCALSPHits_y_FF2  = []
        self.Elec_HCAL_weights_FF2  = []
        self.Elec_ECALSPHits_x_FF2  = []
        self.Elec_ECALSPHits_y_FF2  = []
        self.Elec_ECAL_weights_FF2  = []
        
        self.Phot_HCALSPHits_x_FF2  = []
        self.Phot_HCALSPHits_y_FF2  = []
        self.Phot_HCAL_weights_FF2  = []
        self.Phot_ECALSPHits_x_FF2  = []
        self.Phot_ECALSPHits_y_FF2  = []
        self.Phot_ECAL_weights_FF2  = []
        
         # FF3:
        self.Elec_HCALSPHits_x_FF3 = []
        self.Elec_HCALSPHits_y_FF3  = []
        self.Elec_HCAL_weights_FF3  = []
        self.Elec_ECALSPHits_x_FF3  = []
        self.Elec_ECALSPHits_y_FF3  = []
        self.Elec_ECAL_weights_FF3  = []
        
        self.Phot_HCALSPHits_x_FF3  = []
        self.Phot_HCALSPHits_y_FF3  = []
        self.Phot_HCAL_weights_FF3  = []
        self.Phot_ECALSPHits_x_FF3  = []
        self.Phot_ECALSPHits_y_FF3  = []
        self.Phot_ECAL_weights_FF3  = []
        
         # FF4:
        self.Elec_HCALSPHits_x_FF4 = []
        self.Elec_HCALSPHits_y_FF4  = []
        self.Elec_HCAL_weights_FF4  = []
        self.Elec_ECALSPHits_x_FF4  = []
        self.Elec_ECALSPHits_y_FF4  = []
        self.Elec_ECAL_weights_FF4  = []
        
        self.Phot_HCALSPHits_x_FF4  = []
        self.Phot_HCALSPHits_y_FF4  = []
        self.Phot_HCAL_weights_FF4  = []
        self.Phot_ECALSPHits_x_FF4  = []
        self.Phot_ECALSPHits_y_FF4  = []
        self.Phot_ECAL_weights_FF4  = []
        self.loop();
        self.writeOutHistos();

    def writeOutHistos(self):

        self.fout.cd();
        
    def loop(self):
        partList = []
        nentWAB1 = self.tin1.GetEntriesFast();
        nentWAB2 = self.tin2.GetEntriesFast();
        nentWAB3 = self.tin3.GetEntriesFast();
        nentWAB4 = self.tin4.GetEntriesFast();
        
        for i in range(nentWAB1):
          
            self.tin1.GetEntry(i);
            for ih,hit in enumerate(self.hcalSimHits1):
                [x, y, z] = hit.getPosition();
                self.SimHits_x.append(x)
                self.SimHits_y.append(y)
                self.SimHits_z.append(z)
            
        for i in range(nentWAB1):
            
            self.tin1.GetEntry(i);
            for ih,hit in enumerate(self.hcalHits1):
                 if (hit.isNoise()==0 and hit.getXPos()!=0 and hit.getYPos()!=0 and hit.getZPos()!=0):
                    self.RecHits_x.append(hit.getXPos())
                    self.RecHits_y.append(hit.getYPos())
                    self.RecHits_z.append(hit.getZPos())
                    
        for i in range(nentWAB1):
           
            self.tin1.GetEntry(i);
            for ih,hit in enumerate(self.hcalScoringPlaneHits1):

                if hit.getTrackID()==1:
                    self.Elec_HCALSPHits_x_FF1.append(hit.getPosition()[0])
                    self.Elec_HCALSPHits_y_FF1.append(hit.getPosition()[1])
                 
                    self.Elec_HCAL_weights_FF1.append(factor_wab_FF1)
                if hit.getTrackID()==2:
                    self.Phot_HCALSPHits_x_FF1.append(hit.getPosition()[0])
                    self.Phot_HCALSPHits_y_FF1.append(hit.getPosition()[1])
                    self.Phot_HCAL_weights_FF1.append(factor_wab_FF1)
            
        for i in range(nentWAB1):
           
            self.tin1.GetEntry(i);
            for ih,hit in enumerate(self.ecalScoringPlaneHits1):
                if hit.getTrackID()==1:
                    self.Elec_ECALSPHits_x_FF1.append(hit.getPosition()[0])
                    self.Elec_ECALSPHits_y_FF1.append(hit.getPosition()[1])
                    self.Elec_ECAL_weights_FF1.append(factor_wab_FF1)
                if hit.getTrackID()==2:
                    self.Phot_ECALSPHits_x_FF1.append(hit.getPosition()[0])
                    self.Phot_ECALSPHits_y_FF1.append(hit.getPosition()[1])
                    self.Phot_ECAL_weights_FF1.append(factor_wab_FF1)
                    
                    
                    
        for i in range(nentWAB2):
           
            self.tin2.GetEntry(i);
            for ih,hit in enumerate(self.hcalScoringPlaneHits2):
                
                if hit.getTrackID()==1:
                    self.Elec_HCALSPHits_x_FF2.append(hit.getPosition()[0])
                    self.Elec_HCALSPHits_y_FF2.append(hit.getPosition()[1])
                 
                    self.Elec_HCAL_weights_FF2.append(factor_wab_FF2)
                if hit.getTrackID()==2:
                    self.Phot_HCALSPHits_x_FF2.append(hit.getPosition()[0])
                    self.Phot_HCALSPHits_y_FF2.append(hit.getPosition()[1])
                    self.Phot_HCAL_weights_FF2.append(factor_wab_FF2)
            
        for i in range(nentWAB2):
           
            self.tin2.GetEntry(i);
            for ih,hit in enumerate(self.ecalScoringPlaneHits2):
                if hit.getTrackID()==1:
                    self.Elec_ECALSPHits_x_FF2.append(hit.getPosition()[0])
                    self.Elec_ECALSPHits_y_FF2.append(hit.getPosition()[1])
                    self.Elec_ECAL_weights_FF2.append(factor_wab_FF2)
                if hit.getTrackID()==2:
                    self.Phot_ECALSPHits_x_FF2.append(hit.getPosition()[0])
                    self.Phot_ECALSPHits_y_FF2.append(hit.getPosition()[1])
                    self.Phot_ECAL_weights_FF2.append(factor_wab_FF2)
            
            
        for i in range(nentWAB3):
           
            self.tin3.GetEntry(i);
            for ih,hit in enumerate(self.hcalScoringPlaneHits3):
                
                if hit.getTrackID()==1:
                    self.Elec_HCALSPHits_x_FF3.append(hit.getPosition()[0])
                    self.Elec_HCALSPHits_y_FF3.append(hit.getPosition()[1])
                 
                    self.Elec_HCAL_weights_FF3.append(factor_wab_FF3)
                if hit.getTrackID()==2:
                    self.Phot_HCALSPHits_x_FF3.append(hit.getPosition()[0])
                    self.Phot_HCALSPHits_y_FF3.append(hit.getPosition()[1])
                    self.Phot_HCAL_weights_FF3.append(factor_wab_FF3)
            
        for i in range(nentWAB3):
           
            self.tin3.GetEntry(i);
            for ih,hit in enumerate(self.ecalScoringPlaneHits3):
                if hit.getTrackID()==1:
                    self.Elec_ECALSPHits_x_FF3.append(hit.getPosition()[0])
                    self.Elec_ECALSPHits_y_FF3.append(hit.getPosition()[1])
                    self.Elec_ECAL_weights_FF3.append(factor_wab_FF3)
                if hit.getTrackID()==2:
                    self.Phot_ECALSPHits_x_FF3.append(hit.getPosition()[0])
                    self.Phot_ECALSPHits_y_FF3.append(hit.getPosition()[1])
                    self.Phot_ECAL_weights_FF3.append(factor_wab_FF3)
                    
                    
                    
        for i in range(nentWAB4):
           
            self.tin4.GetEntry(i);
            for ih,hit in enumerate(self.hcalScoringPlaneHits4):
                
                if hit.getTrackID()==1:
                    self.Elec_HCALSPHits_x_FF4.append(hit.getPosition()[0])
                    self.Elec_HCALSPHits_y_FF4.append(hit.getPosition()[1])
                    self.Elec_HCAL_weights_FF4.append(factor_wab_FF4)
                if hit.getTrackID()==2:
                    self.Phot_HCALSPHits_x_FF4.append(hit.getPosition()[0])
                    self.Phot_HCALSPHits_y_FF4.append(hit.getPosition()[1])
                    self.Phot_HCAL_weights_FF4.append(factor_wab_FF4)
            
        for i in range(nentWAB4):
           
            self.tin4.GetEntry(i);
            for ih,hit in enumerate(self.ecalScoringPlaneHits4):
                if hit.getTrackID()==1:
                    self.Elec_ECALSPHits_x_FF4.append(hit.getPosition()[0])
                    self.Elec_ECALSPHits_y_FF4.append(hit.getPosition()[1])
                    self.Elec_ECAL_weights_FF4.append(factor_wab_FF4)
                if hit.getTrackID()==2:
                    self.Phot_ECALSPHits_x_FF4.append(hit.getPosition()[0])
                    self.Phot_ECALSPHits_y_FF4.append(hit.getPosition()[1])
                    self.Phot_ECAL_weights_FF4.append(factor_wab_FF4)
                    
		
def main(options,args) : 
    sc = WabEvent(options.ifile1, options.ifile2, options.ifile3, options.ifile4, options.ofile,options.tag) ;
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

    n, bins, patches = ax.hist(sc.Phot_HCALSPHits_x_FF1,
            bins= 50,histtype='step',
            range=(-800,800), weights=sc.Phot_HCAL_weights_FF1,
            label="FF1")
    n, bins, patches = ax.hist(sc.Phot_HCALSPHits_x_FF2,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Phot_HCAL_weights_FF2,
            label="FF2")
    n, bins, patches = ax.hist(sc.Phot_HCALSPHits_x_FF3,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Phot_HCAL_weights_FF3,
            label="FF3")
    n, bins, patches = ax.hist(sc.Phot_HCALSPHits_x_FF4,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Phot_HCAL_weights_FF4,
            label="FF4")
    ax.set_yscale('log')
    ax.legend()
    ax.set_xlabel('HCAL SP Photons x [mm]')  
    plt.savefig("HCALSPWABPhotons_x.png")
    
    fig, ax = plt.subplots(1,1)

    n, bins, patches = ax.hist(sc.Phot_ECALSPHits_x_FF1,
            bins= 50,histtype='step',
            range=(-800,800), weights=sc.Phot_ECAL_weights_FF1,
            label="FF1")
    n, bins, patches = ax.hist(sc.Phot_ECALSPHits_x_FF2,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Phot_ECAL_weights_FF2,
            label="FF2")
    n, bins, patches = ax.hist(sc.Phot_ECALSPHits_x_FF3,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Phot_ECAL_weights_FF3,
            label="FF3")
    n, bins, patches = ax.hist(sc.Phot_ECALSPHits_x_FF4,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Phot_ECAL_weights_FF4,
            label="FF4")
    ax.set_yscale('log')
    ax.legend()
    ax.set_xlabel('ECAL SP Photons x [mm]')  
    plt.savefig("ECALSPWABPhotons_x.png")
    
    fig, ax = plt.subplots(1,1)

    n, bins, patches = ax.hist(sc.Elec_HCALSPHits_x_FF1,
            bins= 50,histtype='step',
            range=(-800,800), weights=sc.Elec_HCAL_weights_FF1,
            label="FF1")
    n, bins, patches = ax.hist(sc.Elec_HCALSPHits_x_FF2,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Elec_HCAL_weights_FF2,
            label="FF2")
    n, bins, patches = ax.hist(sc.Elec_HCALSPHits_x_FF3,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Elec_HCAL_weights_FF3,
            label="FF3")
    n, bins, patches = ax.hist(sc.Elec_HCALSPHits_x_FF4,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Elec_HCAL_weights_FF4,
            label="FF4")
    ax.set_yscale('log')
    ax.legend()
    ax.set_xlabel('HCAL SP Electrons x [mm]')  
    plt.savefig("HCALSPWABElectrons_x.png")
    
    fig, ax = plt.subplots(1,1)

    n, bins, patches = ax.hist(sc.Elec_ECALSPHits_x_FF1,
            bins= 50,histtype='step',
            range=(-800,800), weights=sc.Elec_ECAL_weights_FF1,
            label="FF1")
    n, bins, patches = ax.hist(sc.Elec_ECALSPHits_x_FF2,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Elec_ECAL_weights_FF2,
            label="FF2")
    n, bins, patches = ax.hist(sc.Elec_ECALSPHits_x_FF3,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Elec_ECAL_weights_FF3,
            label="FF3")
    n, bins, patches = ax.hist(sc.Elec_ECALSPHits_x_FF4,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Elec_ECAL_weights_FF4,
            label="FF4")
    ax.set_yscale('log')
    ax.legend()
    ax.set_xlabel('ECAL SP Electrons x [mm]')  
    plt.savefig("ECALSPWABElectrons_x.png")
    
    fig, ax = plt.subplots(1,1)

    n, bins, patches = ax.hist(sc.Phot_HCALSPHits_y_FF1,
            bins= 50,histtype='step',
            range=(-800,800), weights=sc.Phot_HCAL_weights_FF1,
            label="FF1")
    n, bins, patches = ax.hist(sc.Phot_HCALSPHits_y_FF2,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Phot_HCAL_weights_FF2,
            label="FF2")
    n, bins, patches = ax.hist(sc.Phot_HCALSPHits_y_FF3,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Phot_HCAL_weights_FF3,
            label="FF3")
    n, bins, patches = ax.hist(sc.Phot_HCALSPHits_y_FF4,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Phot_HCAL_weights_FF4,
            label="FF4")
    ax.set_yscale('log')
    ax.legend()
    ax.set_xlabel('HCAL SP Photons y [mm]')  
    plt.savefig("HCALSPWABPhotons_y.png")
    
    fig, ax = plt.subplots(1,1)

    n, bins, patches = ax.hist(sc.Phot_ECALSPHits_y_FF1,
            bins= 50,histtype='step',
            range=(-800,800), weights=sc.Phot_ECAL_weights_FF1,
            label="FF1")
    n, bins, patches = ax.hist(sc.Phot_ECALSPHits_y_FF2,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Phot_ECAL_weights_FF2,
            label="FF2")
    n, bins, patches = ax.hist(sc.Phot_ECALSPHits_y_FF3,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Phot_ECAL_weights_FF3,
            label="FF3")
    n, bins, patches = ax.hist(sc.Phot_ECALSPHits_y_FF4,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Phot_ECAL_weights_FF4,
            label="FF4")
    ax.set_yscale('log')
    ax.legend()
    ax.set_xlabel('ECAL SP Photons y [mm]')  
    plt.savefig("ECALSPWABPhotons_y.png")
    
    fig, ax = plt.subplots(1,1)
   
    n, bins, patches = ax.hist(sc.Elec_HCALSPHits_y_FF1,
            bins= 50,histtype='step',
            range=(-800,800), weights=sc.Elec_HCAL_weights_FF1,
            label="FF1")
    n, bins, patches = ax.hist(sc.Elec_HCALSPHits_y_FF2,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Elec_HCAL_weights_FF2,
            label="FF2")
    n, bins, patches = ax.hist(sc.Elec_HCALSPHits_y_FF3,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Elec_HCAL_weights_FF3,
            label="FF3")
    n, bins, patches = ax.hist(sc.Elec_HCALSPHits_y_FF4,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Elec_HCAL_weights_FF4,
            label="FF4")
    ax.set_yscale('log')
    ax.legend()
    ax.set_xlabel('HCAL SP Electrons y [mm]')  
    plt.savefig("HCALSPWABElectrons_y.png")
    
    fig, ax = plt.subplots(1,1)

    n, bins, patches = ax.hist(sc.Elec_ECALSPHits_y_FF1,
            bins= 50,histtype='step',
            range=(-800,800), weights=sc.Elec_ECAL_weights_FF1,
            label="FF1")
    n, bins, patches = ax.hist(sc.Elec_ECALSPHits_y_FF2,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Elec_ECAL_weights_FF2,
            label="FF2")
    n, bins, patches = ax.hist(sc.Elec_ECALSPHits_y_FF3,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Elec_ECAL_weights_FF3,
            label="FF3")
    n, bins, patches = ax.hist(sc.Elec_ECALSPHits_y_FF4,
            bins=50, alpha=0.5,histtype='step',
            range=(-800,800), weights=sc.Elec_ECAL_weights_FF4,
            label="FF4")
    ax.set_yscale('log')
    ax.legend()
    ax.set_xlabel('ECAL SP Electrons y [mm]')  
    plt.savefig("ECALSPWABElectrons_y.png")
      
    sc.fout.Close();
    
if __name__ == "__main__":
	
    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    # input data files (4)
    parser.add_option('-a','--ifile1', dest='ifile1', default = 'file1.root',help='directory with data1', metavar='idir1')
    parser.add_option('-e','--ifile2', dest='ifile2', default = 'file2.root',help='directory with data2', metavar='idir2')
    parser.add_option('-c','--ifile3', dest='ifile3', default = 'file3.root',help='directory with data3', metavar='idir3')
    parser.add_option('-f','--ifile4', dest='ifile4', default = 'file4.root',help='directory with data4', metavar='idir4')
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
	    
			
			
			
			
			
			
