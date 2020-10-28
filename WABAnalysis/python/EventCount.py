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
        self.HCALEventList_elec = []
        self.HCALEventList_other = []
        self.HCALEventList_phot = []
        self.ECALEventList_elec = []
        self.ECALEventList_phot = []
        self.RecHitsEventList = []
        self.HCALPlaneHitListElectron = []
        self.HCALPlaneHitListPhoton = []
        self.ECALPlaneHitListElectron = []
        self.ECALPlaneHitListPhoton = []
        self.SumHCAL = []
        self.RecHitsESum = []
        #wab:
        self.evHeader1 = ROOT.ldmx.EventHeader()
        self.simParticles1 = ROOT.std.map(int,'ldmx::SimParticle')();
        self.hcalHits1 = ROOT.std.vector('ldmx::HcalHit')();
        self.hcalSimHits1 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.recoilSimHits1 = ROOT.std.vector('ldmx::SimTrackerHit')();
        self.hcalScoringPlaneHits1 = ROOT.std.vector('ldmx::SimTrackerHit')();
        self.ecalScoringPlaneHits1 = ROOT.std.vector('ldmx::SimTrackerHit')();
        self.targetScoringPlaneHits1 = ROOT.std.vector('ldmx::SimTrackerHit')();
        self.tin1.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader1 ));
        self.tin1.SetBranchAddress("HcalSimHits_v12",  ROOT.AddressOf( self.hcalSimHits1 ));
        self.tin1.SetBranchAddress("RecoilSimHits_v12",  ROOT.AddressOf( self.recoilSimHits1 ));
        self.tin1.SetBranchAddress("HcalRecHits_v12",  ROOT.AddressOf( self.hcalHits1 ));
        self.tin1.SetBranchAddress("HcalScoringPlaneHits_v12",  ROOT.AddressOf(self.hcalScoringPlaneHits1));
        self.tin1.SetBranchAddress("EcalScoringPlaneHits_v12",  ROOT.AddressOf(self.ecalScoringPlaneHits1));
        self.tin1.SetBranchAddress("TargetScoringPlaneHits_v12",  ROOT.AddressOf(self.targetScoringPlaneHits1));
        self.tin1.SetBranchAddress("SimParticles_v12",  ROOT.AddressOf( self.simParticles1 ));
        self.RecHits_x = []
        self.RecHits_y = []
        self.RecHits_z = []
        
        self.open_angle = []
        self.dr = []
        self.Electron_Pos = []
        
        self.Elec_HCALSPHits_x = []
        self.Elec_HCALSPHits_y = []
        self.Elec_HCALSPHits_z = []
        self.Elec_HCALSPHits_e = []
        self.Elec_HCALSPHits_phi = []
        
        self.ICR = []
        self.OCR = []
        self.ICR_Rec = []
        self.OCR_Rec = []
        self.Elec_ECALSPHits_x = []
        self.Elec_ECALSPHits_y = []
        self.Elec_ECALSPHits_z = []
        self.Elec_ECALSPHits_e = []
        
        self.Phot_HCALSPHits_x = []
        self.Phot_HCALSPHits_y = []
        self.Phot_HCALSPHits_z = []
        self.Phot_HCALSPHits_e = []
        self.Phot_HCALSPHits_phi = []
        
        self.Sum_HCALSPHits_e = []
        
        self.Phot_ECALSPHits_x = []
        self.Phot_ECALSPHits_y = []
        self.Phot_ECALSPHits_z = []
        self.Phot_ECALSPHits_e = []
        
        self.loop();
        self.writeOutHistos();

    def writeOutHistos(self):

        self.fout.cd();
        
    def loop(self):
 
        nentWAB = self.tin1.GetEntriesFast();
        for i in range(nentWAB):
            sum_sim_e = 0
            self.tin1.GetEntry(i);
            for ih,hit in enumerate(self.hcalSimHits1):
                sum_sim_e += hit.getEdep();
            self.SumHCAL.append(sum_sim_e)
        
        for i in range(nentWAB):
            event = Event()
            self.tin1.GetEntry(i);
            sum_sim_e = 0
            for ih,hit in enumerate(self.hcalHits1):
                 if (hit.isNoise()==0 and hit.getXPos()!=0 and hit.getYPos()!=0 and hit.getZPos()!=0):
                    event.x.append(hit.getXPos())
                    event.y.append(hit.getYPos())
                    event.z.append(hit.getZPos())
                    sum_sim_e += hit.getEnergy()
            self.RecHitsEventList.append(event)    
            self.RecHitsESum.append(sum_sim_e)
        
        k = 300
        for k in range(390,410):
            has_electron_at_HCAL = 0
            has_photon_at_HCAL = 0
            has_both_at_HCAL = 0
            has_electron_in_side = 0
            has_photon_in_side = 0
            has_electron_and_photon = 0
            has_both_at_HCAL_side = 0
            has_electron_and_photon_side = 0
            if k!=400 and k%10 !=0:
                continue
            for i in range(nentWAB):
                has_both_at_HCAL = 0
                has_both_at_HCAL_side=0
                self.tin1.GetEntry(i);
                ElectronHCAL = Event()
                PhotonHCAL = Event()
                ElectronECAL = Event()
                PhotonECAL = Event()
                other = Event()

                for ih,hit in enumerate(self.hcalScoringPlaneHits1):
                    j = hit.getID() & 0xFFF;
                     
                    if j==41 and hit.getTrackID()==1:
                        has_both_at_HCAL +=1
                        has_electron_at_HCAL +=1
                        if abs(hit.getPosition()[0]) > k or abs(hit.getPosition()[1]) > 300:
                            has_electron_in_side +=1
                            has_both_at_HCAL_side +=1
                            ElectronHCAL.x.append(hit.getPosition()[0])
                            ElectronHCAL.y.append(hit.getPosition()[1])
                            ElectronHCAL.z.append(hit.getPosition()[2])
                            self.Elec_HCALSPHits_x.append(hit.getPosition()[0])
                            self.Elec_HCALSPHits_y.append(hit.getPosition()[1])
                            self.Elec_HCALSPHits_z.append(hit.getPosition()[2])
                        
                    if j==41 and hit.getTrackID()==2:
                        has_both_at_HCAL+=1
                        has_photon_at_HCAL +=1
                        
                        if abs(hit.getPosition()[0]) > k or abs(hit.getPosition()[1]) > 300:
                            has_photon_in_side +=1
                            has_both_at_HCAL_side +=1
                            self.Phot_HCALSPHits_x.append(hit.getPosition()[0])
                            self.Phot_HCALSPHits_y.append(hit.getPosition()[1])
                            self.Phot_HCALSPHits_z.append(hit.getPosition()[2])
                            PhotonHCAL.x.append(hit.getPosition()[0])
                            PhotonHCAL.y.append(hit.getPosition()[1])
                            PhotonHCAL.z.append(hit.getPosition()[2]) 
                                
                        
                    if j==41 and hit.getTrackID()!=2 and hit.getTrackID()!=1:
                        if abs(hit.getPosition()[0]) > k or abs(hit.getPosition()[1]) > 300:
                                
                            other.x.append(hit.getPosition()[0])
                            other.y.append(hit.getPosition()[1])
                            other.z.append(hit.getPosition()[2])
                  
                if has_both_at_HCAL ==2:
                    has_electron_and_photon +=1
                if has_both_at_HCAL_side==2:
                    has_electron_and_photon_side +=1
                self.HCALEventList_elec.append(ElectronHCAL)
                self.HCALEventList_phot.append(PhotonHCAL)
                self.HCALEventList_other.append(other)
            print(k,has_electron_in_side,has_photon_in_side,has_electron_and_photon_side)
        
        #print("Has Both", has_electron_and_photon, "Has Electron at HCAL",has_electron_at_HCAL, "Has Photon", has_photon_at_HCAL , "Has E at Side", has_electron_in_side,"Has Photon at side", has_photon_in_side, "both in side", has_electron_and_photon_side)
        for i in range(nentWAB):
            ElectronECAL = Event()
            PhotonECAL = Event()
            self.tin1.GetEntry(i);
            for ih,hit in enumerate(self.ecalScoringPlaneHits1):
                j = hit.getID() & 0xFFF;
               
                if j==34 and hit.getTrackID()==1:
                    self.Elec_ECALSPHits_x.append(hit.getPosition()[0])
                    self.Elec_ECALSPHits_y.append(hit.getPosition()[1])
                    self.Elec_ECALSPHits_z.append(hit.getPosition()[2])
                    self.Elec_ECALSPHits_e.append(hit.getEnergy())
                    ElectronECAL.x.append(hit.getPosition()[0])
                    ElectronECAL.y.append(hit.getPosition()[1])
                    ElectronECAL.z.append(hit.getPosition()[2]) 
                if j==34 and hit.getTrackID()==2:
                    self.Phot_ECALSPHits_x.append(hit.getPosition()[0])
                    self.Phot_ECALSPHits_y.append(hit.getPosition()[1])
                    self.Phot_ECALSPHits_z.append(hit.getPosition()[2])
                    self.Phot_ECALSPHits_e.append(hit.getEnergy())
                    PhotonECAL.x.append(hit.getPosition()[0])
                    PhotonECAL.y.append(hit.getPosition()[1])
                    PhotonECAL.z.append(hit.getPosition()[2]) 
            self.ECALEventList_elec.append(ElectronECAL)
            self.ECALEventList_phot.append(PhotonECAL)

def main(options,args) : 
    sc = WabEvent(options.ifile1, options.ofile,options.tag) ;
    
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
	    
			
			
			
			
			
			
