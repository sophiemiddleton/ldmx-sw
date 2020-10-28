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
        for i in range(nentWAB):
            sum_e = 0
            sum_o = 0
            sum_e_rec = 0
            sum_o_rec = 0
            has_both = 0
            opening_cut = 0
            self.tin1.GetEntry(i);
            ElectronHCAL = Event()
            PhotonHCAL = Event()
            ElectronECAL = Event()
            PhotonECAL = Event()
            other = Event()
            incoming_elec_mom = []  
            incoming_phot_mom = []    
            incoming_elec_pos = []  
            incoming_phot_pos = []    
            self.Electron_Pos.append([])
            dr_cut = 10
            for ih,hit in enumerate(self.hcalScoringPlaneHits1):
                j = hit.getID() & 0xFFF;
                 
                if j==41 and hit.getTrackID()==1:
                    has_both +=1
                    ElectronHCAL.x.append(hit.getPosition()[0])
                    ElectronHCAL.y.append(hit.getPosition()[1])
                    ElectronHCAL.z.append(hit.getPosition()[2])
                    self.Elec_HCALSPHits_x.append(hit.getPosition()[0])
                    self.Elec_HCALSPHits_y.append(hit.getPosition()[1])
                    self.Elec_HCALSPHits_z.append(hit.getPosition()[2])
                    incoming_elec_mom = (hit.getMomentum()) 
                    incoming_elec_pos = (hit.getPosition()) 
                    self.Elec_HCALSPHits_e.append(hit.getEnergy())
                    sum_e +=hit.getEnergy()
                    self.Electron_Pos[i] = incoming_elec_pos
                    
                if j==41 and hit.getTrackID()==2:
                    has_both+=1
                    self.Phot_HCALSPHits_x.append(hit.getPosition()[0])
                    self.Phot_HCALSPHits_y.append(hit.getPosition()[1])
                    self.Phot_HCALSPHits_z.append(hit.getPosition()[2])
                    PhotonHCAL.x.append(hit.getPosition()[0])
                    PhotonHCAL.y.append(hit.getPosition()[1])
                    PhotonHCAL.z.append(hit.getPosition()[2]) 
                    incoming_phot_mom = hit.getMomentum()
                    incoming_phot_pos = hit.getPosition() 
                    self.Phot_HCALSPHits_e.append(hit.getEnergy())
                    if len(incoming_elec_pos) !=0:
                        dr = math.sqrt((hit.getPosition()[0]-incoming_elec_pos[0])*(hit.getPosition()[0]-incoming_elec_pos[0]) + (hit.getPosition()[1]-incoming_elec_pos[1])*(hit.getPosition()[1]-incoming_elec_pos[1]) +(hit.getPosition()[2]-incoming_elec_pos[2])*(hit.getPosition()[2]-incoming_elec_pos[2]))
                        if dr < dr_cut:
                            sum_e += hit.getEnergy()
                        if dr > dr_cut:
                            sum_o += hit.getEnergy()
                    
                if j==41 and hit.getTrackID()!=2 and hit.getTrackID()!=1:
                    other.x.append(hit.getPosition()[0])
                    other.y.append(hit.getPosition()[1])
                    other.z.append(hit.getPosition()[2])
                    """
                    if len(incoming_elec_pos) !=0:
                        dr = math.sqrt((hit.getPosition()[0]-incoming_elec_pos[0])*(hit.getPosition()[0]-incoming_elec_pos[0]) + (hit.getPosition()[1]-incoming_elec_pos[1])*(hit.getPosition()[1]-incoming_elec_pos[1]) +(hit.getPosition()[2]-incoming_elec_pos[2])*(hit.getPosition()[2]-incoming_elec_pos[2]))
                        if dr < dr_cut:
                            sum_e += hit.getEnergy()
                        if dr > dr_cut:
                            sum_o += hit.getEnergy()
                """
                if(len(incoming_phot_mom)!=0 and len(incoming_elec_mom)!=0):
                    
                    mag_g = math.sqrt(incoming_phot_mom[0]*incoming_phot_mom[0] + incoming_phot_mom[1]*incoming_phot_mom[1] +incoming_phot_mom[2]*incoming_phot_mom[2])
                    mag_e = math.sqrt(incoming_elec_mom[0]*incoming_elec_mom[0] + incoming_elec_mom[1]*incoming_elec_mom[1] +incoming_elec_mom[2]*incoming_elec_mom[2])
                    dot = (incoming_elec_mom[0]*incoming_phot_mom[0] + incoming_elec_mom[1]*incoming_phot_mom[1] +incoming_elec_mom[2]*incoming_phot_mom[2])
                    self.dr.append(math.sqrt((incoming_phot_pos[0]-incoming_elec_pos[0])*(incoming_phot_pos[0]-incoming_elec_pos[0]) + (incoming_phot_pos[1]-incoming_elec_pos[1])*(incoming_phot_pos[1]-incoming_elec_pos[1]) +(incoming_phot_pos[2]-incoming_elec_pos[2])*(incoming_phot_pos[2]-incoming_elec_pos[2])))
                    if(mag_g!=0 and mag_e!=0): 
                        opening_cut = math.acos(dot/(mag_g*mag_e))
                        self.open_angle.append(opening_cut)
                 
            self.Sum_HCALSPHits_e.append(sum_e)
            self.HCALEventList_elec.append(ElectronHCAL)
            self.HCALEventList_phot.append(PhotonHCAL)
            self.HCALEventList_other.append(other)
            self.ICR.append(sum_e)
            self.OCR.append(sum_o)
            
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
        """
        for i,j in enumerate(self.RecHitsEventList):
            electron_at = self.Electron_Pos[i]
            print(electron_at, j.x)
            if(len(j.x)!=0 and len(electron_at)!=0):
                for k,l in enumerate(j.x):
                    dr = math.sqrt((j.x[k]-electron_at[0])*(j.x[k]-electron_at[0]) + (j.y[k]-electron_at[1])*(j.y[k] - electron_at[1]) +(j.z[k]-electron_at[2])*(j.z[k]-electron_at[2]))
                    if dr < dr_cut:
                        sum_e_rec += hit.getEdep()
                    if dr > dr_cut:
                        sum_o_rec += hit.getEdep()
        """
        self.ICR_Rec.append(sum_e_rec)
        self.OCR_Rec.append(sum_o_rec)
            
def main(options,args) : 
    sc = WabEvent(options.ifile1, options.ofile,options.tag) ;
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.RecHitsESum,
            bins=50, 
            range=(0,1000), 
            label="Opening Angle")
   
    ax.set_xlabel('Total Energy [MeV] ')  
    plt.savefig("HCALTot.png")
    
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.open_angle,
            bins=50, 
            range=(0,math.pi), 
            label="Opening Angle")
   
    ax.set_xlabel('Opening Angle [rad] ')  
    plt.savefig("OpeningAngle.png")
    ig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.dr,
            bins=50, 
            range=(0,1600), 
            label="Distance")
   
    ax.set_xlabel('Distance between Electron and Photon at HCAL [mm]')  
    plt.savefig("dr.png")
   
    
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.Phot_HCALSPHits_x,
            bins=50, 
            range=(-800,800), 
            label="HCAL SP WAB Photons (x)")
   
    ax.set_xlabel('HCAL SP WAB Photons  x')  
    plt.savefig("HCALSPWABPhotons_x.png")
    
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.Elec_HCALSPHits_phi,
            bins=50, 
            range=(-4,4), 
            label="phi [rad]")
    
    ax.set_xlabel('phi [rad]')  
    plt.savefig("HCALElecPhi.png")
    

    
    fig, ax = plt.subplots(1,1)
    p = ax.scatter(sc.ICR_Rec, sc.OCR_Rec)
    ax.set_xlabel('ICR [MeV]')
    ax.set_ylabel('OCR [MeV]')  

    plt.savefig("ICRvOCR_Reco.png")
    
    
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.Phot_HCALSPHits_phi,
            bins=50, 
            range=(-4,4), 
            label="phi [rad]")
    
    ax.set_xlabel('phi [rad]')  
    plt.savefig("HCALPhotPhi.png")
    
    fig, ax = plt.subplots(1,1)
    p = ax.scatter(sc.ICR, sc.OCR)
    ax.set_xlabel('ICR [MeV]')
    ax.set_ylabel('OCR [MeV]')  

    plt.savefig("ICRvOCR.png")
    
    
    fig, ax = plt.subplots(1,1)
    p = ax.scatter(sc.Phot_HCALSPHits_x, sc.Phot_HCALSPHits_y, c=sc.Phot_HCALSPHits_e, cmap='inferno')
    ax.set_xlabel('HCAL x [mm]')
    ax.set_ylabel('HCAL y [mm]')  
    fig.colorbar(p)
    plt.xlim([-800,800])
    plt.ylim([-800,800])
    plt.savefig("ScoringPlanesXYPhotons.png")
    
    fig, ax = plt.subplots(1,1)
    p = ax.scatter(sc.Elec_HCALSPHits_x, sc.Elec_HCALSPHits_y, c=sc.Elec_HCALSPHits_e, cmap='inferno')
    ax.set_xlabel('HCAL x [mm]')
    ax.set_ylabel('HCAL y [mm]')  
    fig.colorbar(p)
    plt.xlim([-800,800])
    plt.ylim([-800,800])
    plt.savefig("ScoringPlanesXYElectrons.png")
    
    
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.Sum_HCALSPHits_e,
            bins=50, 
            range=(0,5000), 
            label="HCAL SP WAB Photons Energy [MeV]")
   
    ax.set_xlabel('Energy [MeV]')  
    plt.savefig("HCALSPWABTot_e.png")
    
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.Phot_HCALSPHits_e,
            bins=50, 
            range=(0,5000), 
            label="HCAL SP WAB Photons Energy [MeV]")
   
    ax.set_xlabel('Energy [MeV]')  
    plt.savefig("HCALSPWABPhotons_e.png")
    
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.Sum_HCALSPHits_e,
            bins=50, 
            range=(0,5000), 
            label="HCAL SP WAB Photons Energy [MeV]")
   
    ax.set_xlabel('Energy [MeV]')  
    plt.savefig("HCALSPWABTot_e.png")
    
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.Phot_ECALSPHits_e,
            bins=50, 
            range=(0,5000), 
            label="ECAL SP WAB Photons Energy [MeV]")
   
    ax.set_xlabel('Energy [MeV]')  
    plt.savefig("ECALSPWABPhotons_e.png")
    
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.Elec_ECALSPHits_e,
            bins=50, 
            range=(0,5000), 
            label="ECAL SP WAB Electron Energy [MeV]")
   
    ax.set_xlabel('Energy [MeV]')  
    plt.savefig("ECALSPWABElectrons_e.png")
    
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.Elec_HCALSPHits_e,
            bins=50, 
            range=(0,5000), 
            label="HCAL SP WAB Electrons Energy [MeV]")
   
    ax.set_xlabel('Energy [MeV]')  
    plt.savefig("HCALSPWABElectrons_e.png")
    
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
    for i, event in enumerate(sc.HCALEventList_elec):
        if i < 10:
            fig, ax = plt.subplots(1,1)
            p = ax.scatter(sc.HCALEventList_other[i].x, sc.HCALEventList_other[i].y, c='r', label='Secondary')
            p = ax.scatter(sc.HCALEventList_phot[i].x, sc.HCALEventList_phot[i].y,  c='g', label='Photon')
            p = ax.scatter(event.x, event.y,  c='b', label='Electron')
            ax.set_xlabel('HCAL Sim Hit x')
            ax.set_ylabel('HCAL Sim Hit y')
            fig.suptitle('HCAL Sim Hits Event '+str(i))  
            #fig.colorbar(p)
            ax.legend(prop={'size': 10})
            plt.savefig("test"+str(i)+".png")
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
    print("Norm", factor_wab)
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
	    
			
			
			
			
			
			
