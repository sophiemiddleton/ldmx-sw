#!/usr/bin/python
import argparse
import importlib
import ROOT
from ROOT import TTree, TBranch
ROOT.gSystem.Load("/nfs/slac/g/ldmx/users/smidd/NewLDMX/ldmx-sw/install/lib/libEvent.so")
import os
import math
import sys
from array import array
from optparse import OptionParser
import matplotlib.pyplot as plt
sys.path.insert(0, '../')

lumin = 0.22
#wabs
proton_wabs = 74.
xsec_wab = 1795.
n_wab = 10000.
lumin_wab = 0.22
event_weight_wab = xsec_wab/n_wab
factor_wab = event_weight_wab * proton_wabs * lumin

class WabSigCompare:
   
    def __init__(self, fn1, ofn, tag):

        #wab:
        self.fin1 = ROOT.TFile(fn1);
        self.tin1 = self.fin1.Get("LDMX_Events")
     
        self.tag = int(tag);

        self.fn_out = ofn;
        self.fout = ROOT.TFile("hist_"+self.fn_out,"RECREATE");

        #wab:
        self.evHeader1 = ROOT.ldmx.EventHeader()
        self.simParticles1 = ROOT.std.map(int, 'ldmx::SimParticle')();
        self.hcalHits1 = ROOT.std.vector('ldmx::HcalHit')();
        self.hcalSimHits1 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.ecalHits1 = ROOT.std.vector('ldmx::EcalHit')();
        self.ecalSimHits1 = ROOT.std.vector('ldmx::SimCalorimeterHit')();
        self.hcalScoringPlaneHits1 = ROOT.std.vector('ldmx::SimTrackerHit')();
        self.ecalScoringPlaneHits1 = ROOT.std.vector('ldmx::SimTrackerHit')();
        self.targetScoringPlaneHits1 = ROOT.std.vector('ldmx::SimTrackerHit')();
        self.tin1.SetBranchAddress("EventHeader",  ROOT.AddressOf( self.evHeader1 ));
        self.tin1.SetBranchAddress("HcalSimHits_v12",  ROOT.AddressOf( self.hcalSimHits1 ));
        self.tin1.SetBranchAddress("HcalRecHits_v12",  ROOT.AddressOf( self.hcalHits1 ));
        self.tin1.SetBranchAddress("SimParticles_v12",  ROOT.AddressOf( self.simParticles1 ));
        self.tin1.SetBranchAddress("EcalSimHits_v12",  ROOT.AddressOf( self.ecalSimHits1 ));
        self.tin1.SetBranchAddress("EcalRecHits_v12",  ROOT.AddressOf( self.ecalHits1 ));
        self.tin1.SetBranchAddress("HcalScoringPlaneHits_v12",  ROOT.AddressOf(self.hcalScoringPlaneHits1));
        self.tin1.SetBranchAddress("EcalScoringPlaneHits_v12",  ROOT.AddressOf(self.ecalScoringPlaneHits1));
        self.tin1.SetBranchAddress("TargetScoringPlaneHits_v12",  ROOT.AddressOf(self.targetScoringPlaneHits1));
        self.histograms = {};
        self.critical_angle = []
        self.targetSP_energy_wab_electron = []
        self.HCALSP_energy_wab_electron = []
        self.ECALSP_energy_wab_electron = []
        self.target_detector_diff_wab_electron = []
        self.ECAL_HCAL_diff_wab_electron = []
        
        self.targetSP_energy_wab_photon = []
        self.HCALSP_energy_wab_photon = []
        self.ECALSP_energy_wab_photon = []
        self.target_detector_diff_wab_photon = []
        self.ECAL_HCAL_diff_wab_photon = []
        
        self.hasevent = 0
        self.hasphoton = 0
        self.haselectron = 0
       
        self.loop();
        self.writeOutHistos();

    def writeOutHistos(self):

        self.fout.cd();
      
        
    #WABs:
    def loop(self):
        nentWAB = self.tin1.GetEntriesFast();
        for i in range(nentWAB):
            genEnergy = 0;
            self.tin1.GetEntry(i);
            energy_sum_sim_ECAL_photon = 0;
            energy_sum_rec_ECAL_photon = 0;
            energy_sum_sim_HCAL_photon = 0;
            energy_sum_rec_HCAL_photon = 0;
            energy_sum_HCAL_SP_photon = 0;
            energy_sum_ECAL_SP_photon = 0;
            energy_sum_target_SP_photon = 0;
            
            energy_sum_sim_ECAL_electron = 0;
            energy_sum_rec_ECAL_electron = 0;
            energy_sum_sim_HCAL_electron = 0;
            energy_sum_rec_HCAL_electron = 0;
            energy_sum_HCAL_SP_electron = 0;
            energy_sum_ECAL_SP_electron = 0;
            energy_sum_target_SP_electron = 0;
            particles_per_event = 0
            children = []
            nuceli = []
            PN = []
            Incoming_elec = []
            Incoming_phot = []
            HCAL_elec = []
            HCAL_phot = []
            for ih,particle in enumerate(self.simParticles1):
                #self.histograms["simpart_energy_elec"].Fill(particle.second.getEnergy())
                if particle.second.getProcessType() == 9:
                    print("PN", i, nentWAB, particle.second.getProcessType() ,particle.second.getPdgID(), particle.first )
                    PN.append(particle.first)
                if particle.second.getProcessType() == 0 and (particle.second.getPdgID()!=11 and particle.second.getPdgID() !=22):
                    #print("WAB", particle.second.getPdgID(), particle.first )
                    nuceli.append(particle.first)
                if (particle.second.getPdgID() == 22 and particle.second.getProcessType() == 0):
                    children = particle.second.getDaughters()
                    
            
                if (particle.second.getPdgID() == 11 and particle.second.getProcessType() == 0):
                    Incoming_elec = particle.second.getMomentum()
                    
                    
                if (particle.second.getPdgID() == 22 and particle.second.getProcessType() == 0):
                     Incoming_phot = particle.second.getMomentum()
                """
                if(len(Incoming_phot)!=0):
                    mag_g = math.sqrt(Incoming_phot[0]*Incoming_phot[0] + Incoming_phot[1]*Incoming_phot[1] +Incoming_phot[2]*Incoming_phot[2])
                    mag_e = math.sqrt(Incoming_elec[0]*Incoming_elec[0] + Incoming_elec[1]*Incoming_elec[1] +Incoming_elec[2]*Incoming_elec[2])
                    dat = (Incoming_elec[0]*Incoming_phot[0] + Incoming_elec[1]*Incoming_phot[1] +Incoming_elec[2]*Incoming_phot[2])
                    crit_angle = math.acos(dat/(mag_g*mag_e))
                    self.critical_angle.append(crit_angle)
                    print(crit_angle)
                """
            self.hasevent+=1   
            for ih,hit in enumerate(self.hcalScoringPlaneHits1):
                j = hit.getID() & 0xFFF;
                
                if j==41 and hit.getTrackID()==1:
                    energy = hit.getEnergy();
                    energy_sum_HCAL_SP_electron += energy;
                    HCAL_elec = hit.getMomentum();
                if j==41  and  hit.getTrackID()==2:
                    energy = hit.getEnergy();
                    energy_sum_HCAL_SP_photon += energy;
                    HCAL_phot = hit.getMomentum();
                if(len(HCAL_phot)!=0 and len(HCAL_elec)!=0):
                    
                    mag_g = math.sqrt(HCAL_phot[0]*HCAL_phot[0] + HCAL_phot[1]*HCAL_phot[1] +HCAL_phot[2]*HCAL_phot[2])
                    mag_e = math.sqrt(HCAL_elec[0]*HCAL_elec[0] + HCAL_elec[1]*HCAL_elec[1] +HCAL_elec[2]*HCAL_elec[2])
                    dat = (HCAL_elec[0]*HCAL_phot[0] + HCAL_elec[1]*HCAL_phot[1] +HCAL_elec[2]*HCAL_phot[2])
                    cross = [(HCAL_elec[1]*HCAL_phot[2] - HCAL_elec[2]*HCAL_phot[1]),  -1*(HCAL_elec[0]*HCAL_phot[2] - HCAL_elec[2]*HCAL_phot[0]), (HCAL_elec[0]*HCAL_phot[1] - HCAL_elec[1]*HCAL_phot[0])]
                    mag_cross = math.sqrt(cross[0]*cross[0] + cross[1]*cross[1] +cross[2]*cross[2])
                    print(mag_cross, mag_e)
                    if(mag_g!=0 and mag_e!=0): 
                        crit_angle_dot = math.acos(dat/(mag_g*mag_e))
                        crit_angle_cross = math.asin(mag_cross/(mag_e*mag_g))
                    self.critical_angle.append(crit_angle_dot)
                    #print(crit_angle_cross)
                """
                if len(children) !=0:
                     for i, j in enumerate(children):
                        if hit.getTrackID()==children[i]:
                            energy = hit.getEnergy();
                            energy_sum_HCAL_SP_photon += energy;
                """
            
            for ih,hit in enumerate(self.ecalScoringPlaneHits1):
                j = hit.getID() & 0xFFF;
                if j==34 and hit.getTrackID()==1:
                    energy_sum_ECAL_SP_electron += hit.getEnergy();
                if j==34 and hit.getTrackID()==2:
                    energy_sum_ECAL_SP_photon += hit.getEnergy();
                """
                if len(children) !=0:
                     for i, j in enumerate(children):
                        if hit.getTrackID()==children[i]:
                            energy = hit.getEnergy();
                            energy_sum_ECAL_SP_photon += energy;

                """
            if(energy_sum_ECAL_SP_electron!=0):
                self.haselectron+=1
            if(energy_sum_ECAL_SP_photon!=0):
                self.hasphoton+=1
            for ih,hit in enumerate(self.targetScoringPlaneHits1):
                """
                if i ==1 and len(nuceli)!=0:
                    print(hit.getTrackID())
                """
                i = hit.getID() & 0xFFF ;
                if i ==1 and hit.getTrackID()==1:
                    energy_sum_target_SP_electron += hit.getEnergy();
                if i ==1 and hit.getTrackID()==2:
                    energy_sum_target_SP_photon += hit.getEnergy();

            self.HCALSP_energy_wab_electron.append(energy_sum_HCAL_SP_electron);
            self.ECALSP_energy_wab_electron.append(energy_sum_ECAL_SP_electron);
            self.targetSP_energy_wab_electron.append(energy_sum_target_SP_electron);
            self.target_detector_diff_wab_electron.append(energy_sum_target_SP_electron - energy_sum_HCAL_SP_electron);
            self.ECAL_HCAL_diff_wab_electron.append(energy_sum_target_SP_electron - energy_sum_ECAL_SP_electron);
            
            self.HCALSP_energy_wab_photon.append(energy_sum_HCAL_SP_photon);
            self.ECALSP_energy_wab_photon.append(energy_sum_ECAL_SP_photon);
            self.targetSP_energy_wab_photon.append(energy_sum_target_SP_photon);
            self.target_detector_diff_wab_photon.append(energy_sum_target_SP_photon - energy_sum_HCAL_SP_photon);
            self.ECAL_HCAL_diff_wab_photon.append(energy_sum_target_SP_photon - energy_sum_ECAL_SP_photon);
            
			
def main(options,args) : 
    sc = WabSigCompare(options.ifile1,options.ofile,options.tag) ;
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.critical_angle,
            bins=50, 
            range=(0,math.pi), 
            label="Critical Angle ")
    ax.set_xlabel('Angle Between Electron and Photon')  
    ax.set_yscale('log')
    plt.savefig("crit.png")
    
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.targetSP_energy_wab_photon,
            bins=50, 
            range=(0,4000), alpha=0.5,histtype='step',
            label="Target ")
    n, bins, patches = ax.hist(sc.HCALSP_energy_wab_photon,
            bins=50, 
            range=(0,4000),  alpha=0.5,histtype='step',
            label="HCAL")
    n, bins, patches = ax.hist(sc.ECALSP_energy_wab_photon,
            bins=50, 
            range=(0,4000),  alpha=0.5,histtype='step',
            label="ECAL")
    ax.set_xlabel('Photon Energy [MeV]')  
    ax.set_yscale('log')
    ax.legend()
    plt.savefig("TargetEnergy_Phot.png")
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.target_detector_diff_wab_photon,
            bins=100, 
            range=(-4000,4000), 
            label="SP Hit Energy [MeV] ")
    ax.set_xlabel('Photon Energy Change Target - HCAL [MeV]')  
    ax.set_yscale('log')
    plt.savefig("targetdiff_Phot.png")
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.ECAL_HCAL_diff_wab_photon ,
            bins=100, 
            range=(0,4000), 
            label="SP Hit Energy [MeV] ")
    ax.set_xlabel('Photon Energy Change Target - ECAL [MeV]') 
    ax.set_yscale('log') 
    plt.savefig("diffdet_Phot.png")

    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.targetSP_energy_wab_electron,
            bins=50, 
            range=(0,4000), alpha=0.5,histtype='step',
            label="Target ")
    n, bins, patches = ax.hist(sc.HCALSP_energy_wab_electron,
            bins=50, 
            range=(0,4000),  alpha=0.5,histtype='step',
            label="HCAL")
    n, bins, patches = ax.hist(sc.ECALSP_energy_wab_electron,
            bins=50, 
            range=(0,4000),  alpha=0.5,histtype='step',
            label="ECAL")
    ax.set_xlabel('Electron Energy  [MeV]')  
    ax.set_yscale('log')
    ax.legend()
    plt.savefig("TargetEnergy_Elec.png")
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.target_detector_diff_wab_electron,
            bins=100, 
            range=(0,4000),  
            label="SP Hit Energy [MeV] ")
    ax.set_xlabel('Electron Energy  Change Target - HCAL [MeV]')  
    ax.set_yscale('log')
    plt.savefig("difftarget_Elec.png")
    fig, ax = plt.subplots(1,1)
    n, bins, patches = ax.hist(sc.ECAL_HCAL_diff_wab_electron ,
            bins=100, 
            range=(0,4000), 
            label="SP Hit Energy [MeV] ")
    ax.set_xlabel('Electron Energy Change Target - ECAL[MeV]') 
    ax.set_yscale('log') 
    plt.savefig("diffdet_Elec.png")
    print(sc.hasevent, sc.hasphoton, sc.haselectron)
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
			
			
			
			
			
			
			
