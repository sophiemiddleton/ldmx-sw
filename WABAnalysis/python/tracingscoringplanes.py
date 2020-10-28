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
        
         #PlaneHits:
        self.histograms["hcalSP_energy_wab_electron"]= ROOT.TH1F("Scoring Plane Energy Elec",";  Scoring Plane Energy [MeV]; Weighted N", 400, 0,4000);
        self.histograms["ecalSP_energy_wab_electron"]= ROOT.TH1F("Scoring Plane Energy Elec","; Scoring Plane Energy [MeV]; Weighted N", 400, 0,4000);
        self.histograms["targetSP_energy_wab_electron"]= ROOT.TH1F("Scoring Plane Energy Elec","; Scoring Plane Energy [MeV]; Weighted N", 400, 0,4000);
        self.histograms["targetSP_Detectors_diff_electron"]= ROOT.TH1F("Energy Loss Elec",";  Scoring Plane Energy [MeV]; Weighted N", 400, 0,4000);
        self.histograms["targetSP_ECALHCAL_diff_electron"]= ROOT.TH1F("Energy Loss Elec",";  Scoring Plane Energy [MeV]; Weighted N", 400, 0,4000);
   
        self.histograms["hcalSP_energy_wab_photon"]= ROOT.TH1F("Scoring Plane Energy Phot",";  Scoring Plane Energy [MeV]; Weighted N", 400, 0,10000);
        self.histograms["ecalSP_energy_wab_photon"]= ROOT.TH1F("Scoring Plane Energy Phot",";  Scoring Plane Energy [MeV]; Weighted N", 400, 0,10000);
        self.histograms["targetSP_energy_wab_photon"]= ROOT.TH1F("Scoring Plane Energy Phot",";  Scoring Plane Energy [MeV]; Weighted N", 400, 0,10000);
        self.histograms["targetSP_Detectors_diff_photon"]= ROOT.TH1F("Energy Loss Phot",";   Scoring Plane Energy [MeV]; Weighted N", 400, 0,10000);
        self.histograms["targetSP_ECALHCAL_diff_photon"]= ROOT.TH1F("Energy Loss Phot",";   Scoring Plane Energy [MeV]; Weighted N", 400, 0,4000);
       
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
            Incoming_elec = []
            for ih,particle in enumerate(self.simParticles1):
                #self.histograms["simpart_energy_elec"].Fill(particle.second.getEnergy())
                if particle.second.getProcessType() == 0 and (particle.second.getPdgID()!=11 and particle.second.getPdgID() !=22):
                    print("WAB", particle.second.getPdgID(), particle.first )
                    nuceli.append(particle.first)
                if (particle.second.getPdgID() == 22 and particle.second.getProcessType() == 0):
                    children = particle.second.getDaughters()
                    
            """
                if (particle.second.getPdgID() == 11 and particle.second.getProcessType() == 0):
                    Incoming_elec = particle.second.getMomentum()
                    print(Incoming_elec)
            """

            for ih,hit in enumerate(self.hcalScoringPlaneHits1):
                
                if hit.getTrackID()==1:
                    energy = hit.getEnergy();
                    energy_sum_HCAL_SP_electron += energy;
                if hit.getTrackID()==2:
                    energy = hit.getEnergy();
                    energy_sum_HCAL_SP_photon += energy;
               
                if len(children) !=0:
                     for i, j in enumerate(children):
                        if hit.getTrackID()==children[i]:
                            energy = hit.getEnergy();
                            energy_sum_HCAL_SP_photon += energy;
                            
            #print("HCAL SP elec", energy_sum_HCAL_SP_electron)
            #print("HCAL SP phot", energy_sum_HCAL_SP_photon)
           
            
            for ih,hit in enumerate(self.ecalScoringPlaneHits1):
                if hit.getTrackID()==1:
                    energy_sum_ECAL_SP_electron += hit.getEnergy();
                if hit.getTrackID()==2:
                    energy_sum_ECAL_SP_photon += hit.getEnergy();
               
                if len(children) !=0:
                     for i, j in enumerate(children):
                        if hit.getTrackID()==children[i]:
                            energy = hit.getEnergy();
                            energy_sum_ECAL_SP_photon += energy;
            #print("ECAL elec", energy_sum_ECAL_SP_electron);
            #print("ECAL phot", energy_sum_ECAL_SP_photon);
            
            for ih,hit in enumerate(self.targetScoringPlaneHits1):
                if i ==1 and len(nuceli)!=0:
                    print(hit.getTrackID())
                i = hit.getID() & 0xFFF ;
                if i ==1 and hit.getTrackID()==1:
                    energy_sum_target_SP_electron += hit.getEnergy();
                if i ==1 and hit.getTrackID()==2:
                    energy_sum_target_SP_photon += hit.getEnergy();
            #print("target Scoring Plane Energy", energy_sum_target_SP_electron);
            #print("target Scoring Plane Energy", energy_sum_target_SP_photon);
            
            
            self.histograms["hcalSP_energy_wab_electron"].Fill(energy_sum_HCAL_SP_electron);
            self.histograms["ecalSP_energy_wab_electron"].Fill(energy_sum_ECAL_SP_electron);
            self.histograms["targetSP_energy_wab_electron"].Fill(energy_sum_target_SP_electron);
            self.histograms["targetSP_Detectors_diff_electron"].Fill(energy_sum_target_SP_electron - energy_sum_HCAL_SP_electron - energy_sum_ECAL_SP_electron);
            self.histograms["targetSP_ECALHCAL_diff_electron"].Fill(energy_sum_HCAL_SP_electron - energy_sum_ECAL_SP_electron);
            
            self.histograms["hcalSP_energy_wab_photon"].Fill(energy_sum_HCAL_SP_photon);
            self.histograms["ecalSP_energy_wab_photon"].Fill(energy_sum_ECAL_SP_photon);
            self.histograms["targetSP_energy_wab_photon"].Fill(energy_sum_target_SP_photon);
            self.histograms["targetSP_Detectors_diff_photon"].Fill(energy_sum_target_SP_photon - energy_sum_HCAL_SP_photon - energy_sum_ECAL_SP_photon);
            self.histograms["targetSP_ECALHCAL_diff_photon"].Fill(energy_sum_HCAL_SP_photon - energy_sum_ECAL_SP_photon);
            #print("Photon difference", energy_sum_HCAL_SP_photon , energy_sum_ECAL_SP_photon, energy_sum_HCAL_SP_photon - energy_sum_ECAL_SP_photon)
def makeCanvas(hists, tags, title, weights, logy=False, norm=False):

    colors = [1,2,6,7]

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
        if i == 0: hists[i].Draw("P");
        else: 
            hists[i].SetLineColor( colors[i] );
            hists[i].Draw('Psames');
            
    leg.Draw();
    if logy:
        ROOT.gPad.SetLogy();
    c.SaveAs("compare_second_"+str(title)+".root")

			
def main(options,args) : 
    sc = WabSigCompare(options.ifile1,options.ofile,options.tag) ;
    weights = [1,1,1]
    histoE1 = [];
    tags_elec = ["HCAL Plane Hits Elec", "ECAL PlaneHits Elec", "Target PlaneHits Elec"]
   
    histoE1.append(sc.histograms["hcalSP_energy_wab_electron"]);
    histoE1.append(sc.histograms["ecalSP_energy_wab_electron"]);
    histoE1.append(sc.histograms["targetSP_energy_wab_electron"]);
    makeCanvas(histoE1, tags_elec, "PlaneHitsElec", weights, True);
   
    histoE2 = [];
    tags_phot = ["HCAL Plane Hits Phot", "ECAL PlaneHits Phot", "Target PlaneHits Phot"]
    histoE2.append(sc.histograms["hcalSP_energy_wab_photon"]);
    histoE2.append(sc.histograms["ecalSP_energy_wab_photon"]);
    histoE2.append(sc.histograms["targetSP_energy_wab_photon"]);
    makeCanvas(histoE2, tags_phot, "PlaneHitsPhot", weights, True);
    
    histoE3 = [];
    tags_diff = ["Electron:E_{target} - E_{ECAL} - E_{HCAL}", "Photon:E_{target} - E_{ECAL}- E_{HCAL}"]
    weights = [1,1]
    histoE3.append(sc.histograms["targetSP_Detectors_diff_electron"]);
    histoE3.append(sc.histograms["targetSP_Detectors_diff_photon"]);
    makeCanvas(histoE3, tags_diff, "EnergyChanges", weights, True);
    
    
    histoE4 = [];
    tags_diff = ["Electron:E_{ECAL} - E_{HCAL}", "Photon:E_{ECAL}- E_{HCAL}"]
    weights = [1,1]
    histoE4.append(sc.histograms["targetSP_ECALHCAL_diff_electron"]);
    histoE4.append(sc.histograms["targetSP_ECALHCAL_diff_photon"]);
    makeCanvas(histoE4, tags_diff, "EnergyChanges", weights, True);

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
			
			
			
			
			
			
			
