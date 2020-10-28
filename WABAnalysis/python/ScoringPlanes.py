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
        self.histograms["hcalSP_energy_wab"]= ROOT.TH1F("HCAL Scoring Plane Energy",";  HCAL Scoring Plane Energy [MeV]; Weighted N", 400, 0,10000);
        self.histograms["ecalSP_energy_wab"]= ROOT.TH1F("ECAL Scoring Plane Energy",";  ECAL Scoring Plane Energy [MeV]; Weighted N", 400, 0,10000);
        self.histograms["targetSP_energy_wab"]= ROOT.TH1F("Target Scoring Plane Energy",";  Target Scoring Plane Energy [MeV]; Weighted N", 400, 0,10000);
        
        #SimParts:
        self.histograms["simpart_energy_wab"]= ROOT.TH1F("Total Event E SimParts_wab",";  Total Event Energy [MeV]; Weighted N", 400, 0,10000);
        self.histograms["simpart_energy_elec"]= ROOT.TH1F("Total Event E SimHParts_elc","; Total Event Energy [MeV]; Weighted N", 1000, 0,4000);
        self.histograms["simpart_energy_phot"]= ROOT.TH1F("Total Event E SimHParts_phot",";  Total Event Energy [MeV]; Weighted N", 1000, 0,4000);
        self.histograms["simpart_energy_n"]= ROOT.TH1F("Total Event E SimHParts_n",";  Total Event Energy [MeV]; Weighted N", 1000, 0,4000);

        
        #SimHits WAB:
        self.histograms["simhit_edep_wabHCAL"] = ROOT.TH1F("HCAL Energy Dep per SimHit_wab",";HCAL Sim Hit Energy Dep [MeV]; Weighted N", 400, 0,2);
        self.histograms["simhit_totE_wabHCAL"]= ROOT.TH1F("HCAL Total Event E SimHits_wab","; HCAL Total Event Energy [MeV]; Weighted N", 400, 0,4000); 
        self.histograms["simhit_totE_wabECAL"]= ROOT.TH1F("ECAL Total Event E SimHits_wab","; HCAL Total Event Energy [MeV]; Weighted N", 400, 0,4000);
        self.histograms["simhit_totE_wabECALHCAL"]= ROOT.TH1F("ECAL and HCAL Total Event E SimHits_wab","; HCAL and ECAL Total Event Energy [MeV]; Weighted N", 400, 0,4000);
        
        #RecoHits WABS:
        self.histograms["rechit_totE_wabHCAL"]= ROOT.TH1F("HCAL Total Event E RecHits_wab","; HCAL Total Event Energy [MeV]; Weighted N", 400, 0,4000);
        self.histograms["rechit_edep_wabHCAL"] = ROOT.TH1F("HCAL Energy Dep per Rec Hits_wab","; HCAL Rec Hit Energy Dep [MeV]; Weighted N", 400, 0,10);
        self.histograms["rechit_totE_wabECAL"]= ROOT.TH1F("ECAL Total Event E RecHits_wab","; ECAL Total Event Energy [MeV]; Weighted N", 400, 0,4000);
        self.histograms["rechit_totE_wabECALHCAL"]= ROOT.TH1F("ECAL and HCAL Total Event E RecHits_wab","; HCAL and ECAL Total Event Energy [MeV]; Weighted N", 400, 0,4000);
        self.histograms["rechit_diffE_wabECALHCAL"]= ROOT.TH1F("ECAL and HCAL Missed E RecHits_wab","; HCAL and ECAL Missing Energy [MeV]; Weighted N", 400, 0,4000);
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
            energy_sum_sim_ECAL = 0;
            energy_sum_rec_ECAL = 0;
            energy_sum_sim_HCAL = 0;
            energy_sum_rec_HCAL = 0;
            energy_sum_HCAL_SP = 0;
            energy_sum_ECAL_SP = 0;
            energy_sum_target_SP = 0;
            particles_per_event = 0
            for ih,particle in enumerate(self.simParticles1):
                if(abs(particle.second.getPdgID())) == 11:
                    self.histograms["simpart_energy_elec"].Fill(particle.second.getEnergy())
                if(particle.second.getPdgID()) == 22:
                    self.histograms["simpart_energy_phot"].Fill(particle.second.getEnergy())
                if(particle.second.getPdgID()) > 10000:
                    self.histograms["simpart_energy_n"].Fill(particle.second.getEnergy())
                genEnergy += particle.second.getEnergy();
                particles_per_event +=1
                print("Particle", particle.second.getPdgID(), "Process",particle.second.getProcessType())
            print("Sim Energy", genEnergy)
            
            for ih,hit in enumerate(self.hcalScoringPlaneHits1):
                energy = hit.getEnergy();
                energy_sum_HCAL_SP += energy;
                print(ih, "HCAL SP", energy_sum_HCAL_SP)
            print("HCAL Scoring Plane Energy", energy_sum_HCAL_SP);
            
            for ih,hit in enumerate(self.ecalScoringPlaneHits1):
                energy_sum_ECAL_SP += hit.getEnergy();
            print("ECAL Scoring Plane Energy", energy_sum_ECAL_SP);
            
            for ih,hit in enumerate(self.targetScoringPlaneHits1):
                
                i = hit.getID() & 0xFFF;
                if i ==1:
                    energy_sum_target_SP += hit.getEnergy();
            print("target Scoring Plane Energy", energy_sum_target_SP);
            
            for ih,hit in enumerate(self.hcalSimHits1):
                energy_sum_sim_HCAL += hit.getEdep();
            print("HCAL Energy SimHits", energy_sum_sim_HCAL)
            
            for ih,hit in enumerate(self.hcalHits1):		
                if (hit.isNoise()==0):
                    energy_sum_rec_HCAL += abs(hit.getEnergy());
                    
            print("HCAL Energy RecHits", energy_sum_rec_HCAL)
            for ih,hit in enumerate(self.ecalSimHits1):
                energy_sum_sim_ECAL += hit.getEdep();
            print("ECAL Energy SimHits", energy_sum_sim_ECAL)
           
            for ih,hit in enumerate(self.ecalHits1):		
                if (hit.isNoise()==0  and hit.getXPos()!=0 and hit.getYPos()!=0 and hit.getZPos()!=0):
                    energy_sum_rec_ECAL += abs(hit.getAmplitude());
            print("ECAL Energy RecHits", energy_sum_rec_ECAL)
             
            print("Detector Sum SimHits",  energy_sum_sim_ECAL+ energy_sum_sim_HCAL)
            print("Detector Sum RecHits",  energy_sum_rec_ECAL+ energy_sum_rec_HCAL)
            print("========================")

            self.histograms["simhit_totE_wabECAL"].Fill(energy_sum_sim_ECAL);
            self.histograms["rechit_totE_wabECAL"].Fill(energy_sum_rec_ECAL);
            self.histograms["simhit_totE_wabECALHCAL"].Fill(energy_sum_sim_HCAL+energy_sum_sim_ECAL);
          
            self.histograms["rechit_totE_wabECALHCAL"].Fill(energy_sum_rec_ECAL+energy_sum_rec_HCAL);
            self.histograms["simhit_totE_wabHCAL"].Fill(energy_sum_sim_HCAL);
            self.histograms["rechit_totE_wabHCAL"].Fill(energy_sum_rec_HCAL);
            self.histograms["simpart_energy_wab"].Fill(genEnergy);
            diffE = genEnergy - energy_sum_sim_HCAL - energy_sum_sim_ECAL;
            self.histograms["rechit_diffE_wabECALHCAL"].Fill(diffE);
            self.histograms["hcalSP_energy_wab"].Fill(energy_sum_HCAL_SP);
            self.histograms["ecalSP_energy_wab"].Fill(energy_sum_ECAL_SP);
            self.histograms["targetSP_energy_wab"].Fill(energy_sum_target_SP);
        
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
    c.SaveAs("compare_plot_SimPart"+str(title)+".root")

			
def main(options,args) : 
    sc = WabSigCompare(options.ifile1,options.ofile,options.tag) ;
    weights = [1,1,1]#,1]#factor_wab, factor_wab, factor_wab, factor_wab]
    histoE1 = [];
    tags = ["HCAL Plane Hits", "ECAL PlaneHits", "Target PlaneHits"]
    #histoE1.append(sc.histograms["simpart_energy_wab"]);
    histoE1.append(sc.histograms["hcalSP_energy_wab"]);
    histoE1.append(sc.histograms["ecalSP_energy_wab"]);
    histoE1.append(sc.histograms["targetSP_energy_wab"]);
    makeCanvas(histoE1, tags, "SimHitTotE", weights,True);
   
    """
    histoE2 = [];
    histoE2.append(sc.histograms["simpart_energy_wab"]);
    histoE2.append(sc.histograms["rechit_totE_wabECALHCAL"]);
    histoE2.append(sc.histograms["rechit_totE_wabHCAL"]);
    histoE2.append(sc.histograms["rechit_totE_wabECAL"]);
    makeCanvas(histoE2, tags, "RecHitTotE", weights, True);
    
    histoElec = [];
    tags = ["electrons", "photons","nuceli"];
    histoElec.append(sc.histograms["simpart_energy_elec"]);
    histoElec.append(sc.histograms["simpart_energy_phot"]);
    histoElec.append(sc.histograms["simpart_energy_n"]);
    makeCanvas(histoElec, tags, "SimPartEdeps", weights, True);
    
    histoDiff = [];
    tags = ["Incoming - ECAL - HCAL"];
    histoDiff.append(sc.histograms["rechit_diffE_wabECALHCAL"]);
    makeCanvas(histoDiff, tags, "EnergyLoss", weights, True);
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
			
			
			
			
			
			
			
