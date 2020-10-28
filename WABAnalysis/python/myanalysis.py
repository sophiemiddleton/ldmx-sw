#!/usr/bin/python
import argparse
import importlib
import ROOT as r
import os
import math
import sys
import glob
from array import array
from optparse import OptionParser
sys.path.insert(0, '../')

######################################################################

class analyzer:

	def __init__(self, fn, ofn):
                
            # Cut definition
	    self.minVetoHit = 5;
        self.maxTime    = 50;
            
            
	    # histograms
	    self.histograms = {};

            self.histograms["PE_back"]         = r.TH1F("PE_back",     "PE all layers - back;PE;Entries",                200,     0, 1000);
            self.histograms["PE_side"]         = r.TH1F("PE_side",     "PE all layers - side;PE;Entries",                200,     0, 1000);
            self.histograms["Nhit_back"]       = r.TH1F("Nhit_back",   "Hit multiplicity - back;N_{hit};Entries",        100,     0, 100);
            self.histograms["Nhit_side"]       = r.TH1F("NHit_side",   "Hit multiplicity - side;N_{hit};Entries",        100,     0, 100);
            self.histograms["XY_back"]         = r.TH2F("XY_back",     "XY position hits - back;  X (mm) ;  Y (mm)",     640,-1600,1600, 640, -1600,1600);
            self.histograms["XY_side"]         = r.TH2F("XY_side",     "XY position hits - side;  X (mm) ;  Y (mm)",     640,-1600,1600, 640, -1600,1600);
            self.histograms["ZL_back"]         = r.TH2F("ZL_back",     "ZL position hits - back;  Z (mm) ;  L (mm)",    1000,    0,5000, 320,    0,1600);
            self.histograms["ZX_sideTB"]       = r.TH2F("ZX_sideTB",   "ZX position hits - side TB;  Z (mm) ;  X (mm)",  200,    0,1000, 640, -1600,1600);
            self.histograms["ZY_sideTB"]       = r.TH2F("ZY_sideTB",   "ZY position hits - side TB;  Z (mm) ;  Y (mm)",  200,    0,1000, 640, -1600,1600);
            self.histograms["ZX_sideLR"]       = r.TH2F("ZX_sideLR",   "ZX position hits - side LR;  Z (mm) ;  X (mm)",  200,    0,1000, 640, -1600,1600);
            self.histograms["ZY_sideLR"]       = r.TH2F("ZY_sideLR",   "ZY position hits - side LR;  Z (mm) ;  Y (mm)",  200,    0,1000, 640, -1600,1600);
           
            self.histograms["sumPE"]           = r.TH1F("sumPE",          "Sum all PEs;sum PE;Entries",                  200,    0, 10000);
            self.histograms["sumPE_back"]      = r.TH1F("sumPE_back",     "Sum PE all layers - back;PE;Entries",         200,    0, 10000);
            self.histograms["sumPE_side"]      = r.TH1F("sumPE_side",     "Sum PE all layers - side;PE;Entries",         200,    0, 2000);
            self.histograms["maxPE"]           = r.TH1F("maxPE",          "Max PE all layers;PE;Entries",                200,    0, 1000);
            self.histograms["maxPE_back"]      = r.TH1F("maxPE_back",     "Max PE all layers - back;PE;Entries",         200,    0, 1000);
            self.histograms["maxPE_side"]      = r.TH1F("maxPE_side",     "Max PE all layers - side;PE;Entries",         200,    0, 200);
            self.histograms["maxLayer_back"]   = r.TH1F("maxLayer_back",  "First veto layers - back;layer;Entries",      100,    0, 100);
            self.histograms["maxLayer_side"]   = r.TH1F("maxLayer_side",  "First veto layers - side;layer;Entries",      100,    0, 100);

            self.histograms["NoisePE_back"]    = r.TH1F("maxNoisePE_back",  "Noise PE all layers - back;PE;Entries",                 10,    0, 10);
            self.histograms["NoisePE_sideTB"]  = r.TH1F("maxNoisePE_sideTB","Noise PE all layers - side TB;PE;Entries",              10,    0, 10);
            self.histograms["NoisePE_sideLR"]  = r.TH1F("maxNoisePE_sideLR","Noise PE all layers - side LR;PE;Entries",              10,    0, 10);
            self.histograms["NoiseLayer_back"] = r.TH2F("NoiseLayer_back",  "Layer vs Strip noise - back;  Layer ;  Strip",         100, 0, 100, 100, 0, 100);
            self.histograms["NoiseLayer_sideTB"] = r.TH2F("NoiseLayer_sideTB","Layer vs Strip noise - side TB;  Layer ;  Strip",    100, 0, 100,  50, 0, 50);
            self.histograms["NoiseLayer_sideLR"] = r.TH2F("NoiseLayer_sideLR","Layer vs Strip noise - side LR;  Layer ;  Strip",    100, 0, 100,  50, 0, 50);
          
	    #self.histograms["Veto_back"]      = r.TH1F("Veto_back",     "Max PE for single particle - back;PE;Entries",         200,    0, 1000);

            for key in self.histograms: 
                self.histograms[key].GetXaxis().SetTitleSize(0.045);
                self.histograms[key].GetXaxis().SetLabelSize(0.035);
                self.histograms[key].GetYaxis().SetTitleSize(0.045);
                self.histograms[key].GetYaxis().SetLabelSize(0.035);
                if (type(self.histograms[key]) is r.TH2F) : self.histograms[key].GetYaxis().SetTitleOffset(1.6);
            
	    
	    
	    #Process data	    
            for thisFile in glob.glob(fn) : 
                print "Process "+str(thisFile)
                self.processFile(thisFile)
                print "Done looping\n";

	    
	    #Make plots
            plotLog=["sumPE","sumPE_back","sumPE_side","maxPE","maxPE_back","maxPE_side","maxLayer_back","maxLayer_side","PE_back","PE_side",
	             "Nhit_back","Nhit_side","NoisePE_back","NoisePE_side","Veto_back"]
            
            canvas = r.TCanvas("c1","c1",800,800)
            for key in self.histograms: 
                 canvas.SetLogy(0)
                 if (key in plotLog): canvas.SetLogy(1)
                 if (type(self.histograms[key]) is r.TH2F) : self.histograms[key].Draw("colz")
                 else                                      : self.histograms[key].Draw()
                 canvas.SaveAs("plots/"+key+".pdf")                 


	    #Save histos
            self.fn_out = ofn;
	    self.fout = r.TFile(self.fn_out,"RECREATE");
            #self.fout.cd();
	    for key in sorted(self.histograms.iterkeys()): 
		    self.histograms[key].Write();
	    self.fout.Close();

	
	
	####################################################
	### Process files
	####################################################
	def processFile(self,filename):
	    
	    self.fin = r.TFile(filename);
	    self.tin = self.fin.Get("LDMX_Events")

	    self.evHeader     = r.ldmx.EventHeader()
	    #self.ecalSPHits   = r.std.vector('ldmx::SimTrackerHit')()
	    self.hcalHits     = r.std.vector('ldmx::HcalHit')()

	    self.tin.SetBranchAddress("EventHeader",     r.AddressOf( self.evHeader ));
            self.tin.SetBranchAddress("HcalRecHits_sim", r.AddressOf( self.hcalHits ));
            #no scoring plane hits yet
	    #self.tin.SetBranchAddress("EcalScoringPlaneHits_sim", r.AddressOf( self.ecalSPHits ));

            #start looping
            nent = self.tin.GetEntriesFast();
            #nent=100
	    for i in range(nent):
		   if (nent/100 > 0 and i % (1 * nent/100) == 0):
			    sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done")
			    sys.stdout.flush()
		   self.tin.GetEntry(i);
                   self.analyze(i);
	    
	    self.fin.Close()
	    print ' '


	####################################################
	### Analyze
	####################################################
	def analyze(self,ievt):
            
            maxPE_back = 0; maxPE_side=0; sumPE_back=0; sumPE_side=0; maxLayer_back=9999; maxLayer_side = 9999;nHits_back=0;nHits_side=0;
            nVetoBack=0; nVetoSide=0;
            
	    #print str(ievt)+" "+str(self.evHeader.getEventNumber() ) +" "+ str(len(self.hcalHits))

            for hit in self.hcalHits:
                section = hit.getSection()
                layer   = hit.getLayer()
                strip   = hit.getStrip()
                PE      = hit.getPE();
                x       = hit.getX();
                y       = hit.getY();
                z       = hit.getZ();
                l       = max(x,y); 
                time    = hit.getTime();
                
		
		#do noise first, we cut on this after...
		if (hit.getNoise()==1):		
		    if (section==0) : 
		        self.histograms["NoisePE_back"].Fill(PE); 
		        self.histograms["NoiseLayer_back"].Fill(layer,strip); 
		    
		    if (section==1 or section==2): 
		        self.histograms["NoisePE_sideTB"].Fill(PE);	    
		        self.histograms["NoiseLayer_sideTB"].Fill(layer,strip); 

		    if (section==3 or section==4): 
		        self.histograms["NoisePE_sideLR"].Fill(PE);	    
		        self.histograms["NoiseLayer_sideLR"].Fill(layer,strip); 
		
		
                if (time <= -100)          : continue;
                if (time >= self.maxTime)  : continue;
                if (hit.getMinPE()<1)      : continue;
		                
		#back Hcal
		if (section==0):                    
                    
                    nHits_back += 1
		    sumPE_back += PE;
                    self.histograms["PE_back"].Fill(PE);
                    self.histograms["XY_back"].Fill(x,y);
                    self.histograms["ZL_back"].Fill(z,l);
                                                            
                    if (PE > maxPE_back):    
                        maxPE_back = PE;
                        maxLayer_back = layer;  
					    
			
                #Top bottom side HCal                           
                if (section==1 or section==2):
                    nHits_side += 1
                    sumPE_side += PE;
                    self.histograms["PE_side"].Fill(PE);
                    self.histograms["XY_side"].Fill(x,y);
                    self.histograms["ZX_sideTB"].Fill(z,x);
                    self.histograms["ZY_sideTB"].Fill(z,y);
                    
                    if (PE > maxPE_side):    
                        maxPE_side = PE;
                        maxLayer_side = layer;  
                                        
                
		#Left right side Hcal
		if (section==3 or section==4):
                    nHits_side += 1
                    sumPE_side += PE;
                    self.histograms["PE_side"].Fill(PE);
                    self.histograms["XY_side"].Fill(x,y);
                    self.histograms["ZX_sideLR"].Fill(z,x);
                    self.histograms["ZY_sideLR"].Fill(z,y);                   
                         
                    if (PE > maxPE_side):    
                        maxPE_side = PE;
                        maxLayer_side = layer;  
		                       
                    
            maxPE_all = max(maxPE_side,maxPE_back)
            sumPE_all = sumPE_back+sumPE_side
            
	    
	    self.histograms["Nhit_back"].Fill(nHits_back);
	    self.histograms["Nhit_side"].Fill(nHits_side);

            if (sumPE_back>0): self.histograms["sumPE_back"].Fill(sumPE_back);
            if (sumPE_side>0): self.histograms["sumPE_side"].Fill(sumPE_side);
            if (sumPE_all>0) : self.histograms["sumPE"].Fill(sumPE_all);
                
            self.histograms["maxLayer_back"].Fill(maxLayer_back);
            self.histograms["maxLayer_side"].Fill(maxLayer_side);            
            if (maxPE_all>0)     : self.histograms["maxPE"].Fill(maxPE_all);
            if (maxPE_back>0)    : self.histograms["maxPE_back"].Fill(maxPE_back);
            if (maxPE_side>0)    : self.histograms["maxPE_side"].Fill(maxPE_side);
            
	    #if (self.isSingleBack ==1):   self.histograms["Veto_back"].Fill(maxPE_back);  
                        
            #print str(ievt)+" " +str(nNoise)+" "+str(len(self.hcalHits))
            

                 




	####################################################
	### Utils
	####################################################
           
        def isSingleBack(self):            
            
	    nLow=0;nHighBack=0;nHighSide=0;
	    for ih,h in enumerate(self.ecalSPHits):
		sim = h.getSimParticle()
                
		## don't care if stuff goes backward 
		if (h.getLayerID()==1) : continue;
		
		## at least 10 MeV of kinetic energy 
		ke = self.kineMom(h.getMomentum(),sim.getMass());                                
                if (ke<10) : continue;
                
                ##
		if (ke<1000) :
		     nLow +=1;
		     continue
		    
		if (h.getLayerID()==2) : nHighBack += 1
		if (h.getLayerID()>2)  : nHighSide += 1
            
            if (nHighBack==1 and nHighSide==0 and nLow==0) : return 1
            return 0
        

        def dumpSim(self):
            print("             PDG ID    KE       THETA   PHI           START PT                 END PT")
            for ih,h in enumerate(self.simParticles):
                if (self.kine(h)<1.0) : continue;
                print('SIM  {0:2d} {1:11d}  {2:7.2f}   {3:6.2f}  {4:6.2f}    ({5:5.2f},{6:5.2f},{7:5.2f})   ({8:5.2f},{9:5.2f},{10:5.2f})'.format(
                ih,h.getPdgID(),self.kine(h),self.polar(h.getMomentum())*57.2958, self.phi(h.getMomentum())*180/3.1415, h.getVertex()[0], h.getVertex()[1], h.getVertex()[2],
                h.getEndPoint()[0], h.getEndPoint()[1], h.getEndPoint()[2]   ))

        def dumpHcalHits(self):
            for hit in self.hcalHits: 
              if (hit.getPE()<2): continue;
	      hx = hit.getX();
	      hy = hit.getY();
	      hz = hit.getZ();
	      hm = math.sqrt(hx*hx+hy*hy+hz*hz);
	      print "PE=" + str(hit.getPE())+"  phi="+str(math.atan2(hy,hx)*180/3.1415)+"   pos=("+str(hx)+","+str(hy)+","+str(hz)+")";
                
        def polar(self,vec):
            if (self.mag(vec) < 0.001) : return -999;
            return math.acos(vec[2]/math.sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]));

        def phi(self,vec):
            return math.atan2(vec[2],vec[1]);

	def mag(self,vec):
            return math.sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);

	def kineMom(self,vec,mass):
            return math.sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]+mass*mass)-mass;

	def kine(self,sim):
            return sim.getEnergy()-sim.getMass();

######################################################################


def main(options,args) : 
	
	print "Start"
        print options.ifile
	analyzer(options.ifile,options.ofile);





if __name__ == "__main__":
	

	parser = OptionParser()
	parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
	parser.add_option('-i','--ifile', dest='ifile', default = 'file.root',help='directory with data', metavar='idir')
	parser.add_option('-o','--ofile', dest='ofile', default = 'ofile.root',help='directory to write plots', metavar='odir')

	(options, args) = parser.parse_args()

	r.gStyle.SetPadBottomMargin(0.15)
	r.gStyle.SetPadLeftMargin(0.15)
	r.gStyle.SetPadRightMargin(0.12)
	r.gStyle.SetPalette(55)
	r.gStyle.SetOptFit(0)
	r.gStyle.SetOptStat(0)
	r.gStyle.SetFrameBorderMode(0)
	r.gStyle.SetCanvasBorderMode(0)
	r.gStyle.SetPadBorderMode(0)
	r.gStyle.SetPadColor(0)
	r.gStyle.SetCanvasColor(0)
	r.gStyle.SetStatColor(0)
	r.gStyle.SetTitleFillColor(0)
	r.gStyle.SetLineWidth(2)
	r.gStyle.SetHistLineWidth(2)
	r.gROOT.SetBatch()

	# Get the Event library 
	r.gSystem.Load("/nfs/slac/g/ldmx/users/echenard/ldmx-sw/install/lib/libEvent.so");	


main(options,args);

