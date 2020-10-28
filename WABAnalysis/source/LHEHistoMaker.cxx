/**
 * @file LHEHistoMaker.cxx
 * @brief Analyzer used for ECal PN DQM.
 * @author Sophie Middleton, Caltech
 */

#include "LHEAnalysis/include/LHEHistoMaker.h"

//----------//
//   ROOT   //
//----------//
#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"

//----------//
//   LDMX   //
//----------//

#include "Framework/Event.h"
#include "Framework/HistogramPool.h"
#include "Tools/AnalysisUtils.h"

namespace ldmx { 

    LHEHistoMaker::LHEHistoMaker(const std::string &name, Process &process) : 
        Analyzer(name, process) { }

    LHEHistoMaker::~LHEHistoMaker() {}

    void LHEHistoMaker::onProcessStart() {
      

        // Get an instance of the histogram pool  
        histograms_ = HistogramPool::getInstance();    

        std::vector<std::string> labels = {"", 
            "Nothing hard", // 0  
            "1 theta", // 1
            "2 pt", // 2
            
            ""
        };

        std::vector<TH1*> hists = { 
            histograms_->get("theta"),
            histograms_->get("pt")
        };


        for (int ilabel{1}; ilabel < labels.size(); ++ilabel) { 
            for (auto& hist : hists) {
                hist->GetXaxis()->SetBinLabel(ilabel, labels[ilabel-1].c_str());
            }
        }

        labels = {"", 
            "#theta n", // 0
            "p_{t} n", // 1
            ""
        };

        hists = {
            histograms_->get("event_type_compact"),
            histograms_->get("event_type_compact_track_veto"),
            
        };

        for (int ilabel{1}; ilabel < labels.size(); ++ilabel) { 
            for (auto& hist : hists) {
                hist->GetXaxis()->SetBinLabel(ilabel, labels[ilabel-1].c_str());
            }
        }


        getHistoDirectory();

       
    }

    void LHEHistoMaker::configure(const ParameterSet& ps) {
       
    }

    void LHEHistoMaker::analyze(const Event & event) { 
 
        }

} // ldmx

DECLARE_ANALYZER_NS(ldmx, LHEHistoMaker)
