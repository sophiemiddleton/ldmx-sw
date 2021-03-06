/**
 * @file HCalDQM.h
 * @brief Analyzer used for HCal DQM. 
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#ifndef _DQM_HCAL_DQM_H_
#define _DQM_HCAL_DQM_H_

//----------//
//   STL    //
//----------//
#include <algorithm>

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
#include "Framework/EventProcessor.h"
#include "Tools/AnalysisUtils.h"
#include "Event/EventDef.h"

namespace ldmx { 

    class HCalDQM : public Analyzer { 
    
        public: 

            /** Constructor */
            HCalDQM(const std::string &name, Process &process);

            /** Destructor */
            ~HCalDQM();

            /** 
             * Configure the processor using the given user specified parameters.
             * 
             * @param pSet Set of parameters used to configure this processor.
             */
            void configure(const ParameterSet &pSet);
 
            /**
             * Process the event and make histograms ro summaries.
             *
             * @param event The event to analyze.
             */
            void analyze(const Event& event);

            /** Method executed before processing of events begins. */
            void onProcessStart();

        private:

            /** Singleton used to access histograms. */
            HistogramPool* histograms_{nullptr}; 

            /** The maximum PE threshold used for the veto. */
            float maxPEThreshold_{5}; 
            
            /** Name of ECal veto collection. */
            std::string ecalVetoCollectionName_{"EcalVeto"}; 
    };    
    
} // ldmx

#endif // _DQM_HCAL_DQM_H_
