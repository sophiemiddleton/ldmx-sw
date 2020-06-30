/**
 * @file LHEHistoMaker.h
 * @brief Analyzer used for LHE analysis. 
 * @author Sophie Middleton, Caltech
 */

#ifndef _LHEAnalysis_LHEHistMaker_H_
#define _LHEAnalysis_LHEHistMaker_H_

//----------//
//   LDMX   //
//----------//
#include "Framework/EventProcessor.h"

class TH1; 

namespace ldmx { 

    // Forward declarations within the ldmx workspace
    class Event;
    class HistogramPool; 
    class Process;
    class SimParticle;

    class LHEHistoMaker : public Analyzer { 
    
        public: 

            /** Constructor */
            LHEHistoMaker(const std::string &name, Process &process);

            /** Destructor */
            ~LHEHistoMaker();

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

    };    
    
} // ldmx

#endif 
