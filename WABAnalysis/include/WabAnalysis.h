/**
 * @file WabAnalysis.h
 * @brief
 * @author
 */

#ifndef ANALYSIS_WABANALYSIS_H
#define ANALYSIS_WABANALYSIS_H

//LDMX Framework
#include "Event/EventDef.h" // Includes all possible event bus objects
#include "Framework/EventProcessor.h" // Needed to declare processor
#include "Framework/Parameters.h" // Needed to import parameters from configuration file

namespace ldmx {
    
    /**
     * @class WabAnalysis
     * @brief 
     */
    class WabAnalysis : public Analyzer {
        public:

            WabAnalysis(const std::string& name, Process& process) : Analyzer(name, process) {}

            virtual void configure(Parameters& ps);

            virtual void analyze(const Event& event);

            virtual void onFileOpen(EventFile&);

            virtual void onFileClose(EventFile&);

            virtual void onProcessStart(); 

            virtual void onProcessEnd();

        private:

    };
}

#endif /* ANALYSIS_WABANALYSIS_H */
