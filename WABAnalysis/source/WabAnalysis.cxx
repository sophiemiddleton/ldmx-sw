/**
 * @file WabAnalysis.cxx
 * @brief 
 * @author 
 */

#include "WabAnalysis.h"

namespace ldmx {

    void WabAnalysis::configure(Parameters& ps) {

        return;
    }

    void WabAnalysis::analyze(const Event& event) {

        return;
    }
    
    void WabAnalysis::onFileOpen(EventFile&) {

        return;
    }

    void WabAnalysis::onFileClose(EventFile&) {

        return;
    }

    void WabAnalysis::onProcessStart() {

        return;
    }

    void WabAnalysis::onProcessEnd() {

        return;
    }

}

DECLARE_ANALYZER_NS(ldmx, WabAnalysis)
