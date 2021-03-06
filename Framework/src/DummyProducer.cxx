/**
 * @file DummyProducer.cxx
 * @brief Class that defines a dummy Producer implementation that adds some random data to the event
 * @author Jeremy Mans, University of Minnesota
 */

// ROOT
#include "TRandom.h"

// LDMX
#include "Event/EventDef.h"
#include "Framework/EventProcessor.h"
#include "Framework/ParameterSet.h"
#include "Framework/Event.h"

// STL
#include <iostream>

namespace ldmx {

    /**
     * @class DummyProducer
     * @brief Defines a dummy Producer implementation that adds some random data to the event
     */
    class DummyProducer : public Producer {

        public:

            DummyProducer(const std::string& name, Process& process) :
                Producer(name, process) {
            }

            virtual void configure(const ParameterSet& ps) {
                nParticles_ = ps.getInteger("nParticles_");
                aveEnergy_ = ps.getDouble("aveEnergy_");
                direction_ = ps.getVDouble("direction_");
            }

            virtual void produce(Event& event) {
                std::cout << "DummyProducer: Analyzing an event!" << std::endl;

                int iEvent = event.getEventHeader().getEventNumber();
                int np = nParticles_*iEvent;
                std::vector<CalorimeterHit> caloHits;
                for (int i = 0; i < np; i++) {
                    caloHits.emplace_back();
                    caloHits.back().setAmplitude( i );
                    caloHits.back().Print();
                }

                event.add("caloHits", caloHits );

            }

            virtual void onFileOpen() {
                std::cout << "DummyProducer: Opening a file!" << std::endl;
            }

            virtual void onFileClose() {
                std::cout << "DummyProducer: Closing a file!" << std::endl;
            }

            virtual void onProcessStart() {
                std::cout << "DummyProducer: Starting processing!" << std::endl;
            }

            virtual void onProcessEnd() {
                std::cout << "DummyProducer: Finishing processing!" << std::endl;
            }

        private:
            int nParticles_{0};
            double aveEnergy_{0};
            std::vector<double> direction_;
            TRandom random_;
    };

}

DECLARE_ANALYZER_NS(ldmx, DummyProducer);
