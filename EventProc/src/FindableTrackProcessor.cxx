/**
 * @file FindableTrackProcessor.cxx
 * @brief Processor used to find all particles that pass through the recoil
 *        tracker and leave hits consistent with a findable track.
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#include "EventProc/FindableTrackProcessor.h"

namespace ldmx { 

    FindableTrackProcessor::FindableTrackProcessor(const std::string &name, Process &process) :
            Producer(name, process) { }

    FindableTrackProcessor::~FindableTrackProcessor() { }

    void FindableTrackProcessor::configure(const ParameterSet &pset) { }

    void FindableTrackProcessor::produce(Event &event) {

        // Get the collection of sim particles from the event 
        const std::map<int,SimParticle> simParticles = event.getMap<int,SimParticle>("SimParticles");
        if (simParticles.size() == 0) return; 

        // Get the collection of Recoil sim hits from the event
        const std::vector<SiStripHit> siStripHits = event.getCollection<SiStripHit>("SiStripHits");  

        // Create the hit map
        this->createHitMap(siStripHits); 
        
        // Loop through all sim particles and check which are findable 
        int resultCount = 0;
        std::vector<FindableTrackResult> findableTrackResults;
        for (const auto &keyVal : simParticles ) {
            
            const SimParticle *simParticle = &(keyVal.second);

            // If the sim particle is neutral, skip it.
            if (abs(simParticle->getCharge()) != 1) continue;
    
            // Check if the track is findable 
            if (hitMap_.count(keyVal.first) == 1) {
                
                // Create a result instance
                findableTrackResults.emplace_back();
                
                // Set the sim particle associated with the result
                findableTrackResults.back().setParticleTrackID(keyVal.first);

                // Check if the track is findable
                this->isFindable( findableTrackResults.back(), hitMap_[keyVal.first]); 
                resultCount++;
            }      
        }

        //Add the result to the collection
        event.add( "FindableTracks", findableTrackResults );
    }

    void FindableTrackProcessor::createHitMap(const std::vector<SiStripHit> &siStripHits) { 
       
        // Clear the hit map to remove any previous relations 
        hitMap_.clear();

        // Loop over all recoil tracker sim hits and check which layers, if any,
        // the sim particle deposited energy in.
        std::vector<int> recoilHitCount(10, 0);
        for ( const SiStripHit &siStripHit : siStripHits ) {

            // Get the SimTrackerHit associated with this strip hit
            const SimTrackerHit simTrackerHit = siStripHit.getSimTrackerHits().at(0);  

            // Get the MC particle associated with this hit
            int simParticleTrackID = simTrackerHit.getTrackID();

            // If the particle ins't in the hit map, add it.
            if (hitMap_.count(simParticleTrackID) == 0) {
                hitMap_[simParticleTrackID] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            }
            
            // Increment the hit count for the layer
            hitMap_[simParticleTrackID][simTrackerHit.getLayerID() - 1]++;
        }
    }

    void FindableTrackProcessor::isFindable(FindableTrackResult &result, std::vector<int> hitCount) { 
       
        /*std::cout << "[ FindableTrackProcessor ]: Hit Vec [ ";
        for (auto index : hitCount) { 
            std::cout << index << ", ";
        }
        std::cout << " ]" << std::endl;*/

        // Count how many 3D stereo hits are created by this particle
        std::vector<int> hit3dVec{0, 0, 0, 0};
        int hit3dCount{0};
        for (int layerN = 0; layerN  < 8; layerN += 2) { 
            if (hitCount[layerN]*hitCount[layerN+1] != 0) {
                hit3dCount++;
                hit3dVec[layerN/2]++; 
            }
        }
        /*std::cout << "[ FindableTrackProcessor ]: " 
                  << "Total 3D Hits: " << hit3dCount << std::endl;*/
         

        // A track is considered findable if 
        // 1) The first four stereo layers are hit
        // 2) Three of the first four layers are hit and an axial layer is hit
        // 3) Two of the first four layers are hit and both axial layers are hit
        bool trackFound{false};
        if (hit3dCount == 4) { 
            result.setResult(FindableTrackResult::STRATEGY_4S, true); 
            trackFound = true;
            /*std::cout << "[ FindableTrackProcessor ]: " 
                      << "Can be found using 4S strategy." << std::endl;*/
        } 
        
        if ((hit3dVec[0]*hit3dVec[1]*hit3dVec[2] > 0) && (hitCount[8] > 0 || hitCount[9] > 0)) {
            result.setResult(FindableTrackResult::STRATEGY_3S1A, true); 
            trackFound = true;
            /*std::cout << "[ findabletrackprocessor ]: " 
                      << "can be found using 3s1a strategy." << std::endl;*/
        } 
        
        if (hit3dVec[0]*hit3dVec[1] > 0 && (hitCount[8] > 0 && hitCount[9] > 0)) {
            result.setResult(FindableTrackResult::STRATEGY_2S2A, true); 
            trackFound = true;
            /*std::cout << "[ findabletrackprocessor ]: " 
                      << "can be found using 2s2a strategy." << std::endl;*/
        } 
        
        if (hitCount[8] > 0 && hitCount[9] > 0) { 
            result.setResult(FindableTrackResult::STRATEGY_2A, true); 
            trackFound = true;
            /*std::cout << "[ findabletrackprocessor ]: " 
                      << "can be found using 2a strategy." << std::endl;*/
        }

        if (hit3dVec[0]*hit3dVec[1] > 0) {
            result.setResult(FindableTrackResult::STRATEGY_2S, true); 
            trackFound = true;
            /*std::cout << "[ findabletrackprocessor ]: " 
                      << "can be found using 2s strategy." << std::endl;*/
        } 

        if (hit3dVec[0]*hit3dVec[1]*hit3dVec[2] > 0) {
            result.setResult(FindableTrackResult::STRATEGY_3S, true); 
            trackFound = true;
            /*std::cout << "[ findabletrackprocessor ]: " 
                      << "can be found using 3s strategy." << std::endl;*/
        } 

        if (!trackFound) { 
            result.setResult(FindableTrackResult::STRATEGY_NONE, false);    
        }
    }
}

DECLARE_PRODUCER_NS(ldmx, FindableTrackProcessor) 
