/**
 * @file RunManager.h
 * @brief Class providing a Geant4 run manager implementation.
 * @author Jeremy McCormick, SLAC National Accelerator Laboratory
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#ifndef _SIMAPPLICATION_RUNMANAGER_H_
#define _SIMAPPLICATION_RUNMANAGER_H_

//------------//
//   Geant4   //
//------------//
#include "G4RunManager.hh"

//-------------//
//   ldmx-sw   //
//-------------//
#include "Biasing/BiasingMessenger.h"

class G4PhysListFactory; 

namespace ldmx {

    // Forward declare to avoid circular dependency in headers
    class DetectorConstruction; 
    class PluginManager; 
    class PluginMessenger; 

    /**
     * @class RunManager
     * @brief Extension of Geant4 run manager
     */
    class RunManager : public G4RunManager {

        public:

            /**
             * Class constructor.
             */
            RunManager();

            /**
             * Class destructor.
             */
            virtual ~RunManager();

            /**
             * Initialize physics.
             */
            void setupPhysics();

            /**
             * Perform application initialization.
             */
            void Initialize();

            /**
             * Get the user detector construction cast to a specific type.
             * @return The user detector construction.
             */
            DetectorConstruction* getDetectorConstruction(); 

            /** Enable a parallel world. */
            void enableParallelWorld(bool isPWEnabled) { isPWEnabled_ = isPWEnabled; }

            /** Set the path to the GDML description of the parallel world. */
            void setParallelWorldPath(std::string parallelWorldPath) { 
                parallelWorldPath_ = parallelWorldPath; 
            }

            /**
             * Tell RunManager to use the seed from the root file.
             */
            void setUseRootSeed(bool useIt = true) { useRootSeed_ = useIt; }

            /**
             * Should we use the seed from the root file?
             */
            bool useRootSeed() { return useRootSeed_; }

            /**
             * Get access to the plugin manager
             */
            PluginManager* getPluginManager() { return pluginManager_; }

        private:

            /** Plugin messenger. */
            PluginMessenger* pluginMessenger_;

            /** Biasing messenger. */
            BiasingMessenger* biasingMessenger_ {new BiasingMessenger()};

            /**
             * Manager of sim plugins.
             */
            PluginManager* pluginManager_{nullptr};

            /**
             * Factory class for instantiating the physics list.
             */
            G4PhysListFactory* physicsListFactory_{nullptr};

            /** 
             * Flag indicating whether a parallel world should be 
             * registered 
             */
            bool isPWEnabled_{false};

            /** Path to GDML description of parallel world. */
            std::string parallelWorldPath_{""};

            /**
             * Should we use random seed from root file?
             */
            bool useRootSeed_{false};

    }; // RunManager
} // ldmx

#endif // _SIMAPPLICATION_RUNMANAGER_H_
