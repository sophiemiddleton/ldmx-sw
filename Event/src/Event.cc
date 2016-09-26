#include "Event/Event.h"

ClassImp(Event)

int Event::DEFAULT_COLLECTION_SIZE = 100;

std::string Event::SIM_PARTICLES = std::string("SimParticles");

Event::Event() :
        TObject(),
        simParticles(new TClonesArray("SimParticle", Event::DEFAULT_COLLECTION_SIZE)),
        taggerSimHits(new TClonesArray("SimTrackerHit", Event::DEFAULT_COLLECTION_SIZE)),
        recoilSimHits(new TClonesArray("SimTrackerHit", Event::DEFAULT_COLLECTION_SIZE)),
        ecalSimHits(new TClonesArray("SimCalorimeterHit", Event::DEFAULT_COLLECTION_SIZE)),
        hcalSimHits(new TClonesArray("SimCalorimeterHit", Event::DEFAULT_COLLECTION_SIZE)),
        nSimParticles(0),
        nTaggerSimHits(0),
        nRecoilSimHits(0),
        nEcalSimHits(0),
        nHcalSimHits(0) {
}

Event::~Event() {

    Clear();

    delete simParticles;
    delete taggerSimHits;
    delete recoilSimHits;
    delete ecalSimHits;
    delete hcalSimHits;
}

void Event::Clear(Option_t*) {

    TObject::Clear();

    simParticles->Clear("C");
    taggerSimHits->Clear("C");
    recoilSimHits->Clear("C");
    ecalSimHits->Clear("C");
    hcalSimHits->Clear("C");

    nSimParticles = 0;
    nTaggerSimHits = 0;
    nRecoilSimHits = 0;
    nEcalSimHits = 0;
    nHcalSimHits = 0;

    header = EventHeader();
}

EventHeader* Event::getHeader() {
    return &header;
}

int Event::getCollectionSize(const std::string& collectionName) {
    if (collectionName.compare("SimParticles") == 0) {
        return nSimParticles;
    } else if (collectionName.compare("TaggerSimHits") == 0) {
        return nTaggerSimHits;
    } else if (collectionName.compare("RecoilSimHits")) {
        return nRecoilSimHits;
    } else if (collectionName.compare("EcalSimHits")) {
        return nEcalSimHits;
    } else if (collectionName.compare("HcalSimHits")) {
        return nHcalSimHits;
    } else {
        return -1;
    }
}

int Event::nextCollectionIndex(const std::string& collectionName) {
    if (collectionName.compare("SimParticles") == 0) {
        return nSimParticles++;
    } else if (collectionName.compare("TaggerSimHits") == 0) {
        return nTaggerSimHits++;
    } else if (collectionName.compare("RecoilSimHits") == 0) {
        return nRecoilSimHits++;
    } else if (collectionName.compare("EcalSimHits") == 0) {
        return nEcalSimHits++;
    } else if (collectionName.compare("HcalSimHits") == 0) {
        return nHcalSimHits;
    } else {
        return -1;
    }
}

TClonesArray* Event::getCollection(const std::string& collectionName) {
    if (collectionName.compare("SimParticles") == 0) {
        return simParticles;
    } else if (collectionName.compare("TaggerSimHits") == 0) {
        return taggerSimHits;
    } else if (collectionName.compare("RecoilSimHits") == 0) {
        return recoilSimHits;
    } else if (collectionName.compare("EcalSimHits") == 0) {
        return ecalSimHits;
    } else if (collectionName.compare("HcalSimHits") == 0) {
        return hcalSimHits;
    } else {
        return 0;
    }
}

TObject* Event::addObject(const std::string& collectionName) {
    TClonesArray* coll = getCollection(collectionName);
    if (coll != 0) {
        return coll->ConstructedAt(nextCollectionIndex(collectionName));
    } else {
        return 0;
    }
}
