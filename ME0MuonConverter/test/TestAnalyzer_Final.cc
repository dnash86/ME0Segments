#include <ME0Reconstruction/ME0SegmentMatcher/src/ME0SegmentMatcher.h>
//#include <RecoLocalMuon/ME0Segment/src/ME0SegmentBuilder.h>

#include <FWCore/PluginManager/interface/ModuleDef.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 

#include <Geometry/Records/interface/MuonGeometryRecord.h>

#include <ME0Reconstruction/ME0Segment/interface/ME0Segment.h>
#include <ME0Reconstruction/ME0Segment/interface/ME0SegmentCollection.h>

#include <ME0Reconstruction/ME0Segment/interface/ME0Muon.h>
#include <ME0Reconstruction/ME0Segment/interface/ME0MuonCollection.h>

// #include "CLHEP/Matrix/SymMatrix.h"
// #include "CLHEP/Matrix/Matrix.h"
// #include "CLHEP/Vector/ThreeVector.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
//#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TLorentzVector.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
//#include "TRandom3.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"


#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"


class TestAnalyzer_Final : public edm::EDAnalyzer {
public:
  explicit TestAnalyzer_Final(const edm::ParameterSet&);
  ~TestAnalyzer_Final();


  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  void beginJob();

  //protected:
  
  //private:
//Removing this
};

TestAnalyzer_Final::TestAnalyzer_Final(const edm::ParameterSet& iConfig) {}



void TestAnalyzer_Final::beginJob(){}


TestAnalyzer_Final::~TestAnalyzer_Final(){}

void
TestAnalyzer_Final::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

{

  using namespace edm;

  //run_ = (int)iEvent.id().run();
  //event_ = (int)iEvent.id().event();


    //David's functionality
    

  using namespace reco;

  // Handle <ME0MuonCollection > OurMuons;
  // iEvent.getByLabel <ME0MuonCollection> ("me0SegmentMatcher", OurMuons);

  Handle <std::vector<RecoChargedCandidate> > OurCandidates;
  iEvent.getByLabel <std::vector<RecoChargedCandidate> > ("me0MuonConverter", OurCandidates);
  
  //unsigned int recosize=OurCandidates->size();

  //std::cout<<recosize<<std::endl;

  for (std::vector<RecoChargedCandidate>::const_iterator thisCandidate = OurCandidates->begin();
       thisCandidate != OurCandidates->end(); ++thisCandidate){

    std::cout<<"On a muon"<<std::endl;
    std::cout<<thisCandidate->eta()<<std::endl;
  }
  
}

void TestAnalyzer_Final::endJob() {}

DEFINE_FWK_MODULE(TestAnalyzer_Final);
