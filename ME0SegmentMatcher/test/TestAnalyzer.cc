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

class TestAnalyzer : public edm::EDAnalyzer {
public:
  explicit TestAnalyzer(const edm::ParameterSet&);
  ~TestAnalyzer();


  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  void beginJob();

  //protected:
  
  //private:
//Removing this
};

TestAnalyzer::TestAnalyzer(const edm::ParameterSet& iConfig) {}



void TestAnalyzer::beginJob(){}


TestAnalyzer::~TestAnalyzer(){}

void
TestAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

{

  using namespace edm;

  //run_ = (int)iEvent.id().run();
  //event_ = (int)iEvent.id().event();


    //David's functionality
    

  using namespace reco;

  // Handle <ME0MuonCollection > OurMuons;
  // iEvent.getByLabel <ME0MuonCollection> ("me0SegmentMatcher", OurMuons);

  Handle <std::vector<ME0Muon> > OurMuons;
  iEvent.getByLabel <std::vector<ME0Muon> > ("me0SegmentMatcher", OurMuons);
  
  //unsigned int recosize=OurMuons->size();

  //std::cout<<recosize<<std::endl;

  for (std::vector<ME0Muon>::const_iterator thisMuon = OurMuons->begin();
       thisMuon != OurMuons->end(); ++thisMuon){
    std::cout<<thisMuon->pt()<<std::endl;
    //std::cout<<"On a muon"<<std::endl;
  }
  
}

void TestAnalyzer::endJob() {}

DEFINE_FWK_MODULE(TestAnalyzer);
