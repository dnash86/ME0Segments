

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
//#include "Geometry/RPCGeometry/interface/RPCRoll.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"


#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
//#include "CLHEP/Matrix/DiagMatrix.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TLorentzVector.h"
#include "TRandom1.h"

#include "TFile.h"
#include "TTree.h"



#include <map>

//#include "MuProp_v2.h"
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <iostream>
//#include "CMSStyle.C"
#include "TStyle.h"
#include "TMath.h"
//
// class decleration
//

class TestAnalyzer : public edm::EDAnalyzer {
public:
  explicit TestAnalyzer(const edm::ParameterSet&);
  ~TestAnalyzer();


  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  void beginJob();

protected:
  
private:
};

TestAnalyzer::TestAnalyzer(const edm::ParameterSet& iConfig) : {}



void TestAnalyzer::beginJob(){}


TestAnalyzer::~TestAnalyzer(){}

void
TestAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

{

  using namespace edm;

  run_ = (int)iEvent.id().run();
  event_ = (int)iEvent.id().event();


    //David's functionality
    
  std::cout<<"What what"<<std::endl;
  using namespace reco;

  Handle <ME0MuonCollection > OurMuons;
  iEvent.getByLabel <ME0MuonCollection> ("me0SegmentMatcher", OurMuons);
  
  unsigned int recosize=OurMuons->size();

  std::cout<<recosize<<std::endl;cks

  
}

void TestAnalyzer::endJob() {}

DEFINE_FWK_MODULE(TestAnalyzer);
