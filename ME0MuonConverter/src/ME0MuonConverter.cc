/** \file ME0MuonConverter.cc
 *
 * \author David Nash
 */

#include <ME0Reconstruction/ME0MuonConverter/src/ME0MuonConverter.h>
//#include <RecoLocalMuon/ME0Segment/src/ME0SegmentBuilder.h>

#include <FWCore/PluginManager/interface/ModuleDef.h>
#include <FWCore/Framework/interface/MakerMacros.h>

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

ME0MuonConverter::ME0MuonConverter(const edm::ParameterSet& pas) : iev(0) {
	
  produces<std::vector<reco::ME0Muon> >();  //May have to later change this to something that makes more sense, OwnVector, RefVector, etc

}

ME0MuonConverter::~ME0MuonConverter() {

}

void ME0MuonConverter::produce(edm::Event& ev, const edm::EventSetup& setup) {

  using namespace edm;

  using namespace reco;

  Handle <std::vector<ME0Muon> > OurMuons;
  ev.getByLabel <std::vector<ME0Muon> > ("me0SegmentMatcher", OurMuons);
  
  //std::auto_ptr<std::vector<RecoChargedCandidate> > oc( new std::vector<RecoChargedCandidate> ); 

  for (std::vector<ME0Muon>::const_iterator thisMuon = OurMuons->begin();
       thisMuon != OurMuons->end(); ++thisMuon){
    std::cout<<"On a muon:"<<std::endl;
    std::cout<<thisMuon->pt()<<std::endl;
  }
    
  //ev.put(oc);
}


 DEFINE_FWK_MODULE(ME0MuonConverter);
