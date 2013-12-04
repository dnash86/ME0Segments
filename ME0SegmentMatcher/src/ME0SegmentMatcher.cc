/** \file ME0SegmentMatcher.cc
 *
 * \author David Nash
 */

#include <ME0Reconstruction/ME0SegmentMatcher/src/ME0SegmentMatcher.h>
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

ME0SegmentMatcher::ME0SegmentMatcher(const edm::ParameterSet& pas) : iev(0) {
	
  //inputObjectsTag = pas.getParameter<edm::InputTag>("inputObjects");
    //segmentBuilder_ = new ME0SegmentBuilder(pas); // pass on the PS
  //Rand = new TRandom3();
  	// register what this produces
  //produces<ME0MuonCollection>();
  produces<std::vector<reco::ME0Muon> >();  //May have to later change this to something that makes more sense, OwnVector, RefVector, etc

    //Put what we produce here, obviously not what's listed below
  debug_ = pas.getParameter<bool>("debug");
    //produces<std::vector<ME0Segment> >();  
  if (debug_){
    HistFile = new TFile(pas.getParameter<std::string>("DebugHistos").c_str(), "recreate");
    Sigma_R_Total_h = new TH1F("Sigma_R_Total_h"      , "#sigma R Total"   , 100, 0., 25 );
    Sigma_R_Segment_h = new TH1F("Sigma_R_Segment_h"      , "#sigma R Segment"   , 100, 0., 2.5 );
    Sigma_R_Track_h = new TH1F("Sigma_R_Track_h"      , "#sigma R Track"   , 100, 0., 25 );
    Sigma_Phi_Total_h = new TH1F("Sigma_Phi_Total_h"      , "#sigma #phi Total"   , 100, 0., 0.05 );
    Sigma_Phi_Segment_h = new TH1F("Sigma_Phi_Segment_h"      , "#sigma #phi Segment"   , 100, 0., 0.005 );
    Sigma_Phi_Track_h = new TH1F("Sigma_Phi_Track_h"      , "#sigma #phi Track"   , 100, 0., 0.05 );
  }

}

ME0SegmentMatcher::~ME0SegmentMatcher() {

  //LogDebug("ME0SegmentMatcher") << "deleting ME0SegmentMatcher after " << iev << " events w/csc data.";
    //delete segmentBuilder_;
  if (debug_){
    //std::cout<<"On the destructor now"<<std::endl;
    HistFile->cd();
    Sigma_R_Total_h->Write();    Sigma_Phi_Total_h->Write();
    Sigma_R_Segment_h->Write();    Sigma_Phi_Segment_h->Write();
    Sigma_R_Track_h->Write();    Sigma_Phi_Track_h->Write();
    delete HistFile; HistFile = 0;
  }
}

void ME0SegmentMatcher::produce(edm::Event& ev, const edm::EventSetup& setup) {

    LogDebug("ME0SegmentMatcher") << "start producing segments for " << ++iev << "th event ";
	
    // find the geometry (& conditions?) for this event & cache it in the builder
  
    // edm::ESHandle<ME0Geometry> h;
    // setup.get<MuonGeometryRecord>().get(h);
    // const ME0Geometry* pgeom = &*h;


    //segmentBuilder_->setGeometry(pgeom);
	
    // get the collection of ME0RecHit2D

    // edm::Handle<ME0RecHit2DCollection> cscRecHits;
    // ev.getByLabel(inputObjectsTag, cscRecHits);  

    // create empty collection of Segments
    //std::auto_ptr<ME0SegmentCollection> oc( new ME0SegmentCollection );       //Change to collection

    //Do the stuff, put the smeared stuff in to oc

    //Getting the objects we'll need
    
    using namespace edm;
    ESHandle<MagneticField> bField;
    setup.get<IdealMagneticFieldRecord>().get(bField);
    ESHandle<Propagator> shProp;
    setup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", shProp);

    
    using namespace reco;

    Handle<std::vector<ME0Segment> > OurSegments;
    ev.getByLabel<std::vector<ME0Segment> >("me0SegmentProducer", OurSegments);

    Handle <TrackCollection > generalTracks;
    ev.getByLabel <TrackCollection> ("generalTracks", generalTracks);

    std::auto_ptr<std::vector<ME0Muon> > oc( new std::vector<ME0Muon> ); 

    int TrackNumber = 0;
    for (std::vector<Track>::const_iterator thisTrack = generalTracks->begin();
	 thisTrack != generalTracks->end(); ++thisTrack,++TrackNumber){
      //Initializing our plane
      float zSign  = thisTrack->pz()/fabs(thisTrack->pz());
      float zValue = 560. * zSign;
      Plane *plane = new Plane(Surface::PositionType(0,0,zValue),Surface::RotationType());
      //Getting the initial variables for propagation
      int chargeReco = thisTrack->charge(); 
      CLHEP::Hep3Vector p3reco, r3reco;
      p3reco = CLHEP::Hep3Vector(thisTrack->px(), thisTrack->py(), thisTrack->pz());
      r3reco = CLHEP::Hep3Vector(thisTrack->vx(), thisTrack->vy(), thisTrack->vz());
      AlgebraicSymMatrix66 covReco;
        //This is to fill the cov matrix correctly
        AlgebraicSymMatrix55 covReco_curv;
        covReco_curv = thisTrack->covariance();
        FreeTrajectoryState initrecostate = getFromCLHEP(p3reco, r3reco, chargeReco, covReco_curv, &*bField);
        getFromFTS(initrecostate, p3reco, r3reco, chargeReco, covReco);

      //Now we propagate and get the propagated variables from the propagated state
      SteppingHelixStateInfo startrecostate(initrecostate);
      SteppingHelixStateInfo lastrecostate;

      const SteppingHelixPropagator* ThisshProp = 
	dynamic_cast<const SteppingHelixPropagator*>(&*shProp);
	
      lastrecostate = ThisshProp->propagate(startrecostate, *plane);
	
      FreeTrajectoryState finalrecostate;
      lastrecostate.getFreeState(finalrecostate);

      AlgebraicSymMatrix66 covFinalReco;
      CLHEP::Hep3Vector p3FinalReco, r3FinalReco;
      getFromFTS(finalrecostate, p3FinalReco, r3FinalReco, chargeReco, covFinalReco);
    

      //Now we compare for each me0 segment:
      int SegmentNumber = 0;
      for (std::vector<ME0Segment>::const_iterator thisSegment = OurSegments->begin();
	   thisSegment != OurSegments->end(); ++thisSegment,++SegmentNumber){
      
	//ME0Segments actually have globally initialized positions and directions, so lets cast them as global points and vectors
	GlobalPoint thisPosition(thisSegment->localPosition().x(),thisSegment->localPosition().y(),thisSegment->localPosition().z());
	GlobalVector thisDirection(thisSegment->localDirection().x(),thisSegment->localDirection().y(),thisSegment->localDirection().z());
	//The same goes for the error
	AlgebraicMatrix thisCov(4,4,0);    //Load it with parametersError, is there a better way to assign matrices?
	for (int i = 1; i <=4; i++){
	  for (int j = 1; j <=4; j++){
	    thisCov(i,j) = thisSegment->parametersError()(i,j);
	    //std::cout<<thisSegment->parametersError()(i,j)<<"   ";
	  }
	  //std::cout<<std::endl;
	}

	//Computing the sigma for the track
	Double_t rho_track = r3FinalReco.perp();
	Double_t phi_track = r3FinalReco.phi();

	Double_t drhodx_track = r3FinalReco.x()/rho_track;
	Double_t drhody_track = r3FinalReco.y()/rho_track;
	Double_t dphidx_track = -r3FinalReco.y()/(rho_track*rho_track);
	Double_t dphidy_track = r3FinalReco.x()/(rho_track*rho_track);
      
	Double_t sigmarho_track = sqrt( drhodx_track*drhodx_track*covFinalReco(0,0)+
					drhody_track*drhody_track*covFinalReco(1,1)+
					drhodx_track*drhody_track*2*covFinalReco(0,1) );
      
	Double_t sigmaphi_track = sqrt( dphidx_track*dphidx_track*covFinalReco(0,0)+
					dphidy_track*dphidy_track*covFinalReco(1,1)+
					dphidx_track*dphidy_track*2*covFinalReco(0,1) );

	//Computing the sigma for the hit
	Double_t rho_hit = thisPosition.perp();
	Double_t phi_hit = thisPosition.phi();

	Double_t drhodx_hit = thisPosition.x()/rho_hit;
	Double_t drhody_hit = thisPosition.y()/rho_hit;
	Double_t dphidx_hit = -thisPosition.y()/(rho_hit*rho_hit);
	Double_t dphidy_hit = thisPosition.x()/(rho_hit*rho_hit);
      
	Double_t sigmarho_hit = sqrt( drhodx_hit*drhodx_hit*thisCov(2,2)+
				      drhody_hit*drhody_hit*thisCov(3,3)+
				      drhodx_hit*drhody_hit*2*thisCov(2,3) );
      
	Double_t sigmaphi_hit = sqrt( dphidx_hit*dphidx_hit*thisCov(2,2)+
				      dphidy_hit*dphidy_hit*thisCov(3,3)+
				      dphidx_hit*dphidy_hit*2*thisCov(2,3) );

	//Adding the sigmas
	Double_t sigmarho = sqrt(sigmarho_track*sigmarho_track + sigmarho_hit*sigmarho_hit);
	Double_t sigmaphi = sqrt(sigmaphi_track*sigmaphi_track + sigmaphi_hit*sigmaphi_hit);

	//Filling histos for debugging purposes
	if (debug_){
	  Sigma_R_Total_h->Fill(sigmarho);    Sigma_Phi_Total_h->Fill(sigmaphi);
	  Sigma_R_Segment_h->Fill(sigmarho_hit);    Sigma_Phi_Segment_h->Fill(sigmaphi_hit);
	  Sigma_R_Track_h->Fill(sigmarho_track);    Sigma_Phi_Track_h->Fill(sigmaphi_track);
	}
	//Checking if there is a match in rho and in phi, assuming they are pointing in the same direction
	//Try making a histo of these errors to be sure they're small
	bool R_MatchFound = false, Phi_MatchFound = false;
	if ( zSign * thisPosition.z() > 0 ) {                          
	  if ( fabs(rho_hit-rho_track) < 3.0 * sigmarho) R_MatchFound = true;
	  if ( fabs(phi_hit-phi_track) < 3.0 * sigmaphi) Phi_MatchFound = true;
	}
      
	if (R_MatchFound && Phi_MatchFound)
	  {
	    //Stuff we do if a match is found

	    TrackRef thisTrackRef(generalTracks,TrackNumber);
	    ME0SegmentRef thisME0SegmentRef(OurSegments,SegmentNumber);

	    oc->push_back(reco::ME0Muon(thisTrackRef,thisME0SegmentRef));    //Maybe we need a 'new' here?


	    // std::cout<<"Found a Match:"<<std::endl;
	    // std::cout<<"Sig_R = "<<sigmarho<<", Sig_Phi = "<<sigmaphi<<std::endl;
	    // std::cout<<"Track = "<<r3FinalReco.x()<<","<<r3FinalReco.y()<<","<<r3FinalReco.z()<<std::endl;
	    // std::cout<<"Hit = "<<thisSegment->localPosition().x()<<","<<thisSegment->localPosition().y()<<","<<thisSegment->localPosition().z()<<std::endl;

	    // //track updating part, no clue if this is right yet...  
	    // //Modeled after http://cmslxr.fnal.gov/lxr/source/Alignment/CommonAlignmentProducer/plugins/GlobalTrackerMuonAlignment.cc#1385
	    // //------I think this part is declaring and initializing things for our fitter/updator...
	    // KFUpdator* theUpdator = new KFUpdator();
	    // Chi2MeasurementEstimator* theEstimator = new Chi2MeasurementEstimator(100000,100000); //Should check how this should be used
	    // const SteppingHelixPropagator* OurProp = 
	    //   dynamic_cast<const SteppingHelixPropagator*>(&*shProp);
	    // theFitter = new KFTrajectoryFitter(*OurProp, 
	    // 				       *theUpdator, 
	    // 				       *theEstimator);
	    // theSmoother = new KFTrajectorySmoother(*OurProp,
	    // 					   *theUpdator,     
	    // 					   *theEstimator);
	    // //------Now I think we should go on to declaring and initializing the things for updating..

	    // // I guess you need a trajectory state on surface?  Will it need to be persistent?  Maybe not..
	    // //If so, lets hope you can initialize this without a detId...  it should look like (initialTSOS, trackDetId.rawId());
	    // PTrajectoryStateOnDet  PTraj = 
	    //   trajectoryStateTransform::persistentState(initialTSOS);
	    // // We'll need a rechit and direction for the trajectory seed
	    // const TrajectorySeed seedT(PTraj, recHit, direction);
	    // trajVec = theFitter->fit(seedT, recHitMu, initialTSOS);
	  }

      }

    }
  	// fill the collection
    //segmentBuilder_->build(cscRecHits.product(), *oc); //@@ FILL oc

    // put collection in event
    ev.put(oc);
}

FreeTrajectoryState
ME0SegmentMatcher::getFromCLHEP(const CLHEP::Hep3Vector& p3, const CLHEP::Hep3Vector& r3, 
					      int charge, const AlgebraicSymMatrix55& cov,
					      const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CurvilinearTrajectoryError tCov(cov);
  
  return cov.kRows == 5 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

FreeTrajectoryState
ME0SegmentMatcher::getFromCLHEP(const CLHEP::Hep3Vector& p3, const CLHEP::Hep3Vector& r3, 
					      int charge, const AlgebraicSymMatrix66& cov,
					      const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CartesianTrajectoryError tCov(cov);
  
  return cov.kRows == 6 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

void ME0SegmentMatcher::getFromFTS(const FreeTrajectoryState& fts,
				    CLHEP::Hep3Vector& p3, CLHEP::Hep3Vector& r3, 
				    int& charge, AlgebraicSymMatrix66& cov){
  GlobalVector p3GV = fts.momentum();
  GlobalPoint r3GP = fts.position();

  p3.set(p3GV.x(), p3GV.y(), p3GV.z());
  r3.set(r3GP.x(), r3GP.y(), r3GP.z());
  
  charge = fts.charge();
  cov = fts.hasError() ? fts.cartesianError().matrix() : AlgebraicSymMatrix66();

}


 DEFINE_FWK_MODULE(ME0SegmentMatcher);
