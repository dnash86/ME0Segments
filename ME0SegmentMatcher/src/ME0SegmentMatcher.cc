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

ME0SegmentMatcher::ME0SegmentMatcher(const edm::ParameterSet& pas) : iev(0){
	
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
    Sigma_Phi_Total_h = new TH1F("Sigma_Phi_Total_h"      , "#sigma #phi Total"   , 100, 0., 0.08 );
    Sigma_Phi_Segment_h = new TH1F("Sigma_Phi_Segment_h"      , "#sigma #phi Segment"   , 100, 0., 0.01 );
    Sigma_Phi_Track_h = new TH1F("Sigma_Phi_Track_h"      , "#sigma #phi Track"   , 100, 0., 0.08 );
    Muon_Eta = new TH1F("Muon_Eta"      , "Muon Eta"   , 200, -4.5, 4.5 );
    Muon_Pt = new TH1F("Muon_Pt"      , "Muon Pt"   , 200, 0, 100 );
    FailMuon_Pt_HighEta = new TH1F("FailMuon_Pt_HighEta"      , "FailMuon Pt"   , 200, 0, 100 );
    FailMuon_Pt = new TH1F("FailMuon_Pt"      , "FailMuon Pt"   , 200, 0, 100 );
    FailMuon_Eta = new TH1F("FailMuon_Eta"      , "FailMuon Eta"   , 200, -4.5, 4.5 );
    Muon_Pt_HighEta = new TH1F("Muon_Pt_HighEta"      , "Muon Pt"   , 200, 0, 100 );
    Muon_OuterEta = new TH1F("Muon_OuterEta"      , "Muon Pt"   , 200, -4.5, 4.5 );
    MatchedSeg_Eta = new TH1F("MatchedSeg_Eta"      , "MatchedSeg Eta"   , 200, -4.5, 4.5 );
    Segment_Eta = new TH1F("Segment_Eta"      , "Segment Eta"   , 200, -4.5, 4.5 );
    Track_Eta = new TH1F("Track_Eta"      , "Track Eta"   , 200, -4.5, 4.5 );
    PreLoopTrack_Eta = new TH1F("PreLoopTrack_Eta"      , "PreLoopTrack Eta"   , 200, -4.5, 4.5 );
    Sigma_R_Total_prof = new TProfile("Sigma_R_Total_prof"      , "#sigma R Total"   ,     8, 2.4, 4.0, 0., 25, "s");
    Sigma_R_Segment_prof = new TProfile("Sigma_R_Segment_prof"      , "#sigma R Segment"   ,     8, 2.4, 4.0, 0., 2.5, "s");
    Sigma_R_Track_prof = new TProfile("Sigma_R_Track_prof"      , "#sigma R Track"   ,     8, 2.4, 4.0, 0., 25, "s");
    Sigma_Phi_Total_prof = new TProfile("Sigma_Phi_Total_prof"      , "#sigma #phi Total"  ,     8, 2.4, 4.0 , 0., 0.08, "s");
    Sigma_Phi_Segment_prof = new TProfile("Sigma_Phi_Segment_prof"      , "#sigma #phi Segment"  ,     8, 2.4, 4.0 , 0., 0.01, "s");
    Sigma_Phi_Track_prof = new TProfile("Sigma_Phi_Track_prof"      , "#sigma #phi Track"  ,     8, 2.4, 4.0 , 0., 0.08, "s");
    //NumTracks = new int;
    //NumParticles = new int;
    NumTracks = 0;
    NumParticles = 0;
    NumMatches = 0;
  }

}

ME0SegmentMatcher::~ME0SegmentMatcher() {

  //LogDebug("ME0SegmentMatcher") << "deleting ME0SegmentMatcher after " << iev << " events w/csc data.";
    //delete segmentBuilder_;
  if (debug_){
    //std::cout<<"On the destructor now"<<std::endl;
    std::cout<<NumParticles<<" gen particles, "<<NumTracks<<" tracks"<<NumMatches<<" matches"<<std::endl;
    HistFile->cd();
    Sigma_R_Total_h->Write();    Sigma_Phi_Total_h->Write();
    Sigma_R_Segment_h->Write();    Sigma_Phi_Segment_h->Write();
    Sigma_R_Track_h->Write();    Sigma_Phi_Track_h->Write();
    Muon_Eta->Write(); Muon_Pt->Write();   Muon_Pt_HighEta->Write();    Track_Eta->Write(); MatchedSeg_Eta->Write();
    PreLoopTrack_Eta->Write(); Segment_Eta->Write();
    Muon_OuterEta->Write();
    FailMuon_Pt->Write(); FailMuon_Pt_HighEta->Write(); FailMuon_Eta->Write();
    Sigma_R_Total_prof->Write();    Sigma_Phi_Total_prof->Write();
    Sigma_R_Segment_prof->Write();    Sigma_Phi_Segment_prof->Write();
    Sigma_R_Track_prof->Write();    Sigma_Phi_Track_prof->Write();
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
    std::vector<ME0Muon> TempStore; 

    NumTracks += generalTracks->size();
    NumParticles += OurSegments->size();

    int TrackNumber = 0;
    std::vector<int> TkIndices, TkIndex, TkToKeep;
    for (std::vector<Track>::const_iterator thisTrack = generalTracks->begin();
	 thisTrack != generalTracks->end(); ++thisTrack,++TrackNumber){
      //Initializing our plane
      //if ( (abs(thisTrack->eta()) > 4.0) || (abs(thisTrack->eta()) < 2.4) ) continue;
      TkIndices.push_back(TrackNumber);
      float zSign  = thisTrack->pz()/fabs(thisTrack->pz());
      float zValue = 560. * zSign;
      Plane *plane = new Plane(Surface::PositionType(0,0,zValue),Surface::RotationType());
      //Getting the initial variables for propagation
      int chargeReco = thisTrack->charge(); 
      GlobalVector p3reco, r3reco;

      // p3reco = GlobalVector(thisTrack->px(), thisTrack->py(), thisTrack->pz());
      // r3reco = GlobalVector(thisTrack->vx(), thisTrack->vy(), thisTrack->vz());

      p3reco = GlobalVector(thisTrack->outerPx(), thisTrack->outerPy(), thisTrack->outerPz());
      r3reco = GlobalVector(thisTrack->outerX(), thisTrack->outerY(), thisTrack->outerZ());

      AlgebraicSymMatrix66 covReco;
        //This is to fill the cov matrix correctly
        AlgebraicSymMatrix55 covReco_curv;
        //covReco_curv = thisTrack->covariance();
        covReco_curv = thisTrack->outerStateCovariance();
        FreeTrajectoryState initrecostate = getFTS(p3reco, r3reco, chargeReco, covReco_curv, &*bField);
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
      GlobalVector p3FinalReco, r3FinalReco;
      getFromFTS(finalrecostate, p3FinalReco, r3FinalReco, chargeReco, covFinalReco);
    
      //Track_Eta->Fill(thisTrack->eta());      
      //Track_Eta->Fill(thisTrack->outerEta());      
      //Now we compare for each me0 segment:
      int SegmentNumber = 0;


      PreLoopTrack_Eta->Fill(thisTrack->outerEta());      
      for (std::vector<ME0Segment>::const_iterator thisSegment = OurSegments->begin();
	   thisSegment != OurSegments->end(); ++thisSegment,++SegmentNumber){
	if (SegmentNumber==0) Track_Eta->Fill(thisTrack->outerEta()); 
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
	//std::cout<<rho_hit<<", "<<rho_track<<std::endl;
	//std::cout<<phi_hit<<", "<<phi_track<<std::endl;
	bool R_MatchFound = false, Phi_MatchFound = false;
	if ( zSign * thisPosition.z() > 0 ) {                          
	  if ( fabs(rho_hit-rho_track) < 3.0 * sigmarho) R_MatchFound = true;
	  if ( fabs(phi_hit-phi_track) < 3.0 * sigmaphi) Phi_MatchFound = true;
	}

	if (R_MatchFound && Phi_MatchFound)
	  {
	    if (abs(thisTrack->outerEta()) < 1.0){
	      std::cout<<" A match was found with this track and segment info: "<<std::endl;
	      // std::cout<<r3FinalReco.perp()<<", "<<r3FinalReco.phi()<<std::endl;
	      // std::cout<<thisPosition.perp()<<", "<<thisPosition.phi()<<std::endl;
	      std::cout<<rho_hit<<", "<<rho_track<<std::endl;
	      std::cout<<phi_hit<<", "<<phi_track<<std::endl;

	      std::cout<<sigmarho_hit<<", "<<sigmarho_track<<std::endl;
	      std::cout<<sigmaphi_hit<<", "<<sigmaphi_track<<std::endl;

	      std::cout<<thisTrack->outerPt()<<", "<<thisTrack->outerEta()<<std::endl;
	      std::cout<<"Track location: "<<r3FinalReco.perp()<<", "<<r3FinalReco.eta()<<std::endl;
	      //std::cout<<thisPosition.pt()<<", "<<thisPosition.eta()<<std::endl;
	    }
	    NumMatches++;

	    Sigma_R_Total_prof->Fill(thisPosition.eta(),sigmarho);    Sigma_Phi_Total_prof->Fill(thisPosition.eta(),sigmaphi);
	    Sigma_R_Segment_prof->Fill(thisPosition.eta(),sigmarho_hit);    Sigma_Phi_Segment_prof->Fill(thisPosition.eta(),sigmaphi_hit);
	    Sigma_R_Track_prof->Fill(thisPosition.eta(),sigmarho_track);    Sigma_Phi_Track_prof->Fill(thisPosition.eta(),sigmaphi_track);
	    TrackRef thisTrackRef(generalTracks,TrackNumber);
	    ME0SegmentRef thisME0SegmentRef(OurSegments,SegmentNumber);
	    TempStore.push_back(reco::ME0Muon(thisTrackRef,thisME0SegmentRef));
	    TkIndex.push_back(TrackNumber);
	    
	    
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
	else {
	  if (abs(thisTrack->outerEta()) < 1.0) {FailMuon_Pt->Fill(thisTrack->outerPt()); }
	  //else Muon_Pt_HighEta->Fill(thisTrack->pt());
	  else FailMuon_Pt_HighEta->Fill(thisTrack->outerPt());
	  FailMuon_Eta->Fill(thisTrack->outerEta());
	}

      }

    }
    
    // for (std::vector<int>::const_iterator ReferenceMuonNumber = TkIndices.begin();
    // 	 ReferenceMuonNumber != TkIndices.end(); ++ReferenceMuonNumber){

    for (unsigned int i = 0; i < TkIndices.size(); ++i){
      int ReferenceMuonNumber = TkIndices[i];          // The muon number of the track, starts at 0 and increments
      double RefDelR = 99999.9, ComparisonIndex = 0;
      int WhichTrackToKeep=-1;
      for (std::vector<ME0Muon>::const_iterator thisMuon = TempStore.begin();
	   thisMuon != TempStore.end(); ++thisMuon, ++ComparisonIndex){
	
	int thisMuonNumber = TkIndex[ComparisonIndex];    //The track number of the muon we are currently looking at
	if (thisMuonNumber == ReferenceMuonNumber){        //This means we're looking at one track

	  ME0SegmentRef SegRef = thisMuon->me0segment();
	  TrackRef TkRef = thisMuon->innerTrack();
	  LocalPoint SegPos(SegRef->localPosition().x(),SegRef->localPosition().y(),SegRef->localPosition().z());
	  LocalPoint TkPos(TkRef->vx(),TkRef->vy(),TkRef->vz());
	  double delR = reco::deltaR(SegPos,TkPos);

	  if (delR < RefDelR) WhichTrackToKeep = ComparisonIndex;  //Storing a list of the vector indices of tracks to keep
	                                                           //Note: These are not the same as the "Track Numbers"
	}
      }
      if (WhichTrackToKeep != -1) TkToKeep.push_back(WhichTrackToKeep);
    }

    // for (std::vector<int>::const_iterator thisKeepIndex = TkToKeep.begin();
    // 	 thisKeepIndex != TkToKeep.end(); ++thisKeepIndex){
    for (unsigned int i = 0; i < TkToKeep.size(); ++i){
      int thisKeepIndex = TkToKeep[i];
      oc->push_back(TempStore[thisKeepIndex]);    //Maybe we need a 'new' here?


      //Debug plotting stuff
      ME0SegmentRef thisSegment= TempStore[thisKeepIndex].me0segment();
      TrackRef thisTrack = TempStore[thisKeepIndex].innerTrack();
      GlobalPoint thisPosition(TempStore[thisKeepIndex].me0segment()->localPosition().x(),TempStore[thisKeepIndex].me0segment()->localPosition().y(),TempStore[thisKeepIndex].me0segment()->localPosition().z());

      MatchedSeg_Eta->Fill(thisPosition.eta());
      Muon_Eta->Fill(thisTrack->outerEta());

      if (abs(thisTrack->outerEta()) < 1.0) {Muon_Pt->Fill(thisTrack->outerPt()); Muon_OuterEta->Fill(thisTrack->outerEta());}
      else Muon_Pt_HighEta->Fill(thisTrack->outerPt());


      Segment_Eta->Fill(thisPosition.eta());

    }
  	// fill the collection
    //segmentBuilder_->build(cscRecHits.product(), *oc); //@@ FILL oc

    // put collection in event
    ev.put(oc);

}

FreeTrajectoryState
ME0SegmentMatcher::getFTS(const GlobalVector& p3, const GlobalVector& r3, 
			   int charge, const AlgebraicSymMatrix55& cov,
			   const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CurvilinearTrajectoryError tCov(cov);
  
  return cov.kRows == 5 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

FreeTrajectoryState
ME0SegmentMatcher::getFTS(const GlobalVector& p3, const GlobalVector& r3, 
			   int charge, const AlgebraicSymMatrix66& cov,
			   const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CartesianTrajectoryError tCov(cov);
  
  return cov.kRows == 6 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

void ME0SegmentMatcher::getFromFTS(const FreeTrajectoryState& fts,
				    GlobalVector& p3, GlobalVector& r3, 
				    int& charge, AlgebraicSymMatrix66& cov){
  GlobalVector p3GV = fts.momentum();
  GlobalPoint r3GP = fts.position();

  GlobalVector p3T(p3GV.x(), p3GV.y(), p3GV.z());
  GlobalVector r3T(r3GP.x(), r3GP.y(), r3GP.z());
  p3 = p3T;
  r3 = r3T;  //Yikes, was setting this to p3T instead of r3T!?!
  // p3.set(p3GV.x(), p3GV.y(), p3GV.z());
  // r3.set(r3GP.x(), r3GP.y(), r3GP.z());
  
  charge = fts.charge();
  cov = fts.hasError() ? fts.cartesianError().matrix() : AlgebraicSymMatrix66();

}


 DEFINE_FWK_MODULE(ME0SegmentMatcher);
