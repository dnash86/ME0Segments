#ifndef ME0Segment_ME0SegmentMatcher_h
#define ME0Segment_ME0SegmentMatcher_h

/** \class ME0SegmentMatcher 
 * Produces a collection of ME0Segment's in endcap muon ME0s. 
 *
 * $Date: 2010/03/11 23:48:11 $
 *
 * \author David Nash
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"
#include "TH1.h" 
#include "TFile.h"
#include <TProfile.h>
//============================


#include <ME0Reconstruction/ME0Segment/interface/ME0Segment.h>

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
//#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TLorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TRandom3.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"

#include <math.h>

//==========================

//class ME0SegmentBuilder; 

class FreeTrajectoryState;
class MagneticField;
//class TRandom3;
class ME0SegmentMatcher : public edm::EDProducer {
public:
    /// Constructor
    explicit ME0SegmentMatcher(const edm::ParameterSet&);
    /// Destructor
    ~ME0SegmentMatcher();
    /// Produce the ME0Segment collection
    virtual void produce(edm::Event&, const edm::EventSetup&);

    FreeTrajectoryState getFTS(const GlobalVector& , const GlobalVector& , 
				   int , const AlgebraicSymMatrix66& ,
				   const MagneticField* );

    FreeTrajectoryState getFTS(const GlobalVector& , const GlobalVector& , 
				   int , const AlgebraicSymMatrix55& ,
				   const MagneticField* );

    void getFromFTS(const FreeTrajectoryState& ,
		  GlobalVector& , GlobalVector& , 
		  int& , AlgebraicSymMatrix66& );

private:
    //TRandom3 * Rand;
    int iev; // events through
    //edm::InputTag inputObjectsTag; // input tag labelling rechits for input
    //ME0SegmentBuilder* segmentBuilder_;
    bool debug_;
    TFile* HistFile;
    TH1F *Sigma_R_Total_h; TH1F *Sigma_Phi_Total_h;
    TH1F *Sigma_R_Segment_h; TH1F *Sigma_Phi_Segment_h;
    TH1F *Sigma_R_Track_h; TH1F *Sigma_Phi_Track_h;
    TH1F *Muon_Eta;   TH1F *Muon_Pt;   TH1F *Muon_Pt_HighEta;    TH1F *Track_Eta;    TH1F *MatchedSeg_Eta;
    TH1F *Muon_OuterEta;
    TH1F *FailMuon_Pt; TH1F *FailMuon_Pt_HighEta; TH1F *FailMuon_Eta;
    TH1F *PreLoopTrack_Eta; TH1F *Segment_Eta;
    int NumTracks, NumParticles, NumMatches;
    TProfile *Sigma_R_Total_prof; TProfile *Sigma_Phi_Total_prof;
    TProfile *Sigma_R_Segment_prof; TProfile *Sigma_Phi_Segment_prof;
    TProfile *Sigma_R_Track_prof; TProfile *Sigma_Phi_Track_prof;

};

#endif
