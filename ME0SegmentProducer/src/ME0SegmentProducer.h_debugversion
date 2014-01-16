#ifndef ME0Segment_ME0SegmentProducer_h
#define ME0Segment_ME0SegmentProducer_h

/** \class ME0SegmentProducer 
 * Produces a collection of ME0Segment's in endcap muon ME0s. 
 *
 * $Date: 2010/03/11 23:48:11 $
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
class TRandom3;
//class AlgebraicMatrix;
class ME0SegmentProducer : public edm::EDProducer {
public:
    /// Constructor
    explicit ME0SegmentProducer(const edm::ParameterSet&);
    /// Destructor
    ~ME0SegmentProducer();
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

    void RotateCovMatrix(const AlgebraicMatrix&, const AlgebraicSymMatrix&,int, AlgebraicSymMatrix&);


private:
    TRandom3 * Rand;
    int iev; // events through
    //edm::InputTag inputObjectsTag; // input tag labelling rechits for input
    //ME0SegmentBuilder* segmentBuilder_;
    int NumSegs;
};

#endif
