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


//#include "CLHEP/Matrix/SymMatrix.h"
//#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"
#include "TH1.h" 
#include "TFile.h"

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

    FreeTrajectoryState getFromCLHEP(const CLHEP::Hep3Vector& , const CLHEP::Hep3Vector& , 
				   int , const AlgebraicSymMatrix66& ,
				   const MagneticField* );

    FreeTrajectoryState getFromCLHEP(const CLHEP::Hep3Vector& , const CLHEP::Hep3Vector& , 
				   int , const AlgebraicSymMatrix55& ,
				   const MagneticField* );

    void getFromFTS(const FreeTrajectoryState& ,
		  CLHEP::Hep3Vector& , CLHEP::Hep3Vector& , 
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

};

#endif
