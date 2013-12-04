
#include <ME0Reconstruction/ME0Segment/interface/ME0Segment.h>
#include <ME0Reconstruction/ME0Segment/interface/ME0Muon.h>

#include <ME0Reconstruction/ME0Segment/interface/ME0SegmentCollection.h>
#include <ME0Reconstruction/ME0Segment/interface/ME0MuonCollection.h>


#include <DataFormats/Common/interface/Wrapper.h>
#include <vector>

namespace{ 
  struct dictionary {
    ME0Segment seg;
    std::vector<ME0Segment> segs;
    edm::Wrapper< std::vector<ME0Segment> > dwc1;
    
    reco::ME0Muon muon;
    std::vector<reco::ME0Muon> muons;
    edm::Wrapper< std::vector<reco::ME0Muon> > dwc2;

    ME0MuonCollection muoncol;
    edm::Wrapper<ME0MuonCollection> mcw1;
    edm::Ref<ME0MuonCollection> mcr1;

    ME0SegmentCollection segcol;
    edm::Wrapper<ME0SegmentCollection> scw1;
    edm::Ref<ME0SegmentCollection> scr1;
  };
}
