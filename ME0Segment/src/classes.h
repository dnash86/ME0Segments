
#include <ME0Reconstruction/ME0Segment/interface/ME0Segment.h>
//#include <DataFormats/ME0RecHit/interface/ME0SegmentCollection.h>

#include <DataFormats/Common/interface/Wrapper.h>
#include <vector>

namespace{ 
  struct dictionary {
    //ME0SegmentCollection seg;    
    //edm::Wrapper<ME0SegmentCollection> dwc1;
    ME0Segment seg;
    //testingtoseeifitwillfail
    std::vector<ME0Segment> segs;
    edm::Wrapper< std::vector<ME0Segment> > dwc1;
  };
}
