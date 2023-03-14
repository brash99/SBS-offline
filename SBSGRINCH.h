#ifndef ROOT_SBSGRINCH
#define ROOT_SBSGRINCH

//////////////////////////////////////////////////////////////////////////
//
// SBSGRINCH
//
// The Hall A RICH
//
//////////////////////////////////////////////////////////////////////////

//#include "THaPidDetector.h"
#include "SBSCherenkovDetector.h"
#include "SBSCherenkov_ClusterList.h"
#include "TBits.h"
#include "TClonesArray.h"
#include <stdint.h>
#include <map>

class THaTrack;
class THaBenchmark;

//class SBSGRINCH : public THaPidDetector {
class SBSGRINCH : public SBSCherenkovDetector {
  
public:

  explicit SBSGRINCH( const char* name, const char* description="",
	   THaApparatus* apparatus=nullptr );
  virtual ~SBSGRINCH();
  
  virtual void         Clear( Option_t* opt="" );
  virtual Int_t        Decode( const THaEvData& );
  virtual Int_t        CoarseProcess( TClonesArray& tracks );
  virtual Int_t        FineProcess( TClonesArray& tracks );

protected:

  //Double_t fTrackX;      // x pos of Golden Track in RICH plane
  //Double_t fTrackY;      // y pos of Golden Track in RICH plane

  Double_t fMaxSep; // Max separation between PMT and another PMT to count as "neighbors"
  Double_t fMaxSep2; //square of fMaxSep

  Double_t fTrackMatchXslope;
  Double_t fTrackMatchX0;
  Double_t fTrackMatchXsigma; 
  Double_t fTrackMatchYslope;
  Double_t fTrackMatchY0;
  Double_t fTrackMatchYsigma;
  Double_t fTrackMatchNsigmaCut;

  virtual Int_t   FindClusters();
  virtual Int_t   MatchClustersWithTracks( TClonesArray& tracks );
  virtual Int_t   SelectBestCluster(Int_t nmatch=0);

  //  Int_t fBestClusterIndex; //the biggest cluster with a track match, if any matched clusters are found. Perhaps this should go with SBSCherenkov_Detector rather than SBSGRINCH

  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );
  
  ClassDef(SBSGRINCH,0)   //The Hall A RICH
};

#endif









