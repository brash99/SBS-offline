//////////////////////////////////////////////////////////////////////////
//
// SBSRPFarSideHodo
//
// General detector with TDC information
//
//////////////////////////////////////////////////////////////////////////

#ifndef SBSCDet_h
#define SBSCDet_h

#include "SBSGenericDetector.h"
#include "SBSCDet_Hit.h"
#include "TBits.h"
#include "TClonesArray.h"
#include <cstdint>
#include <map>
#include <vector>

class THaTrack;
class THaBenhmark;

class SBSCDet : public SBSGenericDetector {

public:
  SBSCDet( const char* name, const char* description, 
      THaApparatus* apparatus=0);

  virtual ~SBSCDet();

  // Overwritten derived functions
  virtual Int_t   ReadDatabase( const TDatime& date );
  virtual Int_t   DefineVariables( EMode mode = kDefine );
  virtual Int_t   Decode( const THaEvData& evdata );
  virtual Int_t   CoarseProcess( TClonesArray& tracks );
  virtual Int_t   FineProcess( TClonesArray& tracks );
  virtual Int_t   FindGoodHit(SBSElement *element);
  virtual void    Clear( Option_t* opt="" );

  SBSCDet_Hit*		GetHit(Int_t i) const
  {return (SBSCDet_Hit*)fHits->At(i); }

  Int_t 		GetNumHits() const
    { return fHits->GetLast()+1; }

  Int_t                 GetNumPulseCandidates() const
    { return static_cast<Int_t>(fPulsePMT.size()); }
  Int_t                 GetNumLeadingOnlySlots() const
    { return fNLeadingOnlySlots; }
  Int_t                 GetNumTrailingOnlySlots() const
    { return fNTrailingOnlySlots; }
  Int_t                 GetNumInvalidPairSlots() const
    { return fNInvalidPairSlots; }
  Int_t                 GetNumLayerPairs() const
    { return static_cast<Int_t>(fPairIndex.size()); }

  /// Retain every complete decoded CDet TDC pulse (LE, TE, positive ToT) with
  /// stable channel/pulse identity.
  /// Disabled by default so existing replay output is unchanged.
  void                  SetStorePulseCandidates(Bool_t enable = true)
    { fStorePulseCandidates = enable; }

  void ApplyECalTimingCalibration(Double_t ecalTime, Int_t ecalClusterIndex,
                                  Double_t ecalX, Double_t ecalY,
                                  Double_t ecalZ);
  
  //Int_t                GetNumClusters() const
  //  { return fClusters->GetLast()+1; }

protected:

  TClonesArray*		fHits;		// Array of hits for each event
  Double_t	fHit_tmin;
  Double_t 	fHit_tmax;
  Double_t	fHit_totmin;
  Double_t 	fHit_totmax;

  Bool_t        fStorePulseCandidates;
  Int_t         fNLeadingOnlySlots;
  Int_t         fNTrailingOnlySlots;
  Int_t         fNInvalidPairSlots;
  std::vector<Int_t>    fPulsePMT;
  std::vector<Int_t>    fPulseIndex;
  std::vector<Int_t>    fPulseLEIndex;
  std::vector<Int_t>    fPulseTEIndex;
  std::vector<Int_t>    fPulseRow;
  std::vector<Int_t>    fPulseCol;
  std::vector<Int_t>    fPulseLayer;
  std::vector<Double_t> fPulseX;
  std::vector<Double_t> fPulseY;
  std::vector<Double_t> fPulseZ;
  std::vector<Double_t> fPulseLE;
  std::vector<Double_t> fPulseTE;
  std::vector<Double_t> fPulseToT;
  std::vector<Double_t> fPulseLERaw;
  std::vector<Double_t> fPulseTERaw;
  std::vector<Double_t> fPulseToTRaw;
  std::vector<Double_t> fPulseLECorrected;
  std::vector<Double_t> fPulseTECorrected;
  std::vector<Double_t> fPulseToTNs;
  std::vector<Double_t> fPulseECalResidual;
  std::vector<Int_t>    fPulseCalibrationValid;
  std::vector<Double_t> fPulseCorrectedX;
  std::vector<Double_t> fPulseProjectedECalX;
  std::vector<Double_t> fPulseProjectedECalY;
  std::vector<Double_t> fPulseECalXResidual;
  std::vector<Double_t> fPulseECalYResidual;
  std::vector<Int_t>    fPulseBroadQualityPass;
  std::vector<Int_t>    fPulseECalEligibilityPass;
  std::vector<Int_t>    fPulseSpatialPass;

  std::vector<Int_t>    fPairIndex;
  std::vector<Int_t>    fPairPulseIndexL1;
  std::vector<Int_t>    fPairPulseIndexL2;
  std::vector<Int_t>    fPairPMTL1;
  std::vector<Int_t>    fPairPMTL2;
  std::vector<Double_t> fPairTimeL1;
  std::vector<Double_t> fPairTimeL2;
  std::vector<Double_t> fPairTimeMean;
  std::vector<Double_t> fPairDeltaTime;
  std::vector<Double_t> fPairDeltaX;
  std::vector<Double_t> fPairDeltaY;
  std::vector<Double_t> fPairScore;
  std::vector<Double_t> fPairECalResidual;

  Bool_t fTimingCalibrationEnabled;
  Bool_t fTimingCalibrationLoaded;
  Int_t fTimingStatus;
  Int_t fTimingECalClusterIndex;
  Double_t fTimingECalTime;
  Double_t fTimingTDCToNs;
  std::vector<Double_t> fTimingPixelOffset;
  Double_t fTimingECalP0;
  Double_t fTimingECalP1;
  Double_t fTimingECalDelta;
  Double_t fTimingWalkP1L1;
  Double_t fTimingWalkP1L2;
  Double_t fTimingWalkToTRefL1;
  Double_t fTimingWalkToTRefL2;
  Double_t fTimingShift;
  Bool_t fSelectionEnabled;
  Double_t fSelectionLEMin;
  Double_t fSelectionLEMax;
  Double_t fSelectionToTMin;
  Double_t fSelectionToTMax;
  Double_t fSelectionECalTimeMin;
  Double_t fSelectionECalTimeMax;
  Double_t fSelectionECalXMin;
  Double_t fSelectionECalXMax;
  Double_t fSelectionECalYMin;
  Double_t fSelectionECalYMax;
  Double_t fSelectionXScaleL1;
  Double_t fSelectionXScaleL2;
  Double_t fSelectionXAlignmentL1;
  Double_t fSelectionXAlignmentL2;
  Double_t fSelectionXResidualOffset;
  Double_t fSelectionXResidualMax;
  Double_t fSelectionYResidualOffset;
  Double_t fSelectionYResidualMax;
  Bool_t fPairingEnabled;
  Double_t fPairingDeltaTimeCenter;
  Double_t fPairingDeltaTimeMax;
  Double_t fPairingDeltaXMax;
  Double_t fPairingDeltaYMax;
  Double_t fPairingTimeScale;
  Double_t fPairingXScale;
  Bool_t fPairingAllowMultiple;

  void BuildPulseCandidates();
  void BuildLayerPairs();
  void ClearPulseCandidates();
  Int_t GetTimingStatus() const { return fTimingStatus; }
  Int_t GetTimingECalClusterIndex() const { return fTimingECalClusterIndex; }
  Double_t GetTimingECalTime() const { return fTimingECalTime; }

  ClassDef(SBSCDet,9)  // Describes scintillator plane with F1TDC as a detector
};

#endif
