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
  /// Read-only snapshot of one complete decoded CDet pulse candidate.
  /// Indices refer to the event-local pulse arrays and remain valid only for
  /// the event in which the snapshot was requested.
  struct PulseCandidate {
    Int_t pmt;
    Int_t index;
    Int_t leIndex;
    Int_t teIndex;
    Int_t row;
    Int_t col;
    Int_t layer;
    Double_t x;
    Double_t y;
    Double_t z;
    Double_t leadingEdge;
    Double_t trailingEdge;
    Double_t timeOverThreshold;
    Double_t rawLeadingEdge;
    Double_t rawTrailingEdge;
    Double_t rawTimeOverThreshold;
    Double_t correctedLeadingEdge;
    Double_t correctedTrailingEdge;
    Double_t timeOverThresholdNs;
    Double_t ecalResidual;
    Bool_t calibrationValid;
    Double_t correctedX;
    Double_t projectedECalX;
    Double_t projectedECalY;
    Double_t ecalXResidual;
    Double_t ecalYResidual;
    Bool_t broadQualityPass;
    Bool_t ecalEligibilityPass;
    Bool_t spatialPass;
  };

  /// Read-only snapshot of one accepted, one-to-one Layer-1/Layer-2 pair.
  struct LayerPair {
    Int_t index;
    Int_t pulseIndexL1;
    Int_t pulseIndexL2;
    Int_t pmtL1;
    Int_t pmtL2;
    Double_t timeL1;
    Double_t timeL2;
    Double_t meanTime;
    Double_t deltaTime;
    Double_t deltaX;
    Double_t deltaY;
    Double_t score;
    Double_t ecalResidual;
    Double_t trajectoryResidual;
    Double_t ecalScore;
    Int_t yTopology;
    Bool_t greedySelected;
    Int_t selectedPairIndex;
  };

  /// Read-only snapshot of one accepted candidate in an event containing
  /// fully selected pulses in exactly one CDet layer.
  struct SingleLayerCandidate {
    Int_t index;
    Int_t pulseIndex;
    Int_t pmt;
    Int_t layer;
    Double_t correctedLeadingEdge;
    Double_t ecalResidual;
    Double_t xResidual;
    Double_t score;
  };

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
  Int_t                 GetNumLayerPairCandidates() const
    { return static_cast<Int_t>(fPairCandidateIndex.size()); }
  Int_t                 GetNumSingleLayerCandidates() const
    { return static_cast<Int_t>(fSingleLayerCandidateIndex.size()); }
  Int_t                 GetROICandidateStatus() const
    { return fROICandidateStatus; }

  /// Copy one event-local pulse candidate into pulse. Returns false for an
  /// invalid index. This is the supported in-process interface for consumers
  /// such as the global ROI module; callers do not own detector storage.
  Bool_t                GetPulseCandidate(Int_t index,
                                          PulseCandidate& pulse) const;

  /// Copy one event-local accepted pair into pair. Returns false for an
  /// invalid index. These are the existing greedy-selected pairs; all valid
  /// pre-greedy hypotheses are available through GetLayerPairCandidate().
  Bool_t                GetLayerPair(Int_t index, LayerPair& pair) const;

  /// Copy one event-local pair hypothesis that passed all detector-local hard
  /// gates and the configured ECal ellipse, before greedy pulse sharing is
  /// resolved. Returns false for an invalid index. A pulse may therefore occur
  /// in more than one candidate.
  Bool_t                GetLayerPairCandidate(Int_t index,
                                              LayerPair& pair) const;

  /// Copy one accepted exclusive single-layer candidate into candidate.
  /// Returns false for an invalid index. Every candidate in an event belongs
  /// to the same populated layer; events with good pulses in both layers are
  /// excluded from this collection.
  Bool_t                GetSingleLayerCandidate(
      Int_t index, SingleLayerCandidate& candidate) const;

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
  std::vector<Double_t> fPairTrajectoryResidual;
  std::vector<Double_t> fPairECalScore;
  std::vector<Int_t> fPairYTopology;

  std::vector<Int_t>    fPairCandidateIndex;
  std::vector<Int_t>    fPairCandidatePulseIndexL1;
  std::vector<Int_t>    fPairCandidatePulseIndexL2;
  std::vector<Int_t>    fPairCandidatePMTL1;
  std::vector<Int_t>    fPairCandidatePMTL2;
  std::vector<Double_t> fPairCandidateTimeL1;
  std::vector<Double_t> fPairCandidateTimeL2;
  std::vector<Double_t> fPairCandidateTimeMean;
  std::vector<Double_t> fPairCandidateDeltaTime;
  std::vector<Double_t> fPairCandidateDeltaX;
  std::vector<Double_t> fPairCandidateDeltaY;
  std::vector<Double_t> fPairCandidateScore;
  std::vector<Double_t> fPairCandidateECalResidual;
  std::vector<Double_t> fPairCandidateTrajectoryResidual;
  std::vector<Double_t> fPairCandidateECalScore;
  std::vector<Int_t>    fPairCandidateYTopology;
  std::vector<Int_t>    fPairCandidateGreedySelected;
  std::vector<Int_t>    fPairCandidateSelectedPairIndex;

  std::vector<Int_t>    fSingleLayerCandidateIndex;
  std::vector<Int_t>    fSingleLayerCandidatePulseIndex;
  std::vector<Int_t>    fSingleLayerCandidatePMT;
  std::vector<Int_t>    fSingleLayerCandidateLayer;
  std::vector<Double_t> fSingleLayerCandidateTime;
  std::vector<Double_t> fSingleLayerCandidateECalResidual;
  std::vector<Double_t> fSingleLayerCandidateXResidual;
  std::vector<Double_t> fSingleLayerCandidateScore;

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
  Bool_t fPairingECalRankEnabled;
  Double_t fPairingECalTrajectoryCenter;
  Double_t fPairingECalTimingCenter;
  Double_t fPairingECalTrajectoryScale;
  Double_t fPairingECalTimingScale;
  Double_t fPairingECalRadius;
  Bool_t fPairingOppositeSideEnabled;
  Double_t fPairingOppositeDeltaYCenter;
  Double_t fPairingOppositeDeltaYTolerance;
  Double_t fPairingOppositeProjectedYCenter;
  Double_t fPairingOppositeProjectedYMax;
  Bool_t fSingleLayerEnabled;
  Double_t fSingleLayerResidualCenter;
  Double_t fSingleLayerTimingCenter;
  Double_t fSingleLayerResidualScale;
  Double_t fSingleLayerTimingScale;
  Double_t fSingleLayerRadius;
  Int_t fROICandidateStatus;

  void BuildPulseCandidates();
  void BuildLayerPairs();
  void ClearPulseCandidates();
  Int_t GetTimingStatus() const { return fTimingStatus; }
  Int_t GetTimingECalClusterIndex() const { return fTimingECalClusterIndex; }
  Double_t GetTimingECalTime() const { return fTimingECalTime; }

  ClassDef(SBSCDet,14)  // Describes scintillator plane with F1TDC as a detector
};

#endif
