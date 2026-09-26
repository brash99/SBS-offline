//////////////////////////////////////////////////////////////////////////
//
// SBSCDet class implementation
//
//////////////////////////////////////////////////////////////////////////

#include "SBSCDet.h"

#include <algorithm>
#include <cmath>
#include <limits>

ClassImp(SBSCDet);

/*
 * SBSCDet constructor.
 *
 * Use a TDC with trailing edge info, default is no ADC, but available for
 * commissioning only
 */
SBSCDet::SBSCDet( const char* name, const char* description,
    THaApparatus* apparatus ) : SBSGenericDetector(name,description,apparatus),
    fStorePulseCandidates(false), fTimingCalibrationEnabled(false),
    fTimingCalibrationLoaded(false),
    fTimingStatus(0), fTimingECalClusterIndex(-1),
    fTimingECalTime(std::numeric_limits<Double_t>::quiet_NaN()),
    fTimingTDCToNs(0.01), fTimingECalP0(0.0), fTimingECalP1(0.0),
    fTimingECalDelta(0.0), fTimingWalkP1L1(0.0), fTimingWalkP1L2(0.0),
    fTimingWalkToTRefL1(12.0), fTimingWalkToTRefL2(12.0),
    fTimingShift(0.0), fSelectionEnabled(false),
    fSelectionLEMin(0.0), fSelectionLEMax(0.0),
    fSelectionToTMin(0.0), fSelectionToTMax(0.0),
    fSelectionECalTimeMin(0.0), fSelectionECalTimeMax(0.0),
    fSelectionECalXMin(0.0), fSelectionECalXMax(0.0),
    fSelectionECalYMin(0.0), fSelectionECalYMax(0.0),
    fSelectionXScaleL1(1.0), fSelectionXScaleL2(1.0),
    fSelectionXAlignmentL1(0.0), fSelectionXAlignmentL2(0.0),
    fSelectionXResidualOffset(0.0), fSelectionXResidualMax(0.0),
    fSelectionYResidualOffset(0.0), fSelectionYResidualMax(0.0),
    fPairingEnabled(false), fPairingDeltaTimeCenter(0.0),
    fPairingDeltaTimeMax(0.0), fPairingDeltaXMax(0.0),
    fPairingDeltaYMax(0.0), fPairingTimeScale(1.0),
    fPairingXScale(0.01), fPairingAllowMultiple(true),
    fPairingECalRankEnabled(false), fPairingECalTrajectoryCenter(0.0),
    fPairingECalTimingCenter(0.0), fPairingECalTrajectoryScale(0.02),
    fPairingECalTimingScale(5.0), fPairingECalRadius(2.0),
    fPairingOppositeSideEnabled(false),
    fPairingOppositeDeltaYCenter(0.51),
    fPairingOppositeDeltaYTolerance(0.08),
    fPairingOppositeProjectedYCenter(0.0),
    fPairingOppositeProjectedYMax(0.17),
    fSingleLayerEnabled(false), fSingleLayerResidualCenter(0.0),
    fSingleLayerTimingCenter(-26.0), fSingleLayerResidualScale(0.040),
    fSingleLayerTimingScale(5.0), fSingleLayerRadius(2.0),
    fROICandidateStatus(0)
{
  SetModeTDC(SBSModeTDC::kTDC); //  A TDC with leading & trailing edge info
  SetModeADC(SBSModeADC::kNone); // Default is No ADC, but can be re-enabled later

  fHits             = new TClonesArray("SBSCDet_Hit",200);
  fHit_tmin 	    = 100;
  fHit_tmax	    = 10000; // 1ns to 100ns  units of TDC are 10's of picoseconds?
  fHit_totmin 	    = 100;
  fHit_totmax	    = 10000; // 1ns to 100ns  units of TDC are 10's of picoseconds?

  Clear();
}

SBSCDet::~SBSCDet()
{
  // Destructor. Remove variables from global list and free up the memory
  // allocated by us.
  Clear();// so the prgram doesn't complain when deleting clusters
  RemoveVariables();
  delete fHits;
}


///////////////////////////////////////////////////////////////////////////////
/// Read SBSCDet Database
Int_t SBSCDet::ReadDatabase( const TDatime& date )
{
  // We can use this name here for logs
  //static const char* const here = "ReadDatabase()";

  // If we want to add any new variables, uncomment the following and add
  // the new variables we want to read from the database
  FILE* fi = OpenFile( date );
  if( !fi ) return kFileError;
  //Int_t err;

  std::cout<<"SBSCDet::ReadDatabase method"<<std::endl;
  
  Int_t err = SBSGenericDetector::ReadDatabase(date);
  if(err) {
    fclose(fi);
    return err;
  }
  fIsInit = false;

  std::vector<Double_t> xpos,ypos,zpos;
  Int_t timingEnabled = 0;
  Int_t selectionEnabled = 0;
  Int_t pairingEnabled = 0;
  Int_t pairingAllowMultiple = 1;
  Int_t pairingECalRankEnabled = 0;
  Int_t pairingOppositeSideEnabled = 0;
  Int_t singleLayerEnabled = -1;
  std::vector<Double_t> timingPixelOffset;

  DBRequest config_request[] = {
    { "xpos", &xpos,    kDoubleV, 0, 1 },
    { "ypos", &ypos,    kDoubleV, 0, 1 },
    { "zpos", &zpos,    kDoubleV, 0, 1 },
    { "timing.enable", &timingEnabled, kInt, 0, 1 },
    { "timing.tdc_to_ns", &fTimingTDCToNs, kDouble, 0, 1 },
    { "timing.pixel_offset", &timingPixelOffset, kDoubleV, 0, 1 },
    { "timing.ecal_p0", &fTimingECalP0, kDouble, 0, 1 },
    { "timing.ecal_p1", &fTimingECalP1, kDouble, 0, 1 },
    { "timing.ecal_delta", &fTimingECalDelta, kDouble, 0, 1 },
    { "timing.timewalk_p1_l1", &fTimingWalkP1L1, kDouble, 0, 1 },
    { "timing.timewalk_p1_l2", &fTimingWalkP1L2, kDouble, 0, 1 },
    { "timing.timewalk_totref_l1", &fTimingWalkToTRefL1, kDouble, 0, 1 },
    { "timing.timewalk_totref_l2", &fTimingWalkToTRefL2, kDouble, 0, 1 },
    { "timing.shift_ns", &fTimingShift, kDouble, 0, 1 },
    { "selection.enable", &selectionEnabled, kInt, 0, 1 },
    { "selection.le_min", &fSelectionLEMin, kDouble, 0, 1 },
    { "selection.le_max", &fSelectionLEMax, kDouble, 0, 1 },
    { "selection.tot_min", &fSelectionToTMin, kDouble, 0, 1 },
    { "selection.tot_max", &fSelectionToTMax, kDouble, 0, 1 },
    { "selection.ecal_time_min", &fSelectionECalTimeMin, kDouble, 0, 1 },
    { "selection.ecal_time_max", &fSelectionECalTimeMax, kDouble, 0, 1 },
    { "selection.ecal_x_min", &fSelectionECalXMin, kDouble, 0, 1 },
    { "selection.ecal_x_max", &fSelectionECalXMax, kDouble, 0, 1 },
    { "selection.ecal_y_min", &fSelectionECalYMin, kDouble, 0, 1 },
    { "selection.ecal_y_max", &fSelectionECalYMax, kDouble, 0, 1 },
    { "selection.x_scale_l1", &fSelectionXScaleL1, kDouble, 0, 1 },
    { "selection.x_scale_l2", &fSelectionXScaleL2, kDouble, 0, 1 },
    { "selection.x_alignment_l1", &fSelectionXAlignmentL1, kDouble, 0, 1 },
    { "selection.x_alignment_l2", &fSelectionXAlignmentL2, kDouble, 0, 1 },
    { "selection.x_residual_offset", &fSelectionXResidualOffset, kDouble, 0, 1 },
    { "selection.x_residual_max", &fSelectionXResidualMax, kDouble, 0, 1 },
    { "selection.y_residual_offset", &fSelectionYResidualOffset, kDouble, 0, 1 },
    { "selection.y_residual_max", &fSelectionYResidualMax, kDouble, 0, 1 },
    { "pairing.enable", &pairingEnabled, kInt, 0, 1 },
    { "pairing.dt_center", &fPairingDeltaTimeCenter, kDouble, 0, 1 },
    { "pairing.dt_max", &fPairingDeltaTimeMax, kDouble, 0, 1 },
    { "pairing.dx_max", &fPairingDeltaXMax, kDouble, 0, 1 },
    { "pairing.dy_max", &fPairingDeltaYMax, kDouble, 0, 1 },
    { "pairing.time_scale", &fPairingTimeScale, kDouble, 0, 1 },
    { "pairing.x_scale", &fPairingXScale, kDouble, 0, 1 },
    { "pairing.allow_multiple", &pairingAllowMultiple, kInt, 0, 1 },
    { "pairing.ecal_rank_enable", &pairingECalRankEnabled, kInt, 0, 1 },
    { "pairing.ecal_trajectory_center", &fPairingECalTrajectoryCenter, kDouble, 0, 1 },
    { "pairing.ecal_timing_center", &fPairingECalTimingCenter, kDouble, 0, 1 },
    { "pairing.ecal_trajectory_scale", &fPairingECalTrajectoryScale, kDouble, 0, 1 },
    { "pairing.ecal_timing_scale", &fPairingECalTimingScale, kDouble, 0, 1 },
    { "pairing.ecal_radius", &fPairingECalRadius, kDouble, 0, 1 },
    { "pairing.opposite_side_enable", &pairingOppositeSideEnabled, kInt, 0, 1 },
    { "pairing.opposite_dy_center", &fPairingOppositeDeltaYCenter, kDouble, 0, 1 },
    { "pairing.opposite_dy_tolerance", &fPairingOppositeDeltaYTolerance, kDouble, 0, 1 },
    { "pairing.opposite_projected_y_center", &fPairingOppositeProjectedYCenter, kDouble, 0, 1 },
    { "pairing.opposite_projected_y_max", &fPairingOppositeProjectedYMax, kDouble, 0, 1 },
    { "single_layer.enable", &singleLayerEnabled, kInt, 0, 1 },
    { "single_layer.residual_center", &fSingleLayerResidualCenter, kDouble, 0, 1 },
    { "single_layer.timing_center", &fSingleLayerTimingCenter, kDouble, 0, 1 },
    { "single_layer.residual_scale", &fSingleLayerResidualScale, kDouble, 0, 1 },
    { "single_layer.timing_scale", &fSingleLayerTimingScale, kDouble, 0, 1 },
    { "single_layer.radius", &fSingleLayerRadius, kDouble, 0, 1 },
    { 0 } ///< Request must end in a NULL
  };
  err = LoadDB( fi, date, config_request, fPrefix );
  if (err) {
    fclose(fi);
    return err;
  }

  const auto invalidGeometry = [this](const std::vector<Double_t>& values) {
    return !values.empty() && static_cast<Int_t>(values.size()) != fNelem;
  };
  if (invalidGeometry(xpos) || invalidGeometry(ypos) || invalidGeometry(zpos)) {
    Error(Here("ReadDatabase"),
          "CDet geometry vector length does not match number of elements (%d)",
          fNelem);
    fclose(fi);
    return kInitError;
  }

  if (!xpos.empty()) {
    if ((int)xpos.size() == fNelem) {
      for (Int_t ne=0;ne<fNelem;ne++) {
        fElements[ne]->SetX(xpos[ne]);
	//std::cout << "ne = " << ne << " xpos = " << xpos[ne] << std::endl;
      }
    }
  }

  if (!ypos.empty()) {
    if ((int)ypos.size() == fNelem) {
      for (Int_t ne=0;ne<fNelem;ne++) {
        fElements[ne]->SetY(ypos[ne]);
	//std::cout << "ne = " << ne << " ypos = " << ypos[ne] << std::endl;
      }
    }
  }
  
  if (!zpos.empty()) {
    if ((int)zpos.size() == fNelem) {
      for (Int_t ne=0;ne<fNelem;ne++) {
        fElements[ne]->SetZ(zpos[ne]);
	//std::cout << "ne = " << ne << " zpos = " << zpos[ne] << std::endl;
      }
    }
  }

  fTimingCalibrationEnabled = timingEnabled != 0;
  fTimingCalibrationLoaded = false;
  if (fTimingCalibrationEnabled) {
    const Int_t nPhysicalPixels = std::min<Int_t>(2688, fNelem);
    if (static_cast<Int_t>(timingPixelOffset.size()) != nPhysicalPixels ||
        fTimingTDCToNs <= 0.0 || fTimingWalkToTRefL1 <= 0.0 ||
        fTimingWalkToTRefL2 <= 0.0) {
      Error(Here("ReadDatabase"),
            "Invalid CDet timing calibration: expected %d pixel offsets, got %zu",
            nPhysicalPixels, timingPixelOffset.size());
      fclose(fi);
      return kInitError;
    }
    fTimingPixelOffset = timingPixelOffset;
    fTimingCalibrationLoaded = true;
  }

  fSelectionEnabled = selectionEnabled != 0;
  if (fSelectionEnabled &&
      (fSelectionLEMin >= fSelectionLEMax ||
       fSelectionToTMin >= fSelectionToTMax ||
       fSelectionECalTimeMin >= fSelectionECalTimeMax ||
       fSelectionECalXMin >= fSelectionECalXMax ||
       fSelectionECalYMin >= fSelectionECalYMax ||
       fSelectionXResidualMax <= 0.0 || fSelectionYResidualMax <= 0.0)) {
    Error(Here("ReadDatabase"), "Invalid CDet selection database ranges");
    fclose(fi);
    return kInitError;
  }

  fPairingEnabled = pairingEnabled != 0;
  fPairingAllowMultiple = pairingAllowMultiple != 0;
  fPairingECalRankEnabled = pairingECalRankEnabled != 0;
  fPairingOppositeSideEnabled = pairingOppositeSideEnabled != 0;
  fSingleLayerEnabled = singleLayerEnabled < 0 ? fPairingECalRankEnabled
                                               : singleLayerEnabled != 0;
  if (fPairingEnabled &&
      (fPairingDeltaTimeMax <= 0.0 || fPairingDeltaXMax <= 0.0 ||
       fPairingDeltaYMax <= 0.0 || fPairingTimeScale <= 0.0 ||
       fPairingXScale <= 0.0)) {
    Error(Here("ReadDatabase"), "Invalid CDet pairing database parameters");
    fclose(fi);
    return kInitError;
  }
  if (fPairingECalRankEnabled &&
      (fPairingECalTrajectoryScale <= 0.0 ||
       fPairingECalTimingScale <= 0.0 || fPairingECalRadius <= 0.0)) {
    Error(Here("ReadDatabase"),
          "Invalid CDet ECal-informed pairing parameters");
    fclose(fi);
    return kInitError;
  }
  if (fPairingOppositeSideEnabled &&
      (fPairingOppositeDeltaYCenter <= 0.0 ||
       fPairingOppositeDeltaYTolerance <= 0.0 ||
       fPairingOppositeProjectedYMax <= 0.0)) {
    Error(Here("ReadDatabase"),
          "Invalid CDet opposite-side pairing parameters");
    fclose(fi);
    return kInitError;
  }
  if (fSingleLayerEnabled &&
      (fSingleLayerResidualScale <= 0.0 ||
       fSingleLayerTimingScale <= 0.0 || fSingleLayerRadius <= 0.0)) {
    Error(Here("ReadDatabase"),
          "Invalid CDet exclusive single-layer selection parameters");
    fclose(fi);
    return kInitError;
  }

  fIsInit = true;

  fclose(fi);
  return kOK;

  // Make sure to call parent class so that the generic variables can be read
  //return SBSGenericDetector::ReadDatabase(date);

}

//_____________________________________________________________________________
Int_t SBSCDet::DefineVariables( EMode mode )
{
  // Initialize global variables
  Int_t err = SBSGenericDetector::DefineVariables(mode);
  if(err) {
    return err;
  }

  // Uncomment the following to add a ny new variables we want to define
  // as the output

  RVarDef vars[] = {
   { "ngoodhits",       " number of Good PMT hits", "GetNumHits()"  },
   { "hit.pmtnum",  " Hit PMT num",        "fHits.SBSCDet_Hit.GetPMTNum()"},
   { "hit.row",     " PMT hit row",        "fHits.SBSCDet_Hit.GetRow()"   },
   { "hit.col",     " PMT hit column",     "fHits.SBSCDet_Hit.GetCol()"   },
   { "hit.layer",     " PMT hit layer",     "fHits.SBSCDet_Hit.GetLayer()"   },
   { "hit.xhit",    " PMT hit X",          "fHits.SBSCDet_Hit.GetX()"     },
   { "hit.yhit",    " PMT hit Y",          "fHits.SBSCDet_Hit.GetY()"     },
   { "hit.zhit",    " PMT hit Z",          "fHits.SBSCDet_Hit.GetZ()"     },
   { "hit.tdc_le",   " PMT hit TDC LE",  "fHits.SBSCDet_Hit.GetTDC_LE()" },
   { "hit.tdc_te",   " PMT hit TDC TE",   "fHits.SBSCDet_Hit.GetTDC_TE()" },
   { "hit.tdc_tot",   " PMT hit TDC TOT",   "fHits.SBSCDet_Hit.GetToT()" },
   { "pulse.n",       " Number of decoded CDet TDC pulse candidates", "GetNumPulseCandidates()" },
   { "pulse.n_le_only", " Number of decoder slots containing only a leading edge", "GetNumLeadingOnlySlots()" },
   { "pulse.n_te_only", " Number of decoder slots containing only a trailing edge", "GetNumTrailingOnlySlots()" },
   { "pulse.n_invalid_pair", " Number of two-edge slots with nonpositive ToT", "GetNumInvalidPairSlots()" },
   { "pulse.pmtnum",  " Pulse candidate PMT/channel number", "fPulsePMT" },
   { "pulse.index",   " Accepted pulse index within its channel", "fPulseIndex" },
   { "pulse.le_index", " Original decoder-slot index of the accepted LE", "fPulseLEIndex" },
   { "pulse.te_index", " Original decoder-slot index of the accepted TE", "fPulseTEIndex" },
   { "pulse.row",     " Pulse candidate row", "fPulseRow" },
   { "pulse.col",     " Pulse candidate column", "fPulseCol" },
   { "pulse.layer",   " Pulse candidate physical layer (0=Layer 1, 1=Layer 2)", "fPulseLayer" },
   { "pulse.x",       " Pulse candidate channel X position", "fPulseX" },
   { "pulse.y",       " Pulse candidate channel Y position", "fPulseY" },
   { "pulse.z",       " Pulse candidate channel Z position", "fPulseZ" },
   { "pulse.tdc_le",  " Pulse candidate calibrated leading-edge time", "fPulseLE" },
   { "pulse.tdc_te",  " Pulse candidate calibrated trailing-edge time", "fPulseTE" },
   { "pulse.tdc_tot", " Pulse candidate calibrated time over threshold", "fPulseToT" },
   { "pulse.tdc_le_raw",  " Pulse candidate raw leading-edge value", "fPulseLERaw" },
   { "pulse.tdc_te_raw",  " Pulse candidate raw trailing-edge value", "fPulseTERaw" },
   { "pulse.tdc_tot_raw", " Pulse candidate raw time-over-threshold value", "fPulseToTRaw" },
   { "pulse.tdc_le_corr", " Fully timing-corrected pulse leading edge in ns", "fPulseLECorrected" },
   { "pulse.tdc_te_corr", " Fully timing-corrected pulse trailing edge in ns", "fPulseTECorrected" },
   { "pulse.tdc_tot_ns", " Pulse time over threshold in ns", "fPulseToTNs" },
   { "pulse.ecal_residual", " ECal time minus corrected CDet leading edge in ns", "fPulseECalResidual" },
   { "pulse.calib_valid", " Pulse timing-calibration validity flag", "fPulseCalibrationValid" },
   { "pulse.x_corr", " Layer-aligned CDet x coordinate in m", "fPulseCorrectedX" },
   { "pulse.ecal_x_proj", " ECal x projected to the CDet hit z in m", "fPulseProjectedECalX" },
   { "pulse.ecal_y_proj", " ECal y projected to the CDet hit z in m", "fPulseProjectedECalY" },
   { "pulse.ecal_x_residual", " Aligned CDet x minus projected ECal x and configured offset in m", "fPulseECalXResidual" },
   { "pulse.ecal_y_residual", " CDet y minus projected ECal y and configured offset in m", "fPulseECalYResidual" },
   { "pulse.broad_quality_pass", " Broad LE and ToT quality decision", "fPulseBroadQualityPass" },
   { "pulse.ecal_eligible", " ECal cluster position and time eligibility decision", "fPulseECalEligibilityPass" },
   { "pulse.spatial_pass", " ECal-to-CDet projected x and y compatibility decision", "fPulseSpatialPass" },
   { "pair.n", " Number of accepted one-to-one CDet layer pairs", "GetNumLayerPairs()" },
   { "pair.index", " Pair index within the event", "fPairIndex" },
   { "pair.pulse_index_l1", " Index in pulse arrays for Layer-1 member", "fPairPulseIndexL1" },
   { "pair.pulse_index_l2", " Index in pulse arrays for Layer-2 member", "fPairPulseIndexL2" },
   { "pair.pmtnum_l1", " Layer-1 member pixel ID", "fPairPMTL1" },
   { "pair.pmtnum_l2", " Layer-2 member pixel ID", "fPairPMTL2" },
   { "pair.time_l1", " Corrected Layer-1 leading-edge time in ns", "fPairTimeL1" },
   { "pair.time_l2", " Corrected Layer-2 leading-edge time in ns", "fPairTimeL2" },
   { "pair.time_mean", " Mean corrected pair time in ns", "fPairTimeMean" },
   { "pair.dt", " Layer-2 minus Layer-1 corrected time in ns", "fPairDeltaTime" },
   { "pair.dx", " Layer-2 minus Layer-1 aligned x in m", "fPairDeltaX" },
   { "pair.dy", " Layer-2 minus Layer-1 y in m", "fPairDeltaY" },
   { "pair.score", " CDet-only layer-pair ranking score", "fPairScore" },
   { "pair.ecal_residual", " ECal time minus corrected pair mean time in ns", "fPairECalResidual" },
   { "pair.trajectory_residual", " Inter-layer x residual relative to the ECal trajectory in m", "fPairTrajectoryResidual" },
   { "pair.ecal_score", " ECal-informed normalized trajectory-time score", "fPairECalScore" },
   { "pair.y_topology", " Pair y topology: 0 same-side, 1 opposite-side seam", "fPairYTopology" },
   { "pair_candidate.n", " Number of valid CDet layer-pair hypotheses before greedy assignment", "GetNumLayerPairCandidates()" },
   { "pair_candidate.index", " Pair-candidate rank within the event", "fPairCandidateIndex" },
   { "pair_candidate.pulse_index_l1", " Index in pulse arrays for Layer-1 member", "fPairCandidatePulseIndexL1" },
   { "pair_candidate.pulse_index_l2", " Index in pulse arrays for Layer-2 member", "fPairCandidatePulseIndexL2" },
   { "pair_candidate.pmtnum_l1", " Layer-1 member pixel ID", "fPairCandidatePMTL1" },
   { "pair_candidate.pmtnum_l2", " Layer-2 member pixel ID", "fPairCandidatePMTL2" },
   { "pair_candidate.time_l1", " Corrected Layer-1 leading-edge time in ns", "fPairCandidateTimeL1" },
   { "pair_candidate.time_l2", " Corrected Layer-2 leading-edge time in ns", "fPairCandidateTimeL2" },
   { "pair_candidate.time_mean", " Mean corrected candidate time in ns", "fPairCandidateTimeMean" },
   { "pair_candidate.dt", " Layer-2 minus Layer-1 corrected time in ns", "fPairCandidateDeltaTime" },
   { "pair_candidate.dx", " Layer-2 minus Layer-1 aligned x in m", "fPairCandidateDeltaX" },
   { "pair_candidate.dy", " Layer-2 minus Layer-1 y in m", "fPairCandidateDeltaY" },
   { "pair_candidate.score", " CDet-only layer-pair ranking score", "fPairCandidateScore" },
   { "pair_candidate.ecal_residual", " ECal time minus corrected candidate mean time in ns", "fPairCandidateECalResidual" },
   { "pair_candidate.trajectory_residual", " Inter-layer x residual relative to the ECal trajectory in m", "fPairCandidateTrajectoryResidual" },
   { "pair_candidate.ecal_score", " ECal-informed normalized trajectory-time score", "fPairCandidateECalScore" },
   { "pair_candidate.y_topology", " Candidate y topology: 0 same-side, 1 opposite-side seam", "fPairCandidateYTopology" },
   { "pair_candidate.greedy_selected", " One if this hypothesis is retained in pair.*", "fPairCandidateGreedySelected" },
   { "pair_candidate.selected_pair_index", " Corresponding pair.* index, or -1 if not selected", "fPairCandidateSelectedPairIndex" },
   { "single_candidate.n", " Number of accepted exclusive single-layer candidates", "GetNumSingleLayerCandidates()" },
   { "single_candidate.index", " Single-layer candidate index within the event", "fSingleLayerCandidateIndex" },
   { "single_candidate.pulse_index", " Index in pulse arrays for the source pulse", "fSingleLayerCandidatePulseIndex" },
   { "single_candidate.pmtnum", " Single-layer candidate pixel ID", "fSingleLayerCandidatePMT" },
   { "single_candidate.layer", " Single-layer candidate physical layer (0=Layer 1, 1=Layer 2)", "fSingleLayerCandidateLayer" },
   { "single_candidate.time", " Corrected leading-edge time in ns", "fSingleLayerCandidateTime" },
   { "single_candidate.ecal_residual", " ECal time minus corrected leading-edge time in ns", "fSingleLayerCandidateECalResidual" },
   { "single_candidate.x_residual", " Aligned CDet x minus projected ECal x in m", "fSingleLayerCandidateXResidual" },
   { "single_candidate.score", " Normalized single-layer x-time ellipse score", "fSingleLayerCandidateScore" },
   { "roi.status", " ROI candidate status: 0 none, 1 pair hypotheses, 2 Layer-1-only, 3 Layer-2-only", "GetROICandidateStatus()" },
   { "timing.status", " CDet timing status: 0 disabled, 1 missing ECal, 2 applied, -1 invalid calibration", "fTimingStatus" },
   { "timing.ecal_cluster", " ECal cluster index used by CDet timing", "fTimingECalClusterIndex" },
   { "timing.ecal_time", " ECal cluster energy-weighted ADC time used by CDet timing", "fTimingECalTime" },
   { 0 }
  };
  err = DefineVarsFromList( vars, mode );
 
 
  // Finally go back
  return err;
}

//_____________________________________________________________________________
Bool_t SBSCDet::GetPulseCandidate(Int_t index,
                                  PulseCandidate& pulse) const
{
  if (index < 0 || static_cast<size_t>(index) >= fPulsePMT.size())
    return false;

  const size_t i = static_cast<size_t>(index);
  pulse.pmt = fPulsePMT[i];
  pulse.index = fPulseIndex[i];
  pulse.leIndex = fPulseLEIndex[i];
  pulse.teIndex = fPulseTEIndex[i];
  pulse.row = fPulseRow[i];
  pulse.col = fPulseCol[i];
  pulse.layer = fPulseLayer[i];
  pulse.x = fPulseX[i];
  pulse.y = fPulseY[i];
  pulse.z = fPulseZ[i];
  pulse.leadingEdge = fPulseLE[i];
  pulse.trailingEdge = fPulseTE[i];
  pulse.timeOverThreshold = fPulseToT[i];
  pulse.rawLeadingEdge = fPulseLERaw[i];
  pulse.rawTrailingEdge = fPulseTERaw[i];
  pulse.rawTimeOverThreshold = fPulseToTRaw[i];
  pulse.correctedLeadingEdge = fPulseLECorrected[i];
  pulse.correctedTrailingEdge = fPulseTECorrected[i];
  pulse.timeOverThresholdNs = fPulseToTNs[i];
  pulse.ecalResidual = fPulseECalResidual[i];
  pulse.calibrationValid = fPulseCalibrationValid[i] != 0;
  pulse.correctedX = fPulseCorrectedX[i];
  pulse.projectedECalX = fPulseProjectedECalX[i];
  pulse.projectedECalY = fPulseProjectedECalY[i];
  pulse.ecalXResidual = fPulseECalXResidual[i];
  pulse.ecalYResidual = fPulseECalYResidual[i];
  pulse.broadQualityPass = fPulseBroadQualityPass[i] != 0;
  pulse.ecalEligibilityPass = fPulseECalEligibilityPass[i] != 0;
  pulse.spatialPass = fPulseSpatialPass[i] != 0;
  return true;
}

//_____________________________________________________________________________
Bool_t SBSCDet::GetLayerPair(Int_t index, LayerPair& pair) const
{
  if (index < 0 || static_cast<size_t>(index) >= fPairIndex.size())
    return false;

  const size_t i = static_cast<size_t>(index);
  pair.index = fPairIndex[i];
  pair.pulseIndexL1 = fPairPulseIndexL1[i];
  pair.pulseIndexL2 = fPairPulseIndexL2[i];
  pair.pmtL1 = fPairPMTL1[i];
  pair.pmtL2 = fPairPMTL2[i];
  pair.timeL1 = fPairTimeL1[i];
  pair.timeL2 = fPairTimeL2[i];
  pair.meanTime = fPairTimeMean[i];
  pair.deltaTime = fPairDeltaTime[i];
  pair.deltaX = fPairDeltaX[i];
  pair.deltaY = fPairDeltaY[i];
  pair.score = fPairScore[i];
  pair.ecalResidual = fPairECalResidual[i];
  pair.trajectoryResidual = fPairTrajectoryResidual[i];
  pair.ecalScore = fPairECalScore[i];
  pair.yTopology = fPairYTopology[i];
  pair.greedySelected = true;
  pair.selectedPairIndex = fPairIndex[i];
  return true;
}

//_____________________________________________________________________________
Bool_t SBSCDet::GetLayerPairCandidate(Int_t index, LayerPair& pair) const
{
  if (index < 0 ||
      static_cast<size_t>(index) >= fPairCandidateIndex.size())
    return false;

  const size_t i = static_cast<size_t>(index);
  pair.index = fPairCandidateIndex[i];
  pair.pulseIndexL1 = fPairCandidatePulseIndexL1[i];
  pair.pulseIndexL2 = fPairCandidatePulseIndexL2[i];
  pair.pmtL1 = fPairCandidatePMTL1[i];
  pair.pmtL2 = fPairCandidatePMTL2[i];
  pair.timeL1 = fPairCandidateTimeL1[i];
  pair.timeL2 = fPairCandidateTimeL2[i];
  pair.meanTime = fPairCandidateTimeMean[i];
  pair.deltaTime = fPairCandidateDeltaTime[i];
  pair.deltaX = fPairCandidateDeltaX[i];
  pair.deltaY = fPairCandidateDeltaY[i];
  pair.score = fPairCandidateScore[i];
  pair.ecalResidual = fPairCandidateECalResidual[i];
  pair.trajectoryResidual = fPairCandidateTrajectoryResidual[i];
  pair.ecalScore = fPairCandidateECalScore[i];
  pair.yTopology = fPairCandidateYTopology[i];
  pair.greedySelected = fPairCandidateGreedySelected[i] != 0;
  pair.selectedPairIndex = fPairCandidateSelectedPairIndex[i];
  return true;
}

//_____________________________________________________________________________
Bool_t SBSCDet::GetSingleLayerCandidate(
    Int_t index, SingleLayerCandidate& candidate) const
{
  if (index < 0 ||
      static_cast<size_t>(index) >= fSingleLayerCandidateIndex.size())
    return false;

  const size_t i = static_cast<size_t>(index);
  candidate.index = fSingleLayerCandidateIndex[i];
  candidate.pulseIndex = fSingleLayerCandidatePulseIndex[i];
  candidate.pmt = fSingleLayerCandidatePMT[i];
  candidate.layer = fSingleLayerCandidateLayer[i];
  candidate.correctedLeadingEdge = fSingleLayerCandidateTime[i];
  candidate.ecalResidual = fSingleLayerCandidateECalResidual[i];
  candidate.xResidual = fSingleLayerCandidateXResidual[i];
  candidate.score = fSingleLayerCandidateScore[i];
  return true;
}

//_____________________________________________________________________________
Int_t SBSCDet::Decode( const THaEvData& evdata )
{
  //std::cout << "SBSCDet::Decode" << std::endl;
  Int_t err = SBSGenericDetector::Decode(evdata);
  return err;
}
//


/*
 * Clear()
 * called at the end of every event
 */
void SBSCDet::Clear( Option_t* opt )
{
  // If we defined any new variables that we need to clear prior to the next event
  // clear them here:
  // fExample = 0.0;

  // Make sure to call parent class's Clear() also!
  SBSGenericDetector::Clear(opt);
  fHits->Clear("C");
  ClearPulseCandidates();
}

void SBSCDet::ClearPulseCandidates()
{
  fNLeadingOnlySlots = 0;
  fNTrailingOnlySlots = 0;
  fNInvalidPairSlots = 0;
  fPulsePMT.clear();
  fPulseIndex.clear();
  fPulseLEIndex.clear();
  fPulseTEIndex.clear();
  fPulseRow.clear();
  fPulseCol.clear();
  fPulseLayer.clear();
  fPulseX.clear();
  fPulseY.clear();
  fPulseZ.clear();
  fPulseLE.clear();
  fPulseTE.clear();
  fPulseToT.clear();
  fPulseLERaw.clear();
  fPulseTERaw.clear();
  fPulseToTRaw.clear();
  fPulseLECorrected.clear();
  fPulseTECorrected.clear();
  fPulseToTNs.clear();
  fPulseECalResidual.clear();
  fPulseCalibrationValid.clear();
  fPulseCorrectedX.clear();
  fPulseProjectedECalX.clear();
  fPulseProjectedECalY.clear();
  fPulseECalXResidual.clear();
  fPulseECalYResidual.clear();
  fPulseBroadQualityPass.clear();
  fPulseECalEligibilityPass.clear();
  fPulseSpatialPass.clear();
  fPairIndex.clear();
  fPairPulseIndexL1.clear();
  fPairPulseIndexL2.clear();
  fPairPMTL1.clear();
  fPairPMTL2.clear();
  fPairTimeL1.clear();
  fPairTimeL2.clear();
  fPairTimeMean.clear();
  fPairDeltaTime.clear();
  fPairDeltaX.clear();
  fPairDeltaY.clear();
  fPairScore.clear();
  fPairECalResidual.clear();
  fPairTrajectoryResidual.clear();
  fPairECalScore.clear();
  fPairYTopology.clear();
  fPairCandidateIndex.clear();
  fPairCandidatePulseIndexL1.clear();
  fPairCandidatePulseIndexL2.clear();
  fPairCandidatePMTL1.clear();
  fPairCandidatePMTL2.clear();
  fPairCandidateTimeL1.clear();
  fPairCandidateTimeL2.clear();
  fPairCandidateTimeMean.clear();
  fPairCandidateDeltaTime.clear();
  fPairCandidateDeltaX.clear();
  fPairCandidateDeltaY.clear();
  fPairCandidateScore.clear();
  fPairCandidateECalResidual.clear();
  fPairCandidateTrajectoryResidual.clear();
  fPairCandidateECalScore.clear();
  fPairCandidateYTopology.clear();
  fPairCandidateGreedySelected.clear();
  fPairCandidateSelectedPairIndex.clear();
  fSingleLayerCandidateIndex.clear();
  fSingleLayerCandidatePulseIndex.clear();
  fSingleLayerCandidatePMT.clear();
  fSingleLayerCandidateLayer.clear();
  fSingleLayerCandidateTime.clear();
  fSingleLayerCandidateECalResidual.clear();
  fSingleLayerCandidateXResidual.clear();
  fSingleLayerCandidateScore.clear();
  fROICandidateStatus = 0;
  fTimingStatus = fTimingCalibrationEnabled ? 1 : 0;
  fTimingECalClusterIndex = -1;
  fTimingECalTime = std::numeric_limits<Double_t>::quiet_NaN();
}

void SBSCDet::BuildPulseCandidates()
{
  ClearPulseCandidates();
  if (!fStorePulseCandidates)
    return;

  const Int_t nPhysicalPixels = std::min<Int_t>(2688, fNelem);
  for (Int_t i = 0; i < nPhysicalPixels; ++i) {
    SBSElement* element = fElements[i];
    if (!element || !element->TDC() || !element->TDC()->HasData())
      continue;

    struct Edge {
      Double_t raw;
      Double_t val;
      size_t slot;
      Bool_t leading;
    };

    const std::vector<SBSData::TDCHit> decoded = element->TDC()->GetAllHits();
    std::vector<Edge> edges;
    edges.reserve(2 * decoded.size());
    for (size_t islot = 0; islot < decoded.size(); ++islot) {
      if (decoded[islot].le.raw != 0.0)
        edges.push_back({decoded[islot].le.raw, decoded[islot].le.val,
                         islot, true});
      if (decoded[islot].te.raw != 0.0)
        edges.push_back({decoded[islot].te.raw, decoded[islot].te.val,
                         islot, false});
    }
    std::stable_sort(edges.begin(), edges.end(), [](const Edge& left,
                                                     const Edge& right) {
      if (left.raw != right.raw)
        return left.raw < right.raw;
      return left.leading && !right.leading;
    });

    Bool_t pulseOpen = false;
    Edge leadingEdge{};
    Int_t acceptedIndex = 0;
    for (const Edge& edge : edges) {
      if (edge.leading) {
        if (!pulseOpen) {
          leadingEdge = edge;
          pulseOpen = true;
        } else {
          // CDet policy for LE-LE-TE-TE: retain the first LE and ignore
          // additional leading edges until the first trailing edge closes it.
          ++fNLeadingOnlySlots;
        }
        continue;
      }
      if (!pulseOpen) {
        ++fNTrailingOnlySlots;
        continue;
      }

      const Double_t rawToT = edge.raw - leadingEdge.raw;
      const Double_t calibratedToT = edge.val - leadingEdge.val;
      pulseOpen = false;
      if (rawToT <= 0.0 || calibratedToT <= 0.0) {
        ++fNInvalidPairSlots;
        continue;
      }

      fPulsePMT.push_back(element->GetID());
      fPulseIndex.push_back(acceptedIndex++);
      fPulseLEIndex.push_back(static_cast<Int_t>(leadingEdge.slot));
      fPulseTEIndex.push_back(static_cast<Int_t>(edge.slot));
      fPulseRow.push_back(element->GetRow());
      fPulseCol.push_back(element->GetCol());
      // The CDet geometry elements do not currently carry a reliable physical
      // layer number. Pixel numbering is authoritative: 0--1343 is Layer 1
      // and 1344--2687 is Layer 2.
      fPulseLayer.push_back(element->GetID() / 1344);
      fPulseX.push_back(element->GetX());
      fPulseY.push_back(element->GetY());
      fPulseZ.push_back(element->GetZ());
      fPulseLE.push_back(leadingEdge.val);
      fPulseTE.push_back(edge.val);
      fPulseToT.push_back(calibratedToT);
      fPulseLERaw.push_back(leadingEdge.raw);
      fPulseTERaw.push_back(edge.raw);
      fPulseToTRaw.push_back(rawToT);
    }
    if (pulseOpen)
      ++fNLeadingOnlySlots;
  }
}

void SBSCDet::ApplyECalTimingCalibration(Double_t ecalTime,
                                         Int_t ecalClusterIndex,
                                         Double_t ecalX, Double_t ecalY,
                                         Double_t ecalZ)
{
  const Double_t missing = std::numeric_limits<Double_t>::quiet_NaN();
  const size_t nPulses = fPulsePMT.size();
  fPulseLECorrected.assign(nPulses, missing);
  fPulseTECorrected.assign(nPulses, missing);
  fPulseToTNs.assign(nPulses, missing);
  fPulseECalResidual.assign(nPulses, missing);
  fPulseCalibrationValid.assign(nPulses, 0);
  fPulseCorrectedX.assign(nPulses, missing);
  fPulseProjectedECalX.assign(nPulses, missing);
  fPulseProjectedECalY.assign(nPulses, missing);
  fPulseECalXResidual.assign(nPulses, missing);
  fPulseECalYResidual.assign(nPulses, missing);
  fPulseBroadQualityPass.assign(nPulses, 0);
  fPulseECalEligibilityPass.assign(nPulses, 0);
  fPulseSpatialPass.assign(nPulses, 0);
  fTimingECalTime = ecalTime;
  fTimingECalClusterIndex = ecalClusterIndex;

  const Bool_t ecalEligible = fSelectionEnabled &&
      std::isfinite(ecalTime) && std::isfinite(ecalX) &&
      std::isfinite(ecalY) && std::isfinite(ecalZ) && ecalZ != 0.0 &&
      ecalClusterIndex >= 0 &&
      ecalTime > fSelectionECalTimeMin && ecalTime < fSelectionECalTimeMax &&
      ecalX > fSelectionECalXMin && ecalX < fSelectionECalXMax &&
      ecalY > fSelectionECalYMin && ecalY < fSelectionECalYMax &&
      ecalX != 0.0 && ecalY != 0.0;

  // Broad pulse quality is detector-local and remains meaningful even when
  // the event has no usable ECal cluster. Keep it independent of the ECal
  // eligibility and projected-spatial decisions.
  if (fSelectionEnabled) {
    for (size_t i = 0; i < nPulses; ++i) {
      const Int_t layer = fPulsePMT[i] / 1344;
      const Double_t leNs = fPulseLE[i] * fTimingTDCToNs;
      const Double_t totNs = fPulseToT[i] * fTimingTDCToNs;
      fPulseBroadQualityPass[i] =
          leNs >= fSelectionLEMin && leNs <= fSelectionLEMax &&
          totNs >= fSelectionToTMin && totNs <= fSelectionToTMax;
      fPulseECalEligibilityPass[i] = ecalEligible;
      const Double_t xScale = layer == 0 ? fSelectionXScaleL1 : fSelectionXScaleL2;
      const Double_t xAlignment =
          layer == 0 ? fSelectionXAlignmentL1 : fSelectionXAlignmentL2;
      fPulseCorrectedX[i] = fPulseX[i] * xScale - xAlignment;
    }
  }

  if (!fTimingCalibrationLoaded) {
    fTimingStatus = fTimingCalibrationEnabled ? -1 : 0;
    return;
  }
  if (!std::isfinite(ecalTime) || ecalClusterIndex < 0) {
    fTimingStatus = 1;
    return;
  }

  for (size_t i = 0; i < nPulses; ++i) {
    const Int_t pixel = fPulsePMT[i];
    if (pixel < 0 || pixel >= 2688 ||
        pixel >= static_cast<Int_t>(fTimingPixelOffset.size()))
      continue;

    const Double_t totNs = fPulseToT[i] * fTimingTDCToNs;
    if (!std::isfinite(totNs) || totNs <= 0.0)
      continue;

    const Int_t layer = pixel / 1344;
    const Double_t walkP1 = layer == 0 ? fTimingWalkP1L1 : fTimingWalkP1L2;
    const Double_t walkReference =
        layer == 0 ? fTimingWalkToTRefL1 : fTimingWalkToTRefL2;
    const Double_t timeWalk = walkP1 *
        (1.0 / std::sqrt(totNs) - 1.0 / std::sqrt(walkReference));
    const Double_t commonCorrection = fTimingPixelOffset[pixel] -
        (fTimingECalP0 + fTimingECalP1 * ecalTime) +
        fTimingECalDelta - timeWalk + fTimingShift;

    fPulseLECorrected[i] = fPulseLE[i] * fTimingTDCToNs + commonCorrection;
    fPulseTECorrected[i] = fPulseTE[i] * fTimingTDCToNs + commonCorrection;
    fPulseToTNs[i] = totNs;
    fPulseECalResidual[i] = ecalTime - fPulseLECorrected[i];
    fPulseCalibrationValid[i] = 1;

    if (fSelectionEnabled && ecalEligible) {
        const Double_t correctedX = fPulseCorrectedX[i];
        const Double_t projectedX = ecalX * fPulseZ[i] / ecalZ;
        const Double_t projectedY = ecalY * fPulseZ[i] / ecalZ;
        const Double_t xResidual = correctedX - projectedX - fSelectionXResidualOffset;
        const Double_t yResidual = fPulseY[i] - projectedY - fSelectionYResidualOffset;
        fPulseCorrectedX[i] = correctedX;
        fPulseProjectedECalX[i] = projectedX;
        fPulseProjectedECalY[i] = projectedY;
        fPulseECalXResidual[i] = xResidual;
        fPulseECalYResidual[i] = yResidual;
        fPulseSpatialPass[i] = std::fabs(xResidual) <= fSelectionXResidualMax &&
            std::fabs(yResidual) <= fSelectionYResidualMax;
    }
  }
  fTimingStatus = 2;
  BuildLayerPairs();
}

void SBSCDet::BuildLayerPairs()
{
  if (!fPairingEnabled || fTimingStatus != 2)
    return;

  std::vector<Int_t> layer1;
  std::vector<Int_t> layer2;
  for (size_t i = 0; i < fPulsePMT.size(); ++i) {
    if (!fPulseCalibrationValid[i] || !fPulseBroadQualityPass[i] ||
        !fPulseECalEligibilityPass[i] || !fPulseSpatialPass[i])
      continue;
    if (fPulsePMT[i] < 1344)
      layer1.push_back(static_cast<Int_t>(i));
    else
      layer2.push_back(static_cast<Int_t>(i));
  }

  // Preserve the hydrogen-study definition of an exclusive single-layer
  // event: fully selected pulses populate exactly one layer. Such an event
  // cannot form a two-layer pair, so every pulse passing the separately
  // validated x-time ellipse remains available to the downstream ROI logic.
  if (fSingleLayerEnabled && (layer1.empty() != layer2.empty())) {
    const std::vector<Int_t>& populatedLayer = layer1.empty() ? layer2 : layer1;
    for (const Int_t pulseIndex : populatedLayer) {
      const Double_t xResidual = fPulseCorrectedX[pulseIndex] -
          fPulseProjectedECalX[pulseIndex];
      const Double_t timingResidual = fPulseECalResidual[pulseIndex];
      if (!std::isfinite(xResidual) || !std::isfinite(timingResidual))
        continue;
      const Double_t xPull =
          (xResidual - fSingleLayerResidualCenter) /
          fSingleLayerResidualScale;
      const Double_t timingPull =
          (timingResidual - fSingleLayerTimingCenter) /
          fSingleLayerTimingScale;
      const Double_t score = xPull*xPull + timingPull*timingPull;
      if (score > fSingleLayerRadius*fSingleLayerRadius)
        continue;
      fSingleLayerCandidateIndex.push_back(
          static_cast<Int_t>(fSingleLayerCandidateIndex.size()));
      fSingleLayerCandidatePulseIndex.push_back(pulseIndex);
      fSingleLayerCandidatePMT.push_back(fPulsePMT[pulseIndex]);
      fSingleLayerCandidateLayer.push_back(fPulseLayer[pulseIndex]);
      fSingleLayerCandidateTime.push_back(fPulseLECorrected[pulseIndex]);
      fSingleLayerCandidateECalResidual.push_back(timingResidual);
      fSingleLayerCandidateXResidual.push_back(xResidual);
      fSingleLayerCandidateScore.push_back(score);
    }
    if (!fSingleLayerCandidateIndex.empty())
      fROICandidateStatus = layer1.empty() ? 3 : 2;
  }
  if (layer1.empty() || layer2.empty())
    return;

  struct Candidate {
    Int_t pulse1;
    Int_t pulse2;
    Double_t dt;
    Double_t dx;
    Double_t dy;
    Double_t cdetScore;
    Double_t ecalResidual;
    Double_t trajectoryResidual;
    Double_t ecalScore;
    Int_t yTopology;
  };
  std::vector<Candidate> candidates;
  candidates.reserve(layer1.size() * layer2.size());
  for (const Int_t pulse1 : layer1) {
    for (const Int_t pulse2 : layer2) {
      const Double_t dt = fPulseLECorrected[pulse2] -
          fPulseLECorrected[pulse1];
      const Double_t dx = fPulseCorrectedX[pulse2] -
          fPulseCorrectedX[pulse1];
      const Double_t dy = fPulseY[pulse2] - fPulseY[pulse1];
      const Bool_t sameSide = std::fabs(dy) <= fPairingDeltaYMax;
      const Double_t alignedProjectedY = 0.5 *
          (fPulseProjectedECalY[pulse1] + fPulseProjectedECalY[pulse2]) +
          fSelectionYResidualOffset;
      const Bool_t oppositeSide = fPairingOppositeSideEnabled &&
          std::fabs(std::fabs(dy) - fPairingOppositeDeltaYCenter) <=
              fPairingOppositeDeltaYTolerance &&
          std::fabs(alignedProjectedY - fPairingOppositeProjectedYCenter) <=
              fPairingOppositeProjectedYMax;
      if (std::fabs(dt - fPairingDeltaTimeCenter) > fPairingDeltaTimeMax ||
          std::fabs(dx) > fPairingDeltaXMax || (!sameSide && !oppositeSide))
        continue;
      const Double_t timePull =
          (dt - fPairingDeltaTimeCenter) / fPairingTimeScale;
      const Double_t xPull = dx / fPairingXScale;
      const Double_t cdetScore = timePull*timePull + xPull*xPull;

      // Compare the inter-layer displacement with the displacement expected
      // from the ECal-to-target ray. Since projected_x = x_ECal*z/z_ECal,
      // this difference is equivalent to comparing the two trajectory slopes
      // without dividing by the small layer separation.
      const Double_t projectedDeltaX = fPulseProjectedECalX[pulse2] -
          fPulseProjectedECalX[pulse1];
      const Double_t trajectoryResidual = dx - projectedDeltaX;
      const Double_t ecalResidual = 0.5 *
          (fPulseECalResidual[pulse1] + fPulseECalResidual[pulse2]);
      const Double_t trajectoryPull =
          (trajectoryResidual - fPairingECalTrajectoryCenter) /
          fPairingECalTrajectoryScale;
      const Double_t timingPull =
          (ecalResidual - fPairingECalTimingCenter) /
          fPairingECalTimingScale;
      const Double_t ecalScore =
          trajectoryPull*trajectoryPull + timingPull*timingPull;
      if (fPairingECalRankEnabled &&
          ecalScore > fPairingECalRadius*fPairingECalRadius)
        continue;
      candidates.push_back({pulse1, pulse2, dt, dx, dy, cdetScore,
                            ecalResidual, trajectoryResidual, ecalScore,
                            oppositeSide && !sameSide ? 1 : 0});
    }
  }

  // Preserve candidate-generation order for the vanishingly rare exact-score
  // tie. This makes the greedy matching deterministic and auditable.
  std::stable_sort(candidates.begin(), candidates.end(),
      [this](const Candidate& left, const Candidate& right) {
        if (fPairingECalRankEnabled && left.ecalScore != right.ecalScore)
          return left.ecalScore < right.ecalScore;
        return left.cdetScore < right.cdetScore;
      });

  // Preserve every ranked hypothesis that passed the detector-local gates.
  // Unlike the legacy pair collection below, candidates may share pulses so
  // that a downstream global ROI algorithm can resolve the ambiguity using
  // information from other detectors.
  for (const Candidate& candidate : candidates) {
    const Double_t time1 = fPulseLECorrected[candidate.pulse1];
    const Double_t time2 = fPulseLECorrected[candidate.pulse2];
    fPairCandidateIndex.push_back(
        static_cast<Int_t>(fPairCandidateIndex.size()));
    fPairCandidatePulseIndexL1.push_back(candidate.pulse1);
    fPairCandidatePulseIndexL2.push_back(candidate.pulse2);
    fPairCandidatePMTL1.push_back(fPulsePMT[candidate.pulse1]);
    fPairCandidatePMTL2.push_back(fPulsePMT[candidate.pulse2]);
    fPairCandidateTimeL1.push_back(time1);
    fPairCandidateTimeL2.push_back(time2);
    fPairCandidateTimeMean.push_back(0.5 * (time1 + time2));
    fPairCandidateDeltaTime.push_back(candidate.dt);
    fPairCandidateDeltaX.push_back(candidate.dx);
    fPairCandidateDeltaY.push_back(candidate.dy);
    fPairCandidateScore.push_back(candidate.cdetScore);
    fPairCandidateECalResidual.push_back(candidate.ecalResidual);
    fPairCandidateTrajectoryResidual.push_back(candidate.trajectoryResidual);
    fPairCandidateECalScore.push_back(candidate.ecalScore);
    fPairCandidateYTopology.push_back(candidate.yTopology);
    fPairCandidateGreedySelected.push_back(0);
    fPairCandidateSelectedPairIndex.push_back(-1);
  }
  if (!fPairCandidateIndex.empty())
    fROICandidateStatus = 1;

  std::vector<Bool_t> used(fPulsePMT.size(), false);
  for (size_t candidateIndex = 0; candidateIndex < candidates.size();
       ++candidateIndex) {
    const Candidate& candidate = candidates[candidateIndex];
    if (used[candidate.pulse1] || used[candidate.pulse2])
      continue;
    used[candidate.pulse1] = true;
    used[candidate.pulse2] = true;
    const Double_t time1 = fPulseLECorrected[candidate.pulse1];
    const Double_t time2 = fPulseLECorrected[candidate.pulse2];
    const Double_t meanTime = 0.5 * (time1 + time2);
    const Int_t selectedPairIndex = static_cast<Int_t>(fPairIndex.size());
    fPairIndex.push_back(selectedPairIndex);
    fPairPulseIndexL1.push_back(candidate.pulse1);
    fPairPulseIndexL2.push_back(candidate.pulse2);
    fPairPMTL1.push_back(fPulsePMT[candidate.pulse1]);
    fPairPMTL2.push_back(fPulsePMT[candidate.pulse2]);
    fPairTimeL1.push_back(time1);
    fPairTimeL2.push_back(time2);
    fPairTimeMean.push_back(meanTime);
    fPairDeltaTime.push_back(candidate.dt);
    fPairDeltaX.push_back(candidate.dx);
    fPairDeltaY.push_back(candidate.dy);
    fPairScore.push_back(candidate.cdetScore);
    fPairECalResidual.push_back(candidate.ecalResidual);
    fPairTrajectoryResidual.push_back(candidate.trajectoryResidual);
    fPairECalScore.push_back(candidate.ecalScore);
    fPairYTopology.push_back(candidate.yTopology);
    fPairCandidateGreedySelected[candidateIndex] = 1;
    fPairCandidateSelectedPairIndex[candidateIndex] = selectedPairIndex;
    if (!fPairingAllowMultiple)
      break;
  }
}

/*
 * FindGoodHit()
 */
Int_t SBSCDet::FindGoodHit(SBSElement *element)
{
  return SBSGenericDetector::FindGoodHit(element);
}

Int_t SBSCDet::CoarseProcess( TClonesArray& tracks )
{
  if(fCoarseProcessed)
    return 0;

  //std::cout << "SBSCDet::CoarseProcess ... " << std::endl;

  // Call the parent class so that it can prepare the data structure on the
  // event it just read from file

  //std::cout << "fNGoodTDChits = " << fNGoodTDChits << std::endl;
  SBSGenericDetector::CoarseProcess(tracks);
  BuildPulseCandidates();
  //std::cout << "Back SBSGenericDetector::CoarseProcess ... fNGoodTDChits = " << fNGoodTDChits << std::endl;

  double x, y, z;
  //double tmin, tmax;

  Int_t nHit = 0;
  SBSCDet_Hit* the_hit = nullptr;

  //fNtrackMatch = 0;

  // Loop through all of the "good" hits coming out of SBSGenericDetector, which has time cuts wide open
  for(int k = 0; k<fNGoodTDChits; k++){
    //tmin = -fElements[fGood.TDCelemID[k]]->TDC()->GetGoodTimeCut();
    //tmax = +fElements[fGood.TDCelemID[k]]->TDC()->GetGoodTimeCut();
    //std::cout << "GetGoodTimeCut() -> " << tmin << " " << tmax << std::endl;

    //double t0 = fElements[fGood.TDCelemID[k]]->TDC()->GetGoodTimeCut();

    //    if(tmin<=fGood.t[k] && fGood.t[k]<=tmax){
    //if (fGood.TDCelemID[k] >= 2688) {
	//std::cout << "Processing good hit " << k << " fHit_tmin = " << fHit_tmin << " fHit_tmax = " << 
    	//fHit_tmax << " le time = " << fGood.t[k] << " te time = " << fGood.t_te[k] << " PMT = " << fGood.TDCelemID[k] << std::endl; 
    //}
    if( fHit_tmin <= fGood.t[k] && fGood.t[k] <= fHit_tmax && fHit_totmin <= fGood.t_ToT[k] && fGood.t_ToT[k] <= fHit_totmax){
      the_hit = new( (*fHits)[nHit++] ) SBSCDet_Hit();

      the_hit->SetPMTNum(fGood.TDCelemID[k]);
      the_hit->SetRow(fGood.TDCrow[k]);
      the_hit->SetCol(fGood.TDCcol[k]);
      the_hit->SetLayer(fGood.TDClayer[k]);
      the_hit->SetTDC_LE(fGood.t[k]);
      the_hit->SetTDC_TE(fGood.t_te[k]);
      the_hit->SetToT(fGood.t_ToT[k]);

      x = (fElements[fGood.TDCelemID[k]])->GetX();
      y = (fElements[fGood.TDCelemID[k]])->GetY();
      z = (fElements[fGood.TDCelemID[k]])->GetZ();
      //std::cout << "PMT = " << fGood.TDCelemID[k] << " X = " << x << " Y = " << y << " Z = " << z << " LE = " << fGood.t[k] << " TE = " << fGood.t_te[k] << " ToT = " << fGood.t_ToT[k] << std::endl;

      the_hit->SetX(x);
      the_hit->SetY(y);
      the_hit->SetZ(z);
    }
  }
  //clustering to be done by dereived class...




// EJB ---- old code below - March 29, 2025 --------------------------------
  // All good hits now defined.  Now determine the position based on
  // time differences between paddles
  // For example:
  //
  //for(int row = 0; row < fNrows; row++) {
  //{
  //  SBSData::TDCHit lPMT = fElements[row][0][0]->TDC()->GetGoodHit();
  //  SBSData::TDCHit rPMT = fElements[row][1][0]->TDC()->GetGoodHit();
  //  diff = lPMT.le.val - rPMT.le.val; // Leading edge difference
  //}

  fCoarseProcessed = 1;

  //std::cout << "End Coarse Process "  << std::endl;
  return 0;
}

Int_t SBSCDet::FineProcess( TClonesArray& tracks )
{

  if(fFineProcessed)
    return 0;

  // Do more detailed processing here.  Parent class does nothing, so no need
  // to call it.
  // We can prepare more detailed output if we want.

  fFineProcessed = 1;
  return 0;
}
