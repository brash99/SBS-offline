//*-- Author :    Andrew Puckett 2025-03-23:

//////////////////////////////////////////////////////////////////////////
//
// SBSGEPRegionOfInterestModule
//
// Grab  ECAL cluster position (and possibly other information)
// Store it  
//
//////////////////////////////////////////////////////////////////////////

#include "SBSGEPRegionOfInterestModule.h"
#include "InterStageModule.h"
#include "THaGlobals.h"
#include "SBSGEPEArm.h" //For the electron arm:
#include "SBSEArm.h" //For the proton arm
#include "SBSECal.h"
#include "SBSCDet.h"
#include "SBSHCal.h"
#include "SBSGEMSpectrometerTracker.h"
#include "SBSGEMPolarimeterTracker.h"
#include "TClonesArray.h"
#include "THaTrack.h"
#include "TMath.h"
#include "TList.h"
#include <cmath>
#include <limits>

namespace {
struct WeightedLineFit {
  double intercept = 0.0;
  double slope = 0.0;
  double chi2 = 0.0;
  int ndf = -1;
  bool valid = false;
  std::vector<double> standardizedResidual;
};

WeightedLineFit FitCoordinate(const std::vector<double>& z,
                              const std::vector<double>& value,
                              const std::vector<double>& sigma)
{
  WeightedLineFit fit;
  if (z.size() < 2 || z.size() != value.size() || z.size() != sigma.size())
    return fit;

  double sw = 0.0, swz = 0.0, swzz = 0.0, swv = 0.0, swzv = 0.0;
  for (size_t i = 0; i < z.size(); ++i) {
    if (!std::isfinite(z[i]) || !std::isfinite(value[i]) ||
        !std::isfinite(sigma[i]) || sigma[i] <= 0.0)
      return fit;
    const double weight = 1.0 / (sigma[i] * sigma[i]);
    sw += weight;
    swz += weight * z[i];
    swzz += weight * z[i] * z[i];
    swv += weight * value[i];
    swzv += weight * z[i] * value[i];
  }
  const double determinant = sw * swzz - swz * swz;
  if (!(determinant > 0.0))
    return fit;

  fit.intercept = (swzz * swv - swz * swzv) / determinant;
  fit.slope = (sw * swzv - swz * swv) / determinant;
  const double varIntercept = swzz / determinant;
  const double covariance = -swz / determinant;
  const double varSlope = sw / determinant;
  fit.standardizedResidual.reserve(z.size());
  for (size_t i = 0; i < z.size(); ++i) {
    const double residual = value[i] - fit.intercept - fit.slope * z[i];
    const double rawPull = residual / sigma[i];
    fit.chi2 += rawPull * rawPull;
    const double weight = 1.0 / (sigma[i] * sigma[i]);
    const double leverage = weight * (varIntercept +
        2.0 * z[i] * covariance + z[i] * z[i] * varSlope);
    const double residualVarianceFraction = 1.0 - leverage;
    fit.standardizedResidual.push_back(residualVarianceFraction > 1.0e-12 ?
        rawPull / std::sqrt(residualVarianceFraction) :
        std::numeric_limits<double>::quiet_NaN());
  }
  fit.ndf = static_cast<int>(z.size()) - 2;
  fit.valid = true;
  return fit;
}

WeightedLineFit FitCoordinateThroughPoint(const std::vector<double>& z,
                                          const std::vector<double>& value,
                                          const std::vector<double>& sigma,
                                          double zAnchor,
                                          double valueAnchor)
{
  WeightedLineFit fit;
  if (z.empty() || z.size() != value.size() || z.size() != sigma.size())
    return fit;

  double numerator = 0.0;
  double denominator = 0.0;
  for (size_t i = 0; i < z.size(); ++i) {
    if (!std::isfinite(z[i]) || !std::isfinite(value[i]) ||
        !std::isfinite(sigma[i]) || sigma[i] <= 0.0)
      return fit;
    const double dz = z[i] - zAnchor;
    const double weight = 1.0 / (sigma[i] * sigma[i]);
    numerator += weight * dz * (value[i] - valueAnchor);
    denominator += weight * dz * dz;
  }
  if (!(denominator > 0.0))
    return fit;

  fit.slope = numerator / denominator;
  fit.intercept = valueAnchor - fit.slope * zAnchor;
  for (size_t i = 0; i < z.size(); ++i) {
    const double pull =
        (value[i] - fit.intercept - fit.slope * z[i]) / sigma[i];
    fit.chi2 += pull * pull;
  }
  fit.ndf = static_cast<int>(z.size()) - 1;
  fit.valid = true;
  return fit;
}
}
//_____________________________________________________________________________
SBSGEPRegionOfInterestModule::SBSGEPRegionOfInterestModule( const char *name, const char *description, Int_t stage ) : InterStageModule(name,description,stage){
  //Constructor; for now, does nothing other than instantiate
  //Default Earm and Parm names:
  fEarmName = "earm";
  fParmName = "sbs";
  fEarmDetName = "ecal";
  fEarmCDetName = "cdet";
  fParmDetName = "gemFT";
  fParmDetNamePol = "gemFPP";
  fParmDetNameCalo = "hcal";
  
  fTestTracks = new TClonesArray("THaTrack",1);

  fTargZ0 = 0.0;
  fSigmaXECal = 0.006;
  fSigmaYECal = 0.006;
  fSigmaXCDet = 0.017973;
  fSigmaYCDet = 0.255;
  
  fDataValid = false; 
}
//_____________________________________________________________________________
SBSGEPRegionOfInterestModule::~SBSGEPRegionOfInterestModule()
{
  //Destructor; for now, does nothing except call THaAnalysisObject::RemoveVariables();
  RemoveVariables();

  delete fTestTracks;
};

//_____________________________________________________________________________
//Clear method: Invoke standard InterStageModule::Clear():
void SBSGEPRegionOfInterestModule::Clear( Option_t *opt )
{
  InterStageModule::Clear(opt);
  //Clear out any other event-level variables here:
  fECALclusterpos_global.SetXYZ(kBig,kBig,kBig);

  fECAL_energy = kBig;
  
  fetheta_central = kBig;
  fephi_central = kBig;
  fEprime_central = kBig;
  fptheta_central = kBig;
  fpphi_central = kBig;
  fPp_central = kBig;
  fECAL_energy = kBig;

  fTestTracks->Clear("C");

  fCDetFound = 0;
  fCDetTimingStatus = 0;
  fCDetROICandidateStatus = 0;
  fCDetNumPulseCandidates = 0;
  fCDetNumPairCandidates = 0;
  fCDetNumSingleCandidates = 0;
  fCDetHypIndex.clear();
  fCDetHypSourceType.clear();
  fCDetHypSourceIndex.clear();
  fCDetHypPulseIndexL1.clear();
  fCDetHypPulseIndexL2.clear();
  fCDetHypNPoints.clear();
  fCDetHypSourceScore.clear();
  fCDetHypX0.clear();
  fCDetHypXSlope.clear();
  fCDetHypXChi2.clear();
  fCDetHypXNDF.clear();
  fCDetHypY0.clear();
  fCDetHypYSlope.clear();
  fCDetHypYChi2.clear();
  fCDetHypYNDF.clear();
  fCDetHypXECalPull.clear();
  fCDetHypXL1Pull.clear();
  fCDetHypXL2Pull.clear();
  fCDetHypYECalPull.clear();
  fCDetHypYL1Pull.clear();
  fCDetHypYL2Pull.clear();
  fCDetHypThetaGlobal.clear();
  fCDetHypPhiGlobal.clear();
  fCDetHypXAtNominalTarget.clear();
  fCDetHypYAtNominalTarget.clear();
  fCDetVertexHypIndex.clear();
  fCDetVertexBin.clear();
  fCDetVertexZ.clear();
  fCDetVertexX.clear();
  fCDetVertexY.clear();
  fCDetVertexXSlope.clear();
  fCDetVertexXChi2.clear();
  fCDetVertexXNDF.clear();
  fCDetVertexYSlope.clear();
  fCDetVertexYResidualL1.clear();
  fCDetVertexYResidualL2.clear();
  fCDetVertexYCompatibleL1.clear();
  fCDetVertexYCompatibleL2.clear();
  fCDetVertexYCompatible.clear();
  fCDetVertexThetaGlobal.clear();
  fCDetVertexPhiGlobal.clear();
}

//_____________________________________________________________________________
Int_t SBSGEPRegionOfInterestModule::DefineVariables( THaAnalysisObject::EMode mode )
{
  // Define/delete event-by-event global variables

  Int_t ret = InterStageModule::DefineVariables(mode);// exports fDataValid etc.
  if( ret )
    return ret;

  RVarDef vars[] = {
    { "xECAL_global", "Global ECAL cluster x (m)", "fECALclusterpos_global.X()" },
    { "yECAL_global", "Global ECAL cluster y (m)", "fECALclusterpos_global.Y()" },
    { "zECAL_global", "Global ECAL cluster z (m)", "fECALclusterpos_global.Z()" },
    { "ECAL_energy", "ECAL best cluster energy (GeV)", "fECAL_energy" },
    { "etheta",  "electron polar angle (rad)", "fetheta_central" },
    { "ephi",  "electron azimuthal angle (rad)", "fephi_central" },
    { "Eprime", "electron expected energy (GeV)", "fEprime_central" },
    { "ptheta", "proton expected polar angle (rad)", "fptheta_central" },
    { "pphi", "proton expected azimuthal angle (rad)", "fpphi_central" },
    { "pp", "proton expected momentum (GeV/c)", "fPp_central" },
    { "xfp0", "predicted X at fp (assuming point target at origin)", "fxfp_central" },
    { "yfp0", "predicted Y at fp (assuming point target at origin)", "fyfp_central" },
    { "xpfp0", "predicted X' at fp (assuming point target at origin)", "fxpfp_central" },
    { "ypfp0", "predicted Y' at fp (assuming point target at origin)", "fypfp_central" },
    { "cdet.found", "CDet detector found by the ROI module", "fCDetFound" },
    { "cdet.timing_status", "CDet timing-calibration status", "fCDetTimingStatus" },
    { "cdet.roi_status", "CDet detector-local ROI candidate status", "fCDetROICandidateStatus" },
    { "cdet.npulse", "Number of complete CDet pulse candidates", "fCDetNumPulseCandidates" },
    { "cdet.npair_candidate", "Number of pre-greedy CDet pair hypotheses", "fCDetNumPairCandidates" },
    { "cdet.nsingle_candidate", "Number of exclusive single-layer CDet hypotheses", "fCDetNumSingleCandidates" },
    { "cdet.hyp.n", "Number of diagnostic CDet electron-ray hypotheses", "GetNumCDetHypotheses()" },
    { "cdet.hyp.index", "Event-local diagnostic hypothesis index", "fCDetHypIndex" },
    { "cdet.hyp.source_type", "Hypothesis source: 1 pair, 2 Layer-1-only, 3 Layer-2-only", "fCDetHypSourceType" },
    { "cdet.hyp.source_index", "Source pair_candidate or single_candidate index", "fCDetHypSourceIndex" },
    { "cdet.hyp.pulse_index_l1", "Source Layer-1 pulse index, or -1", "fCDetHypPulseIndexL1" },
    { "cdet.hyp.pulse_index_l2", "Source Layer-2 pulse index, or -1", "fCDetHypPulseIndexL2" },
    { "cdet.hyp.npoints", "Number of ECal/CDet points in the fit", "fCDetHypNPoints" },
    { "cdet.hyp.source_score", "Detector-local source-candidate score", "fCDetHypSourceScore" },
    { "cdet.hyp.x0", "Electron-ray x intercept at electron-arm z=0 (m)", "fCDetHypX0" },
    { "cdet.hyp.xslope", "Electron-ray dx/dz in the electron-arm frame", "fCDetHypXSlope" },
    { "cdet.hyp.xchi2", "Resolution-weighted x fit chi-square", "fCDetHypXChi2" },
    { "cdet.hyp.xndf", "x fit number of degrees of freedom", "fCDetHypXNDF" },
    { "cdet.hyp.y0", "Electron-ray y intercept at electron-arm z=0 (m)", "fCDetHypY0" },
    { "cdet.hyp.yslope", "Electron-ray dy/dz in the electron-arm frame", "fCDetHypYSlope" },
    { "cdet.hyp.ychi2", "Resolution-weighted y fit chi-square", "fCDetHypYChi2" },
    { "cdet.hyp.yndf", "y fit number of degrees of freedom", "fCDetHypYNDF" },
    { "cdet.hyp.xpull_ecal", "Leverage-corrected ECal x residual, or NaN", "fCDetHypXECalPull" },
    { "cdet.hyp.xpull_l1", "Leverage-corrected Layer-1 CDet x residual, or NaN", "fCDetHypXL1Pull" },
    { "cdet.hyp.xpull_l2", "Leverage-corrected Layer-2 CDet x residual, or NaN", "fCDetHypXL2Pull" },
    { "cdet.hyp.ypull_ecal", "Leverage-corrected ECal y residual, or NaN", "fCDetHypYECalPull" },
    { "cdet.hyp.ypull_l1", "Leverage-corrected Layer-1 CDet y residual, or NaN", "fCDetHypYL1Pull" },
    { "cdet.hyp.ypull_l2", "Leverage-corrected Layer-2 CDet y residual, or NaN", "fCDetHypYL2Pull" },
    { "cdet.hyp.theta_global", "Free detector-only fitted global polar angle (diagnostic, rad)", "fCDetHypThetaGlobal" },
    { "cdet.hyp.phi_global", "Free detector-only fitted global azimuth (diagnostic, rad)", "fCDetHypPhiGlobal" },
    { "cdet.hyp.x_at_ztarg0", "Free detector-only fitted x at nominal target z (diagnostic, m)", "fCDetHypXAtNominalTarget" },
    { "cdet.hyp.y_at_ztarg0", "Free detector-only fitted y at nominal target z (diagnostic, m)", "fCDetHypYAtNominalTarget" },
    { "cdet.vertex.n", "Number of diagnostic hypothesis/target-z associations", "GetNumCDetVertexAssociations()" },
    { "cdet.vertex.hyp_index", "Associated diagnostic hypothesis index", "fCDetVertexHypIndex" },
    { "cdet.vertex.bin", "Existing target-z scan bin", "fCDetVertexBin" },
    { "cdet.vertex.z", "Target-z scan coordinate (m)", "fCDetVertexZ" },
    { "cdet.vertex.x", "Free detector-only ray x at target-z bin (diagnostic, m)", "fCDetVertexX" },
    { "cdet.vertex.y", "Free detector-only ray y at target-z bin (diagnostic, m)", "fCDetVertexY" },
    { "cdet.vertex.xslope", "Best dx/dz constrained through x=0 at this target-z bin", "fCDetVertexXSlope" },
    { "cdet.vertex.xchi2", "x chi-square for the target-z-constrained ray", "fCDetVertexXChi2" },
    { "cdet.vertex.xndf", "x NDF for the target-z-constrained ray", "fCDetVertexXNDF" },
    { "cdet.vertex.yslope", "ECal-to-target dy/dz used only for CDet y compatibility", "fCDetVertexYSlope" },
    { "cdet.vertex.yresidual_l1", "Layer-1 half-bar-center y residual, or NaN (m)", "fCDetVertexYResidualL1" },
    { "cdet.vertex.yresidual_l2", "Layer-2 half-bar-center y residual, or NaN (m)", "fCDetVertexYResidualL2" },
    { "cdet.vertex.ycompatible_l1", "Layer-1 y compatibility: -1 missing, 0 fail, 1 pass", "fCDetVertexYCompatibleL1" },
    { "cdet.vertex.ycompatible_l2", "Layer-2 y compatibility: -1 missing, 0 fail, 1 pass", "fCDetVertexYCompatibleL2" },
    { "cdet.vertex.ycompatible", "All available CDet layers pass half-bar y compatibility", "fCDetVertexYCompatible" },
    { "cdet.vertex.theta_global", "Target-constrained electron-ray global polar angle (rad)", "fCDetVertexThetaGlobal" },
    { "cdet.vertex.phi_global", "Target-constrained electron-ray global azimuth (rad)", "fCDetVertexPhiGlobal" },
    { nullptr }
  };

  return DefineVarsFromList( vars, mode );
}
//_____________________________________________________________________________
Int_t SBSGEPRegionOfInterestModule::ReadRunDatabase( const TDatime &date ){
  //Load beam energy:

  FILE* file = OpenRunDBFile( date );
  if( !file ) return kFileError;

  double ebeamtemp;
  
  const DBRequest req[] = {
    { "ebeam", &ebeamtemp, kDouble, 0, 0, 1 },
    { nullptr }
  };
  Int_t err = LoadDB( file, date, req );
  fclose(file);
  if( err )
    return kInitError;

  //We're neglecting electron mass here (for the purposes of this module it won't matter):
  fBeam4Vect.SetPxPyPzE( 0.0, 0.0, ebeamtemp, ebeamtemp);
  
  return kOK; 
}

Int_t SBSGEPRegionOfInterestModule::ReadDatabase( const TDatime &date ){
  //Load vertex z bin definitions, earm name and parm name:
  FILE* file = OpenFile( date );
  if( !file ){
    std::cerr << "SBSGEPRegionOfInterestModule::ReadDatabase(): database not found!"<< std::endl;
    return kFileError;
  }

  //load vertex z bins for scanning target length:
  const DBRequest request[] = {
    { "nbins_zvertex", &fNbinsVertexZ, kInt, 0, 1, 1 },
    { "zvertex_min", &fVertexZmin, kDouble, 0, 1, 1 },
    { "zvertex_max", &fVertexZmax, kDouble, 0, 1, 1 },
    { "earm_name", &fEarmName, kString, 0, 1, 1 },
    { "parm_name", &fParmName, kString, 0, 1, 1 },
    { "edet_name", &fEarmDetName, kString, 0, 1, 1 },
    { "cdet_name", &fEarmCDetName, kString, 0, 1, 1 },
    { "pdet_name", &fParmDetName, kString, 0, 1, 1 },
    { "pdetpol_name", &fParmDetNamePol, kString, 0, 1, 1 },
    { "pdetcalo_name", &fParmDetNameCalo, kString, 0, 1, 1 },
    { "z0targ", &fTargZ0, kDouble, 0, 1, 1 },
    { "sigma_x_ecal", &fSigmaXECal, kDouble, 0, 1, 1 },
    { "sigma_y_ecal", &fSigmaYECal, kDouble, 0, 1, 1 },
    { "sigma_x_cdet", &fSigmaXCDet, kDouble, 0, 1, 1 },
    { "sigma_y_cdet", &fSigmaYCDet, kDouble, 0, 1, 1 },
    { nullptr }
  };
  
  Int_t status = LoadDB( file, date, request, fPrefix, 1 ); //The "1" after fPrefix means search up the tree
  fclose(file);
  if( status != 0 ){
    return status;
  }

  if (fNbinsVertexZ <= 0 || fSigmaXECal <= 0.0 || fSigmaYECal <= 0.0 ||
      fSigmaXCDet <= 0.0 || fSigmaYCDet <= 0.0) {
    Error(Here("ReadDatabase"),
          "invalid target-z binning or ECal/CDet position uncertainty");
    return kInitError;
  }

  fIsInit = true;

  return kOK;
  
}

//_____________________________________________________________________________
Int_t SBSGEPRegionOfInterestModule::Process( const THaEvData &evdata ){
  //Okay here we go: we've written the code needed to start writing the code.

  THaApparatus *app = 0;

  bool gotEarm = false;
  bool gotParm = false;
  bool gotEdet = false;
  bool gotCDet = false;
  bool gotPdet = false;
  bool gotPdetPol = false;
  bool gotPdetCalo = false;
  
  TIter aiter(gHaApps);

  SBSGEPEArm *Earm = nullptr;
  SBSEArm *Parm = nullptr;

  SBSECal *Edet = nullptr;
  SBSCDet *CDet = nullptr;
  SBSGEMSpectrometerTracker *Pdet = nullptr;
  SBSGEMPolarimeterTracker *PdetPol = nullptr;

  SBSHCal *PdetCalo = nullptr;
  
  while( (app = (THaApparatus*) aiter()) ){
    std::string appname = app->GetName();
    if( app->InheritsFrom("SBSGEPEArm") ){
      if( appname == fEarmName ){
	Earm = dynamic_cast<SBSGEPEArm*>(app);
	gotEarm = true;

	Edet = dynamic_cast<SBSECal*>(Earm->GetDetector(fEarmDetName.c_str()));

	if( Edet ) gotEdet = true;

	CDet = dynamic_cast<SBSCDet*>(Earm->GetDetector(fEarmCDetName.c_str()));
	if( CDet ) gotCDet = true;
      }
    }
    if( app->InheritsFrom("SBSEArm") ){
      if( appname == fParmName ){
	Parm = dynamic_cast<SBSEArm*>(app);
	gotParm = true;

	Pdet = dynamic_cast<SBSGEMSpectrometerTracker*>(Parm->GetDetector(fParmDetName.c_str()));
	if( Pdet ) gotPdet = true;

	PdetPol = dynamic_cast<SBSGEMPolarimeterTracker*>(Parm->GetDetector(fParmDetNamePol.c_str()));

	if( PdetPol ) gotPdetPol = true;
	
	PdetCalo = dynamic_cast<SBSHCal*>(Parm->GetDetector(fParmDetNameCalo.c_str()));
	if( PdetCalo ) gotPdetCalo = true;
	
      }
    }
  }

  fCDetFound = gotCDet ? 1 : 0;
  if (CDet) {
    fCDetTimingStatus = CDet->GetTimingStatus();
    fCDetROICandidateStatus = CDet->GetROICandidateStatus();
    fCDetNumPulseCandidates = CDet->GetNumPulseCandidates();
    fCDetNumPairCandidates = CDet->GetNumLayerPairCandidates();
    fCDetNumSingleCandidates = CDet->GetNumSingleLayerCandidates();
  }

  if( !gotParm || !gotEarm || !gotEdet || !gotPdet || !gotPdetPol || !gotPdetCalo ){
    std::cout << "Error: missing Earm and/or Parm and/or Edet and/or Pdet and/or PdetPol and/or PdetCalo! (gotEarm, gotParm, gotEdet, gotPdet, gotPdetPol, gotPdetCalo)=(" << gotEarm << ", " << gotParm << ", "
	      << gotEdet << ", " << gotPdet << ", " << gotPdetPol
	      << ", " << gotPdetCalo << ")" << std::endl;
    fDataValid = false;
    return 0;
  }

  if( !Parm->GetGEPtrackingMode() ){ //proton arm must have "GEP tracking mode" set to use this module
    fDataValid = false;
    return 0;
  }
  
  //If we reached this point, then we can grab E arm and P arm info:

  //Always clear out proton arm front tracker constraint points before evaluating ROI:
  
  Pdet->ClearConstraints();
  PdetPol->ClearConstraints();

  //Require at least an HCAL cluster and an E Arm "track":
  if( PdetCalo->GetNclust() <= 0 || Earm->GetNTracks() <= 0 ) return 0;

  auto HCalClusters = PdetCalo->GetClusters();

  int ibest_hcal = PdetCalo->GetBestClusterIndex();

  //The following are HCAL cluster positions in the GEM FP coordinate system!
  double xHCAL = HCalClusters[ibest_hcal]->GetX() + PdetCalo->GetOrigin().X();
  double yHCAL = HCalClusters[ibest_hcal]->GetY() + PdetCalo->GetOrigin().Y();
  double zHCAL = PdetCalo->GetOrigin().Z();
  
  double ThetaEarm = Earm->GetThetaGeo(); //E arm is ordinarily on beam left, so this angle SHOULD be positive
  double ThetaParm = Parm->GetThetaGeo(); //P arm is ordinarily on beam right, so this angle SHOULD be negative

  // The following lines assume ThetaEarm > 0 for beam left:
  TVector3 Earm_zaxis( sin(ThetaEarm), 0.0, cos(ThetaEarm) );
  TVector3 Earm_xaxis(0,-1,0); //TRANSPORT system; +x = down
  TVector3 Earm_yaxis = Earm_zaxis.Cross( Earm_xaxis ).Unit();
  
  //The following lines assume ThetaParm < 0 for beam right:
  TVector3 Parm_zaxis( sin(ThetaParm), 0.0, cos(ThetaParm) );
  TVector3 Parm_xaxis( 0, -1, 0 ); //TRANSPORT system; +x = down
  TVector3 Parm_yaxis = Parm_zaxis.Cross( Parm_xaxis ).Unit();

  //Set HCAL global cluster position. Not yet clear whether and/or how we will use this:
  fHCALclusterpos_global = Parm->GetHCALdist() * Parm_zaxis +
    HCalClusters[ibest_hcal]->GetX() * Parm_xaxis +
    HCalClusters[ibest_hcal]->GetY() * Parm_yaxis; 

  
  
  //Grab the E arm "track" 
  //  if( Earm->GetNTracks() >= 1 ){
  TClonesArray *EarmTracks = Earm->GetTracks();
  
  THaTrack *EarmTrack = ( (THaTrack*) (*EarmTracks)[0] );

  double xclust = EarmTrack->GetX();
  double yclust = EarmTrack->GetY();
  double ECALdist = Earm->GetECalDist() + Edet->GetOrigin().Z();

  // Note: The "E arm track" X and Y positions have already been offset by earm.ecal.position!
  // However, "GetECalDist()" returns the "ecaldist" parameter defined via the run DB, and it is NOT corrected according to earm.ecal.position!
  // Thus, the "track" angles here would be inconsistent with those defined in SBSGEPEArm unless we also offset ECALdist here!
  // For consistency, let's ALSO offset ECALdist by the Z position above!
  // NOTE also, this means that earm.ecal.position is always to be interpreted as a "small" offset from nominal or "ideal"
  
  fECAL_energy = EarmTrack->GetEnergy();
  
  //TVector3 ECALpos(xclust,yclust,ECALdist);
  TVector3 ECALpos_global = xclust * Earm_xaxis + yclust * Earm_yaxis + ECALdist * Earm_zaxis;
  
  fECALclusterpos_global = ECALpos_global;

  // Build read-only diagnostic electron-ray hypotheses from the complete
  // detector-local CDet candidate collections. This deliberately does not
  // modify the GEM constraints below.
  if (CDet && CDet->GetTimingStatus() == 2) {
    const auto appendHypothesis = [&](Int_t sourceType, Int_t sourceIndex,
                                      Int_t pulseIndexL1, Int_t pulseIndexL2,
                                      Double_t sourceScore) {
      std::vector<double> z{ECALdist};
      std::vector<double> x{xclust};
      std::vector<double> y{yclust};
      std::vector<double> sigmaX{fSigmaXECal};
      std::vector<double> sigmaY{fSigmaYECal};
      const double missing = std::numeric_limits<double>::quiet_NaN();
      double zL1 = missing, xL1 = missing, yL1 = missing;
      double zL2 = missing, xL2 = missing, yL2 = missing;
      int pointIndexL1 = -1, pointIndexL2 = -1;

      const auto appendPulse = [&](Int_t pulseIndex) {
        if (pulseIndex < 0)
          return true;
        SBSCDet::PulseCandidate pulse;
        if (!CDet->GetPulseCandidate(pulseIndex, pulse) ||
            !pulse.calibrationValid || !std::isfinite(pulse.correctedX) ||
            !std::isfinite(pulse.y) || !std::isfinite(pulse.z))
          return false;
        z.push_back(pulse.z);
        x.push_back(pulse.correctedX);
        y.push_back(pulse.y);
        sigmaX.push_back(fSigmaXCDet);
        sigmaY.push_back(fSigmaYCDet);
        if (pulse.layer == 0) {
          pointIndexL1 = static_cast<int>(z.size()) - 1;
          zL1 = pulse.z;
          xL1 = pulse.correctedX;
          yL1 = pulse.y;
        } else if (pulse.layer == 1) {
          pointIndexL2 = static_cast<int>(z.size()) - 1;
          zL2 = pulse.z;
          xL2 = pulse.correctedX;
          yL2 = pulse.y;
        }
        return true;
      };

      if (!appendPulse(pulseIndexL1) || !appendPulse(pulseIndexL2))
        return;

      const WeightedLineFit xfit = FitCoordinate(z, x, sigmaX);
      const WeightedLineFit yfit = FitCoordinate(z, y, sigmaY);
      if (!xfit.valid || !yfit.valid)
        return;

      const Int_t hypothesisIndex =
          static_cast<Int_t>(fCDetHypIndex.size());
      const TVector3 directionGlobal =
          (xfit.slope * Earm_xaxis + yfit.slope * Earm_yaxis + Earm_zaxis)
              .Unit();

      fCDetHypIndex.push_back(hypothesisIndex);
      fCDetHypSourceType.push_back(sourceType);
      fCDetHypSourceIndex.push_back(sourceIndex);
      fCDetHypPulseIndexL1.push_back(pulseIndexL1);
      fCDetHypPulseIndexL2.push_back(pulseIndexL2);
      fCDetHypNPoints.push_back(static_cast<Int_t>(z.size()));
      fCDetHypSourceScore.push_back(sourceScore);
      fCDetHypX0.push_back(xfit.intercept);
      fCDetHypXSlope.push_back(xfit.slope);
      fCDetHypXChi2.push_back(xfit.chi2);
      fCDetHypXNDF.push_back(xfit.ndf);
      fCDetHypY0.push_back(yfit.intercept);
      fCDetHypYSlope.push_back(yfit.slope);
      fCDetHypYChi2.push_back(yfit.chi2);
      fCDetHypYNDF.push_back(yfit.ndf);
      fCDetHypXECalPull.push_back(xfit.standardizedResidual[0]);
      fCDetHypXL1Pull.push_back(pointIndexL1 >= 0 ?
          xfit.standardizedResidual[pointIndexL1] : missing);
      fCDetHypXL2Pull.push_back(pointIndexL2 >= 0 ?
          xfit.standardizedResidual[pointIndexL2] : missing);
      fCDetHypYECalPull.push_back(yfit.standardizedResidual[0]);
      fCDetHypYL1Pull.push_back(pointIndexL1 >= 0 ?
          yfit.standardizedResidual[pointIndexL1] : missing);
      fCDetHypYL2Pull.push_back(pointIndexL2 >= 0 ?
          yfit.standardizedResidual[pointIndexL2] : missing);
      fCDetHypThetaGlobal.push_back(directionGlobal.Theta());
      fCDetHypPhiGlobal.push_back(directionGlobal.Phi());
      fCDetHypXAtNominalTarget.push_back(
          xfit.intercept + xfit.slope * fTargZ0);
      fCDetHypYAtNominalTarget.push_back(
          yfit.intercept + yfit.slope * fTargZ0);

      const double zbinwidth =
          (fVertexZmax - fVertexZmin) / double(fNbinsVertexZ);
      for (Int_t ibin = 0; ibin < fNbinsVertexZ; ++ibin) {
        const double zvertex = fVertexZmin + (ibin + 0.5) * zbinwidth;
        const WeightedLineFit constrainedX =
            FitCoordinateThroughPoint(z, x, sigmaX, zvertex, 0.0);
        const double ySlope = yclust / (ECALdist - zvertex);
        const double yResidualL1 = std::isfinite(zL1) ?
            yL1 - ySlope * (zL1 - zvertex) : missing;
        const double yResidualL2 = std::isfinite(zL2) ?
            yL2 - ySlope * (zL2 - zvertex) : missing;
        const Int_t yCompatibleL1 = std::isfinite(yResidualL1) ?
            (std::abs(yResidualL1) <= fSigmaYCDet ? 1 : 0) : -1;
        const Int_t yCompatibleL2 = std::isfinite(yResidualL2) ?
            (std::abs(yResidualL2) <= fSigmaYCDet ? 1 : 0) : -1;
        const Int_t yCompatible =
            (yCompatibleL1 != 0 && yCompatibleL2 != 0) ? 1 : 0;
        const TVector3 constrainedDirectionGlobal =
            (constrainedX.slope * Earm_xaxis + ySlope * Earm_yaxis +
             Earm_zaxis).Unit();

        fCDetVertexHypIndex.push_back(hypothesisIndex);
        fCDetVertexBin.push_back(ibin);
        fCDetVertexZ.push_back(zvertex);
        fCDetVertexX.push_back(xfit.intercept + xfit.slope * zvertex);
        fCDetVertexY.push_back(yfit.intercept + yfit.slope * zvertex);
        fCDetVertexXSlope.push_back(constrainedX.slope);
        fCDetVertexXChi2.push_back(constrainedX.chi2);
        fCDetVertexXNDF.push_back(constrainedX.ndf);
        fCDetVertexYSlope.push_back(ySlope);
        fCDetVertexYResidualL1.push_back(yResidualL1);
        fCDetVertexYResidualL2.push_back(yResidualL2);
        fCDetVertexYCompatibleL1.push_back(yCompatibleL1);
        fCDetVertexYCompatibleL2.push_back(yCompatibleL2);
        fCDetVertexYCompatible.push_back(yCompatible);
        fCDetVertexThetaGlobal.push_back(constrainedDirectionGlobal.Theta());
        fCDetVertexPhiGlobal.push_back(constrainedDirectionGlobal.Phi());
      }
    };

    for (Int_t i = 0; i < CDet->GetNumLayerPairCandidates(); ++i) {
      SBSCDet::LayerPair pair;
      if (CDet->GetLayerPairCandidate(i, pair))
        appendHypothesis(1, pair.index, pair.pulseIndexL1,
                         pair.pulseIndexL2, pair.ecalScore);
    }

    for (Int_t i = 0; i < CDet->GetNumSingleLayerCandidates(); ++i) {
      SBSCDet::SingleLayerCandidate single;
      if (!CDet->GetSingleLayerCandidate(i, single))
        continue;
      appendHypothesis(single.layer == 0 ? 2 : 3, single.index,
                       single.layer == 0 ? single.pulseIndex : -1,
                       single.layer == 1 ? single.pulseIndex : -1,
                       single.score);
    }
  }
  
  Pdet->SetECALpos( ECALpos_global ); //For implementation of "elastic constraint" within track-finding in the FT
  
  TVector3 vertex_central(0,0,fTargZ0);
  
  // Central ECAL direction:
  TVector3 ECALdir_global = (ECALpos_global - vertex_central).Unit();
  
  double ebeam = fBeam4Vect.E();
  double Mp = fmass_proton_GeV;
  
  Pdet->SetBeamE( ebeam );
  //Unclear whether and how we'll use this for the back tracker:
  PdetPol->SetBeamE( ebeam );
  PdetPol->SetECALpos( ECALpos_global );
    
  fetheta_central = ECALdir_global.Theta();
  fephi_central = ECALdir_global.Phi();
  fEprime_central = ebeam/(1.0+ebeam/Mp*(1.0-cos(fetheta_central)));

  double Q2 = 2.0*ebeam*fEprime_central*(1.0-cos(fetheta_central));
  double tau = Q2/(4.0*Mp*Mp);
  fPp_central = sqrt(Q2*(1.0+tau)); // = sqrt(nu^2 + 2M nu)
  fptheta_central = acos( (ebeam-fEprime_central*cos(fetheta_central))/fPp_central );
  fpphi_central = fephi_central + TMath::Pi();

  // std::cout << "SBSGEMRegionOfInterestModule: multi tracks enabled = "
  // 	      << Pdet->MultiTracksEnabled() << std::endl;
  //Get constraint point offsets for centering:

  //Calculate "central" expected proton track regardless of whether we're using the z-vertex binning:

  // This calculation assumes a point target at the origin:
  TVector3 pnhat_central( sin(fptheta_central)*cos(fpphi_central),sin(fptheta_central)*sin(fpphi_central),cos(fptheta_central));
  TVector3 ProtonMomentum = fPp_central * pnhat_central;
  double raytemp[6];
    
  TVector3 vdummy(0,0,fTargZ0);
  TVector3 dummy;
  Parm->LabToTransport( vdummy, ProtonMomentum, dummy, raytemp );
    
  double xptar = raytemp[1];
  double yptar = raytemp[3];
  double xtar = raytemp[0];
  double ytar = raytemp[2];
    
  int itrack = fTestTracks->GetLast()+1;
  THaTrack *Ttemp = new( (*fTestTracks)[itrack] ) THaTrack();
    
  Ttemp->SetTarget( xtar, ytar, xptar, yptar );
  Ttemp->SetMomentum( fPp_central );
    
  Parm->CalcFpCoords( Ttemp );
    
  Ttemp->Set( Ttemp->GetDX(), Ttemp->GetDY(), Ttemp->GetDTheta(), Ttemp->GetDPhi() );
    
  double xfp = Ttemp->GetX();
  double yfp = Ttemp->GetY();
  double xpfp = Ttemp->GetTheta();
  double ypfp = Ttemp->GetPhi();
    
  //set output variables:
  fxfp_central = xfp;
  fyfp_central = yfp;
  fxpfp_central = xpfp;
  fypfp_central = ypfp;

  //Now calculate "central" values of front and back constraints:
  
  
  double x0fcp_front = Parm->GetFrontConstraintX0(0);
  double y0fcp_front = Parm->GetFrontConstraintY0(0);
  double x0bcp_front = Parm->GetBackConstraintX0(0);
  double y0bcp_front = Parm->GetBackConstraintY0(0);

  double x0fcp_back = Parm->GetFrontConstraintX0(1);
  double y0fcp_back = Parm->GetFrontConstraintY0(1);
  double x0bcp_back = Parm->GetBackConstraintX0(1);
  double y0bcp_back = Parm->GetBackConstraintY0(1);

  // Calculate all these once, then use them to set front and back tracker
  // constraints IFF the "multi-tracks" flag is set: 
  std::vector<TVector3> FCPfront(fNbinsVertexZ), FCPback(fNbinsVertexZ),
    BCPfront(fNbinsVertexZ), BCPback(fNbinsVertexZ);

  TVector3 vertex;
  double zbinwidth = (fVertexZmax - fVertexZmin)/double(fNbinsVertexZ);
  
  //if( Pdet->MultiTracksEnabled() ){ //use multiple constraint points:
      
  //      std::cout << "GEP region of interest vertex scan:" << std::endl;
  for( int ibin=0; ibin<fNbinsVertexZ; ibin++ ){
    
    //Set vertex assumption:
    vertex.SetXYZ(0.0,0.0,fVertexZmin + (ibin+0.5)*zbinwidth);
    
    // std::cout << "(ibin, zvertex)=(" << ibin << ", " << vertex.Z() << ")"
    // 	  << std::endl;
    
    //Now calculate electron scattering angle and the rest of e and p kinematic variables:
    TVector3 enhat = (ECALpos_global-vertex).Unit();
    
    double etheta = enhat.Theta();
    double ephi = enhat.Phi();
    double eprime = ebeam/(1.0+ebeam/Mp*(1.0-cos(etheta)));
    Q2 = 2.0*ebeam*eprime*(1.0-cos(etheta));
    tau = Q2/(4.0*Mp*Mp);
    
    double pp = sqrt(Q2*(1.0+tau));
    double ptheta = acos( (ebeam-eprime*cos(etheta))/pp );
    double pphi = ephi + TMath::Pi();
    
    //Proton direction (unit vector):
    TVector3 pnhat(sin(ptheta)*cos(pphi),sin(ptheta)*sin(pphi),cos(ptheta));
    
    ProtonMomentum = pp*pnhat;
    //Next we need to calculate this in SBS (Parm) transport coordinates:
    
    //	double raytemp[6];
    
    TVector3 tvert; 
    
    Parm->LabToTransport( vertex, ProtonMomentum, tvert, raytemp );
    //	TVector3 pnhat_SBS( pnhat.Dot( Parm_xaxis ), pnhat.Dot( Parm_yaxis ), pnhat.Dot( Parm_zaxis ) );
    //double xptar_p = pnhat_SBS.X()/pnhat_SBS.Z();
    //double yptar_p = pnhat_SBS.Y()/pnhat_SBS.Z();
    
    xptar = raytemp[1];
    yptar = raytemp[3];
    xtar = raytemp[0];
    ytar = raytemp[2];
    
    //Now with these quantities calculated, we are able to use the forward optics matrix to predict the FP track:
    
    itrack = fTestTracks->GetLast() + 1;
    
    //Initialize with the default constructor (no arguments);
    Ttemp = new( (*fTestTracks)[itrack] ) THaTrack();
    
    Ttemp->SetTarget( xtar, ytar, xptar, yptar );
    Ttemp->SetMomentum( pp );
    
    //Track to focal plane:
    Parm->CalcFpCoords( Ttemp );
    
    //Although it doesn't matter that much, set the "regular" FP coordinates based on the "detector" coordinates:
    Ttemp->Set( Ttemp->GetDX(), Ttemp->GetDY(), Ttemp->GetDTheta(), Ttemp->GetDPhi() );
    //After the line above, the "det" coordinates are the same as the "regular" coordinates.
    //We could, of course, just grab the "det" coordinates directly:
    
    xfp = Ttemp->GetX();
    yfp = Ttemp->GetY();
    xpfp = Ttemp->GetTheta();
    ypfp = Ttemp->GetPhi();
    
    //Now we add front and back constraint points based on these calculated track parameters.
    // What Z value should we assume for the back constraint point? For the front it's easy.
    // For the back 
    
    double zfcp_front = 0.0;
    double zbcp_front = Pdet->GetZmaxLayer() + 0.05; //the 0.05 here (5 cm past the back GEM layer) is arbitrary... may need adjustment

    FCPfront[ibin].SetXYZ( xfp + xpfp*zfcp_front + x0fcp_front, yfp + ypfp*zfcp_front + y0fcp_front, zfcp_front );
    BCPfront[ibin].SetXYZ( xfp + xpfp*zbcp_front + x0bcp_front, yfp + ypfp*zbcp_front + y0bcp_front, zbcp_front );
    
    double zfcp_back = Parm->GetAnalyzerZ0(); //project to midpoint of analyzer
    double zbcp_back = zHCAL; //This is the HCAL "origin" z coordinate (about 6.4 m downstream from GEM/optics origin)

    FCPback[ibin].SetXYZ( xfp + xpfp*zfcp_back + x0fcp_back, yfp + ypfp*zfcp_back + y0fcp_back, zfcp_back );
    BCPback[ibin].SetXYZ( xHCAL + x0bcp_back, yHCAL + y0bcp_back, zHCAL );
    
    // Pdet->SetFrontConstraintPoint( xfp + x0fcp, yfp + y0fcp, 0.0 );
    // Pdet->SetBackConstraintPoint( xfp + xpfp * zback + x0bcp, yfp + ypfp*zback + y0bcp, zback );
    
    //   //Also add constraint points to diagnostic outputs:
    // Parm->AddFrontConstraintPoint( xfp + x0fcp, yfp + y0fcp, 0.0 );
    // Parm->AddBackConstraintPoint( xfp + xpfp * zback + x0bcp, yfp + ypfp*zback + y0bcp, zback );
    
    //   fDataValid = true;
    //Now, IF the multi-track search is enabled (and if the code is written correctly), this should be sufficient
  } //end loop over z vertex bins 

  if( Pdet->MultiTracksEnabled() ){

    for( int ibin=0; ibin<fNbinsVertexZ; ibin++ ){
      Pdet->SetFrontConstraintPoint( FCPfront[ibin] );
      Pdet->SetBackConstraintPoint( BCPfront[ibin] );
      //The following lines are just filling diagnostic output variables:
      //Parm->AddFrontConstraintPoint( FCPfront[ibin] );
      //Parm->AddBackConstraintPoint( BCPfront[ibin] );
    }
    
  } else { // end if multi-tracks enabled)
      
    double zfront = 0.0;
    double zback = Pdet->GetZmaxLayer() + 0.05;

    TVector3 fcptemp( fxfp_central + fxpfp_central*zfront + x0fcp_front, fyfp_central + fypfp_central*zfront + y0fcp_front, zfront );
    TVector3 bcptemp( fxfp_central + fxpfp_central*zback + x0bcp_front, fyfp_central+fypfp_central*zback + y0bcp_front, zback );
    
    Pdet->SetFrontConstraintPoint( fcptemp );
    Pdet->SetBackConstraintPoint( bcptemp );

    // Parm->AddFrontConstraintPoint( fcptemp );
    // Parm->AddBackConstraintPoint( bcptemp );
  }

  if( PdetPol->MultiTracksEnabled() ){ //Then do a z-vertex scan:
    for( int ibin=0; ibin<fNbinsVertexZ; ibin++ ){
      PdetPol->SetFrontConstraintPoint( FCPback[ibin] );
      PdetPol->SetBackConstraintPoint( BCPback[ibin] );

      //Save these for later; after all tracking decisions are finalized
      //Parm->AddFrontConstraintPoint( FCPback[ibin] );
      //Parm->AddBackConstraintPoint( BCPback[ibin] );
    }
  } else {
    double zfront = Parm->GetAnalyzerZ0();
    double zback = zHCAL;

    TVector3 fcptemp( fxfp_central + fxpfp_central*zfront + x0fcp_back,
		      fyfp_central + fypfp_central*zfront + y0fcp_back,
		      zfront );
    TVector3 bcptemp( fxfp_central + fxpfp_central*zback + x0bcp_back,
		      fyfp_central + fypfp_central*zback + y0bcp_back,
		      zback );
    PdetPol->SetFrontConstraintPoint( fcptemp ); 
    PdetPol->SetBackConstraintPoint( bcptemp );

    //Save these for later; after all tracking decisions are finalized.
    //Parm->AddFrontConstraintPoint( fcptemp );
    //Parm->AddBackConstraintPoint( bcptemp );
    
  }

  fDataValid = true;
  
  return 0;
}

ClassImp(SBSGEPRegionOfInterestModule);
