#ifndef SBSGEPRegionOfInterestModule_h_
#define SBSGEPRegionOfInterestModule_h_

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is adapted from Podd_TimeCorrectionModule example:
//
// Its role is to grab cluster position (and possibly other variables) from SBSGEPEArm
// after the CoarseReconstruct stage for all detectors.
// It then populates a list of front and back constraint points and widths for the GEP
// front tracker
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "InterStageModule.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <vector>

class TClonesArray;
class THaTrack;
//class InterStageModule;

using namespace Podd;

class SBSGEPRegionOfInterestModule : public InterStageModule {
public:
  SBSGEPRegionOfInterestModule( const char *name, const char *description, Int_t stage );
  virtual ~SBSGEPRegionOfInterestModule();

  virtual void Clear( Option_t *opt="" );
  virtual Int_t Process( const THaEvData & );

  //I don't think we need a custom Init method here:
  //virtual EStatus Init( const TDatime& date );

  Double_t GetXfpCentral() const { return fxfp_central; }
  Double_t GetYfpCentral() const { return fyfp_central; }
  Double_t GetThfpCentral() const { return fxpfp_central; }
  Double_t GetPhfpCentral() const { return fypfp_central; }
  Int_t GetNumCDetHypotheses() const
    { return static_cast<Int_t>(fCDetHypIndex.size()); }
  Int_t GetNumCDetVertexAssociations() const
    { return static_cast<Int_t>(fCDetVertexHypIndex.size()); }
  
protected:

  virtual Int_t  DefineVariables( EMode mode = kDefine );
  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  ReadRunDatabase( const TDatime& date ); //This is to load beam energy (redundant, I know, but whatever)

  //no need for this yet
  //  virtual Int_t   Begin( THaRunBase* r=0 );
  
  //constant (per-run) parameters: 
  //-----------------------------------------------------------------------------------------------------------------------
  //This should be loaded from gHaRun->GetParameters()->GetBeamE(); alternatively we could just load it from the database?  
  TLorentzVector fBeam4Vect; //Where is the most convenient place to get the beam energy from? --> Run database 

  const double fmass_proton_GeV = 0.93827208816;
  //Define z vertex bins for scanning target extent (these should be loaded from regular DB):
  Int_t fNbinsVertexZ;
  Double_t fVertexZmin;
  Double_t fVertexZmax;

  Double_t fTargZ0;

  // Working measurement uncertainties for diagnostic electron-ray fits.
  Double_t fSigmaXECal;
  Double_t fSigmaYECal;
  Double_t fSigmaXCDet;
  Double_t fSigmaYCDet;
  
  // Names of Earm and Parm: read from DB; I don't have a strong preference for
  // how to store these; might as well use std::string 
  std::string fEarmName;
  std::string fParmName;

  std::string fEarmDetName;
  std::string fEarmCDetName;
  std::string fParmDetName;
  std::string fParmDetNamePol;
  std::string fParmDetNameCalo;
  
  //We might as well store spectrometer 3-vectors here, or would that be redundant with the ones in the spectrometer classes? 

  //variable (per-event) parameters:
  //------------------------------------------------------------------------------------------------------------------------------
  TVector3 fECALclusterpos_global; //ECAL cluster position in "global" Hall A Coordinates (+x to beam left, +y up, +z along beam) 
  TVector3 fHCALclusterpos_global; //HCAL cluster position in "global" Hall A Coordinates (unclear as of yet whether and how we will directly use this)
  
  Double_t fECAL_energy;
  
  //Electron and proton final-state kinematics from ECAL cluster pos for point-target assumption:
  Double_t fetheta_central;
  Double_t fephi_central;
  Double_t fEprime_central;
  Double_t fptheta_central;
  Double_t fpphi_central;
  Double_t fPp_central;

  //Add variables to define the "central" elastically scattered proton ray (assuming point target at the origin):
  Double_t fxfp_central;
  Double_t fyfp_central;
  Double_t fxpfp_central;
  Double_t fypfp_central;

  // Read-only CDet discovery and diagnostic electron-ray hypotheses. These
  // outputs do not alter the existing GEM constraint families.
  Int_t fCDetFound;
  Int_t fCDetTimingStatus;
  Int_t fCDetROICandidateStatus;
  Int_t fCDetNumPulseCandidates;
  Int_t fCDetNumPairCandidates;
  Int_t fCDetNumSingleCandidates;

  std::vector<Int_t> fCDetHypIndex;
  std::vector<Int_t> fCDetHypSourceType;
  std::vector<Int_t> fCDetHypSourceIndex;
  std::vector<Int_t> fCDetHypPulseIndexL1;
  std::vector<Int_t> fCDetHypPulseIndexL2;
  std::vector<Int_t> fCDetHypYTopology;
  std::vector<Int_t> fCDetHypNPoints;
  std::vector<Double_t> fCDetHypSourceScore;
  std::vector<Double_t> fCDetHypX0;
  std::vector<Double_t> fCDetHypXSlope;
  std::vector<Double_t> fCDetHypXChi2;
  std::vector<Int_t> fCDetHypXNDF;
  std::vector<Double_t> fCDetHypY0;
  std::vector<Double_t> fCDetHypYSlope;
  std::vector<Double_t> fCDetHypYChi2;
  std::vector<Int_t> fCDetHypYNDF;
  std::vector<Double_t> fCDetHypXECalPull;
  std::vector<Double_t> fCDetHypXL1Pull;
  std::vector<Double_t> fCDetHypXL2Pull;
  std::vector<Double_t> fCDetHypYECalPull;
  std::vector<Double_t> fCDetHypYL1Pull;
  std::vector<Double_t> fCDetHypYL2Pull;
  std::vector<Double_t> fCDetHypThetaGlobal;
  std::vector<Double_t> fCDetHypPhiGlobal;
  std::vector<Double_t> fCDetHypXAtNominalTarget;
  std::vector<Double_t> fCDetHypYAtNominalTarget;

  // Flattened hypothesis/target-z-bin associations.
  std::vector<Int_t> fCDetVertexHypIndex;
  std::vector<Int_t> fCDetVertexBin;
  std::vector<Double_t> fCDetVertexZ;
  std::vector<Double_t> fCDetVertexX;
  std::vector<Double_t> fCDetVertexY;
  std::vector<Double_t> fCDetVertexXSlope;
  std::vector<Double_t> fCDetVertexXChi2;
  std::vector<Int_t> fCDetVertexXNDF;
  std::vector<Double_t> fCDetVertexYSlope;
  std::vector<Double_t> fCDetVertexYResidualL1;
  std::vector<Double_t> fCDetVertexYResidualL2;
  std::vector<Int_t> fCDetVertexYCompatibleL1;
  std::vector<Int_t> fCDetVertexYCompatibleL2;
  std::vector<Int_t> fCDetVertexYSeamCompatible;
  std::vector<Int_t> fCDetVertexYCompatible;
  std::vector<Double_t> fCDetVertexThetaGlobal;
  std::vector<Double_t> fCDetVertexPhiGlobal;
  
  TClonesArray *fTestTracks;
  
  // We may want to add some CDET-related info here once we understand what's going on there:
  
  
  ClassDef(SBSGEPRegionOfInterestModule,0)
};

#endif
