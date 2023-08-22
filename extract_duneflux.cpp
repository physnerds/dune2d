
#include "extract_duneflux.h"



//Rotation Functions......
TVector3 RotatePhi(TVector3 det_loc){
  double cosPhi = std::cos(phi_det_to_beam);
  double sinPhi = std::sin(phi_det_to_beam);

  return{ cosPhi*det_loc(0) - sinPhi*det_loc(1),
          sinPhi*det_loc(0)+cosPhi*det_loc(1),
          det_loc(2)

  };
}

TVector3 RotateTheta(TVector3 det_loc){
  double cosTheta = std::cos(theta_det_to_beam);
  double sinTheta = std::sin(theta_det_to_beam);

  TVector3 temp_vec(cosTheta*det_loc(0)+sinTheta*det_loc(2),
          det_loc(1),
          -sinTheta*det_loc(0)+cosTheta*det_loc(2));
  //std::cout<<theta_det_to_beam<<" "<<det_loc(0)<<" "<<det_loc(1)<<" "<<det_loc(2)<<" "<<temp_vec(0)<<" "<<temp_vec(1)<<" "<<temp_vec(2)<<std::endl;
  return temp_vec;
}

TVector3 RotateToBeam(TVector3 det_loc){
 // TVector temp_vec = RotatePhi(det_loc);
  return { RotatePhi(RotateTheta(det_loc))};
  }

//=============================================================================
// CalcEnuWgt Wrapper
//=============================================================================
// Wrapper calcEnuWgt function
//   uses User2BeamPos, which converts it to Beam coords (and into cm)
//   Feeds beam pos (in cm) to bsim::calcEnuWgt which returns Enu and wgt_xy
//
// bsim::calcEnuWgt expects a position in Beam coords (cm)
//
// Re: User2BeamPos converting (m) to (cm), see 
// GDk2NuFluxXMLHelper::ParseParamSet + XML file "setting user units"
void calcEnuWgt(double x, double y, double z, 
                  bsim::Dk2Nu* dk2nu, 
                  double& Enu, 
                  double& wgt_xy)
{
  TLorentzVector X4det, X4beam;
  TVector3              Xbeam;
  X4det.SetXYZT(x*100, y*100, z*100+57400.0, 0); //note that z is basically distance from MC 0 to ND location
 // gdk2nu->User2BeamPos(X4det, X4beam); //X4beam will be in cm. Need to do rotations etc from Detector to Beam Frame of Reference 
  Xbeam = X4det.Vect();  // only want spatial part
  Xbeam =  RotateToBeam (Xbeam);
  bsim::calcEnuWgt(dk2nu, Xbeam, Enu, wgt_xy);
}

//=============================================================================
// GetPOTfromMeta
//=============================================================================
int GetPOTfromMeta(std::string fluxpath){
//this is modified because we will only read one rootfile per job....
  int pot   = 0;
  int added = 0;
 TChain *cmeta = new TChain("dkmetaTree");
 cmeta->Add(fluxpath.c_str());
 bsim::DkMeta* dkmeta = new bsim::DkMeta;
 cmeta->SetBranchAddress("dkmeta",&dkmeta);
 Long64_t meta_entries = cmeta->GetEntries();
 for(Long64_t jentry = 0;jentry<meta_entries;jentry++){
 cmeta->GetEntry(jentry);
  pot += dkmeta->pots;
 }
 delete cmeta;
 delete dkmeta;
 return pot;
 }

bool InsideDetector(int x, int y){
    //basically the origin should be the center of the detector.
    //You can add the complications related to the DUNE ND later on...I will assume the following:

    if(abs(x)>dunexdim && abs(y)>duneydim )return false;
    return true; 
}
