#ifndef extract_duneflux_h
#define extract_duneflux_h

#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector.h"
//#include "CommonIncludes.h"

#include "TSystemDirectory.h"
#include "TChain.h"

#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/calcLocationWeights.h"
#include "dk2nu/tree/dkmeta.h"


//==================================================================
//AB: DUNE Near Detector Dimensions
//==================================================================
    //basically the origin should be the center of the detector.
    //You can add the complications related to the DUNE ND later on...I will assume the following:
    //Based on Info provided by Deepika Jena (except for Z)
    //length along X (left-right) = 6m/2
    //length along Y (top-bottom) = 2m/2
    //I will just assume the detector is 5 m along Z
    constexpr double dunexdim = 3.0; 
    constexpr double duneydim = 1.0;
    constexpr double dunezdim = 2.5;
//AB: Angle of rotations from Detector To Beam Frame of Reference in radians....
//Note that the beam is slanted wrto the detector....
    constexpr double theta_det_to_beam = 0.0; //3 degrees...convert to radians
    constexpr double phi_det_to_beam = 0.0;

//AB: Rotations to go from Detector Frame of Reference To Beam Frame Of Reference
//AB: Phi to rotate along Z and Theta to Rotate Along Y axis....

TVector3 RotatePhi(TVector3 det_loc); // AB: Give the co-ordinates in detector frame of reference
TVector3 RotateTheta(TVector3 det_loc); //AB: Give the co-ordinates in detector frame of reference
TVector3 RotateToBeam(TVector3 thet_rot); //AB: Give the co-ordinates in detector frame of reference

bool InsideDetector(int x, int y);

void calcEnuWgt(double x, double y, double z, 
                bsim::Dk2Nu* dk2nu, 
                double& Enu, 
                double& wgt_xy);

//=============================================================================
// GetPOTfromMeta
//=============================================================================
int GetPOTfromMeta(std::string fluxpath);
  
#endif //ifndef extract_duneflux_h
