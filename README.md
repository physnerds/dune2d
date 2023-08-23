## Purpose of the Code
1. This code attempts to estimate the neutrino flux at the DUNE ND assuming the detector has physical dimensions throughout a physical detector volume (X(left,right with respect to the beam co-ordinates), Y (top,bottom with respect to beam co-ordinates),Z(front, back with respect to beam co-ordinates)).
  
2. This test suite was built on dunegpvm machine.
3. Code is not committed in dune repo since it needs to go through extensive tests. But locally, it lives here: /dune/app/users/bashyal8/g4lbne_dev/v3r5p9/g4lbne/BeamSimStudies/dune2d/
4. To build the code, follow the below instructions:
5. 
    a. $mkdir wordir

    b. $cd workdir
   
    c. In your dunegpvm area, $git clone https://github.com/physnerds/dune2d
   
    d. Edit output directory:
      ``std::string TOPDIR_BLUEARC = "/dune/data/users/bashyal8/flux/dune2D/";`` in `CommonIncludes.h`
      Note that the code is committed with `debug` flag on. To get this thing running in grid, you need to turn off the `debug` flag and edit the `pnfs` area accordingly. I have not tested the grid mode of this code yet.
   
      Make sure that "${USER}/flux/dune2D" exists. ${USER} is your username. Otherwise do it on your own....
   
    e. $source set_env.sh
   
    f. $mkdir build; cd build
    
    g. $cmake ..
   
    h. If the build is succesful you can try: `$./dune2Dflux 1 /pnfs/dune/persistent/users/djena/DUNE2023/20230728TgtNeutrino/v3r5p10/QGSP_BERT/OfficialEngDesignSept2021/neutrino/flux/g4lbne_v3r5p10_QGSP_BERT_OfficialEngDesignSept2021_neutrino_00028.dk2nu.root`
   
    i. Output will appear in TOPDIR_PNFS directory.
   
    j. `1` in h is the hstat and you can change it to any number you want to reuse the pion decay info. For example `1` means a pion decay info is used only once to get neutrino flux and energy. `2` would mean, it is used two times to get the flux and neutrino flux and so on. The multiplier is based on a uniform random number generator across the geometry of the Near Detector. The multiplier info is saved in the histogram `hStatMultiplier`
   

## TESTS DONE IN DUNEGPVM
 The main function to calculate weight is calcEnuWgt in extract_duneflux. The function takes x,y,z (with respect to the g4lbnf MC 0) co-ordinates to calculate the flux at the DUNE detector location. The relies on dk2nu function. dk2nu also pre-calculates the flux at the center of the front fact of the detector location in the ntuple dk2nu->nuray[1].wgt and dk2nu->nuray[1].Enu. I verified that the flux over fiducial volume is implemented correctly by forcing all the neutrinos to project over the center of the detector and then comparing against the above variables. 
 Besides this, I think a proper calculation of flux will need the rotation of co-ordinates from beam frame (in which the neutrino info is saved in the dk2nu file) to the detector frame (where we want to find the detector). 

 ## RotateToBeam(det_loc) 
 This function is meant to rotate the detector location from detector frame of reference to beam frame of reference. Right now the rotation angles theta and phi are set at 0 (i.e no rotation required) but could update the values. I have not done the cross-check with the function and hope that the users do these cross-checks before showing these results in wider audience.

### Angles in RotateToBeam(det_loc)
The angles are set in "etract_duneflux.h".
The function that actually does the rotation is "RotateToBeam(TVector3 theta_rot)" which then calls the functions "RotatePhi and RotateTheta" as consequitive linear operation. Both Theta and Phi are meant to be in radians (AND NOT DEGREES) and assumes that the axis of rotation between the Beam Frame and Detector Frame is MC 0. 



