#include <time.h>
#include <stdio.h>
#include <sstream>

#include "extract_duneflux.h"
#include "CommonIncludes.h"

#include "TH1.h"
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/calcLocationWeights.h"
#include "dk2nu/tree/dkmeta.h"

void LoopEntries( TChain* cflux, bool grid, bool debug, int StatsMultiplier);
void calcEnuWgt(double x, double y, double z, bsim::Dk2Nu* dk2nu, double& Enu, double& wgt_xy);
int GetPOTfromMeta(const char *dirname, const char *extension);

  bool do_1D_hists = true;
  bool do_2D_hists= true;
//arguments required in original jobsub python script might be ....
//shift_name, shift_sigma and Statsmultiplier and a boolean to know if we want to run the extract flux command or not
//i guess. only 1D histogram for now....

int dune2Dflux( int StatsMultiplier = 1,std::string influx="*.root")
{

  bool grid = false;
  bool debug = true;


  //outdirectory might not need different naming convention....
  time_t timer;
  std::ostringstream ss;
  ss<<timer;
  std::string time_stamp = ss.str();
  
  std::cout<<" Making Output files with directory "<<std::endl;
   std::string fout_dir  = GetOutDir(grid, debug);
    std::string fout_name = Form("outfile_fluxes_%s.root", 
                                     time_stamp.c_str() );
    
    if(debug) fout_name = "test.root";
    fout_name = fout_dir+fout_name;
    TFile *fout = new TFile( fout_name.c_str(), "RECREATE"); 

   //============================================================================
  // Initialize Hists
  //============================================================================
    TH1::SetDefaultSumw2();
    InitializeHists(fout);

  std::cout<<"Searching ntuples in "<<influx<<std::endl;    
  
  TChain *cflux = new TChain("dk2nuTree");
  cflux->Add(influx.c_str());
  
  std::cout<<" Added files to TChain "<<influx<<std::endl;
    
  //============================================================================
  // Get POT info from meta
  //============================================================================
    double pot = GetPOTfromMeta(influx);
    hPOT->Fill(0.5,pot);

    std::cout << "\n  " << pot << " POT found" << std::endl;

    //fill a stats multiplier branch
    hStatsMultiplier->Fill((double)(StatsMultiplier+0.5));

  //============================================================================
  // Loop Entries
  //============================================================================
    LoopEntries ( cflux, grid, debug, StatsMultiplier);

  fout->Write();

  fout->Close();
  
  std::cout << "\n  flux file made: " << fout_name <<"\n"<< std::endl;

  return 0;
}

//==============================================================================
// Loop Entries
//==============================================================================
void LoopEntries( TChain* cflux,bool grid, bool debug, int StatsMultiplier)
{
  std::cout << "\n  Beginning to loop over ntuples and fill histos "<<std::endl;


  //============================================================================
  // Load trees
  //============================================================================
    bsim::Dk2Nu*  dk2nu  = new bsim::Dk2Nu;
    cflux->SetBranchAddress("dk2nu", &dk2nu);

  //============================================================================
  // The actual loop
  //============================================================================
    Long64_t nentries = cflux->GetEntries();
    std::cout << "\n  nentries = " << nentries << std::endl;
    nentries = debug ? 25000 : nentries;
    std::cout << "\n  Looping over " << nentries << " entries." << std::endl;
    
    for (Long64_t i=0; i < nentries; ++i ) {
      if (i % (nentries/10000) == 0) std::cout << "  entry " << i << std::endl;

      cflux->GetEntry(i);

      //neutrino pid
      int decayntype = 0;
      decayntype  = dk2nu->decay.ntype;
      
      int nutype = -1;
      if (decayntype == 14)       nutype = 0;
      else if (decayntype == -14) nutype = 1;
      else if (decayntype == 12)  nutype = 2;
      else if (decayntype == -12) nutype = 3;

      //default g4numi energy and weights
      //for ntuples processed with g4numi before 03/2017, the position weights will
      //be NOT for the center of minerva (-0.24, -0.24, 1303.) but for a different position
      //(-0.55, -0.55, 1302.) which is wrong.
      //check the dk2nu/etc/locations.txt to make sure you're getting the right position.
      double decaynimpwt=-999.;
      decaynimpwt = dk2nu->decay.nimpwt;  //importance weighting


        double Enu, wgt_xy;

        //center of the Front face of the DUNE ND
        calcEnuWgt(0, 0, 0.0, dk2nu, Enu, wgt_xy);
        hFluxCenter[nutype]->Fill(Enu, wgt_xy*decaynimpwt);
      //  std::cout<<Enu<<" "<<wgt_xy<<" "<<dk2nu->nuray[1].E<<" "<<dk2nu->nuray[1].wgt<<std::endl;
      //========================================================================
      // 1DHists -- flux at rndm x-y-z point in MINERvA tracker
      //         -- single x-y bin AND daisy flower binning
      //========================================================================
      if(do_1D_hists){
        for( int N = 0; N < StatsMultiplier; ++N){
          double xin = 0.0;
          double yin = 0.0;
          //first pick a random point inside of a minerva hex, centered at 0,0
           xin = generator.Uniform(-dunexdim,dunexdim);
           yin = generator.Uniform(-duneydim,duneydim);


          // random z position
          double zin = generator.Uniform(-dunezdim, dunezdim);

          // single x-y bin
          calcEnuWgt(xin, yin, zin, dk2nu, Enu, wgt_xy);
         // std::cout<<xin<<" "<<yin<<" "<<zin<<" "<<wgt_xy<<std::endl;
          hFluxRndDetXYZ[nutype]->Fill(Enu, wgt_xy*decaynimpwt);


        }//end stats multiplier
      }//end 1d hists

      //========================================================================
      //2DHists -- single z position
      //           downstream, center, or upstream.
      //========================================================================
      if(do_2D_hists){
        //perform the loops in (cm)...
        for (int j = loopxmin+1; j <= loopxmax; j = j+loopxbinwidth){
          for (int k = loopymin+1; k <= loopymax; k = k+loopybinwidth){
            double x = double(j); //(cm)
            double y = double(k); //(cm)

            //...but provide calcEnuWgt with (m)
            x /= 100.;     //(m)
            y /= 100.;     //(m)
            //std::cout << x << ", " << y << std::endl;

            // random z position
            double zin = generator.Uniform(-dunezdim/2, dunezdim/2);
            
            //Again: provide this function with positions in (m)
            calcEnuWgt(x, y, zin, dk2nu, Enu, wgt_xy);
            hXYFluxRndZ[nutype]->Fill(x, y, wgt_xy*decaynimpwt);
            // fill histos with (m)
            if (Enu<4)       hXYFlux0_4GeVRndZ[nutype]->Fill(x,y,wgt_xy*decaynimpwt);
            else if (Enu<5)  hXYFlux4_5GeVRndZ[nutype]->Fill(x,y,wgt_xy*decaynimpwt);
            else if (Enu<6)  hXYFlux5_6GeVRndZ[nutype]->Fill(x,y,wgt_xy*decaynimpwt);
            else if (Enu<7)  hXYFlux6_7GeVRndZ[nutype]->Fill(x,y,wgt_xy*decaynimpwt);
            else if (Enu<8)  hXYFlux7_8GeVRndZ[nutype]->Fill(x,y,wgt_xy*decaynimpwt);
            else if (Enu<10) hXYFlux8_10GeVRndZ[nutype]->Fill(x,y,wgt_xy*decaynimpwt);
            else if (Enu<12) hXYFlux10_12GeVRndZ[nutype]->Fill(x,y,wgt_xy*decaynimpwt);
            else if (Enu<14) hXYFlux12_14GeVRndZ[nutype]->Fill(x,y,wgt_xy*decaynimpwt);
            else if (Enu<20) hXYFlux14_20GeVRndZ[nutype]->Fill(x,y,wgt_xy*decaynimpwt);
            else             hXYFlux20_50GeVRndZ[nutype]->Fill(x,y,wgt_xy*decaynimpwt);
          }
        }//end x-y loop
      }//end 2d hists
   
    }//end entry loop
  
  std::cout << "\n  Done filling histos for "<< std::endl;
}

int main(int argc, char* argv[])
{
  int stat_multiplier = atoi(argv[1]);
  std::string influx = argv[2];
dune2Dflux( stat_multiplier, influx);
}
