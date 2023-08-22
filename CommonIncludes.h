#ifndef CommonIncludes_h
#define CommonIncludes_h
#include <string>
#include <iostream>
#include <map>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TSystem.h"

#include "TLegend.h"
#include "TLegendEntry.h"




//==================================================================
//WHERE TO SAVE THE OUTPUT FILES
//==================================================================

  //const std::string TOPDIR_PNFS    = "/pnfs/minerva/persistent/users/bmesserl/flux/ME_uncerts/";
  std::string TOPDIR_PNFS = "/pnfs/dune/persistent/users/bashyal8/grid_flux_test/";
  const std::string TOPDIR_BLUEARC = "/dune/data/users/bashyal8/flux/dune2D/";


  const std::string INDIR_GRID     = "$CONDOR_DIR_INPUT/";
  const std::string OUTDIR_GRID    = "$CONDOR_DIR_EXTRACTFLUX/";
  const std::string FLUXFILE_GRID  = INDIR_GRID + "/sample_dk2nu.root";

  //INDIR_NOGRID depends on specifics. See GetInDir().
  const std::string OUTDIR_NOGRID       = TOPDIR_BLUEARC;
  const std::string OUTDIR_NOGRID_DEBUG = OUTDIR_NOGRID + "test/";

  TRandom3 generator(0);
//=============================================================================
// Declare Histograms
//=============================================================================
TH1D* hPOT;
TH1D* hStatsMultiplier;

//1D hists in particular xyz positions
TH1F* hFluxCenter[4];

//1D hists integrated over det
TH1F* hFluxRndDetXYZ[4];

//2D hists at various z positions, in various bins of Enu
TH2F* hXYFluxRndZ[4];
TH2F* hXYFlux0_4GeVRndZ[4],  *hXYFlux4_5GeVRndZ[4],   *hXYFlux5_6GeVRndZ[4],   *hXYFlux6_7GeVRndZ[4],   *hXYFlux7_8GeVRndZ[4];
TH2F* hXYFlux8_10GeVRndZ[4], *hXYFlux10_12GeVRndZ[4], *hXYFlux12_14GeVRndZ[4], *hXYFlux14_20GeVRndZ[4], *hXYFlux20_50GeVRndZ[4];

 static const Double_t bins[] = { 0, .5,  1,  1.5, 2, 2.5, 3, 3.5,  4,  4.5, 5, 5.5,
                                       6,  6.5,  7,  7.5,  8,  8.5,  9,  9.5, 10, 10.5,
                                      11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5,
                                      16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5,
                                      21, 21.5, 22, 22.5, 23, 23.5, 24, 24.5, 25, 25.5,
                                      26, 26.5, 27, 27.5, 28, 28.5, 29, 29.5, 30.,
                                      31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 
                                      41., 42., 43., 44., 45., 46., 47., 48., 49., 50.};

  const Int_t nbins = sizeof(bins)/sizeof(bins[0]) - 1;

  //Binning for 2D histograms
  const double xmin = -4.; //m
  const double xmax =  4.; //m
  const int xbins = 32;
  const double xbinwidth = (xmax - xmin)/double(xbins); //0.05m

  const double ymin = -4; //m
  const double ymax =  4; //m
  const int ybins =  32;
  const double ybinwidth = (ymax - ymin)/double(ybins); //0.05m

  //Loop position bins in cm, but fill histos in m
  const int loopxmin = int(xmin*100);         //-400cm
  const int loopxmax = int(xmax*100);         // 400cm
  const int loopymin = int(ymin*100);         //-400cm
  const int loopymax = int(ymax*100);         // 400cm
  const double loopxbinwidth = xbinwidth*100; //   5cm
  const double loopybinwidth = ybinwidth*100; //   5cm


void InitializeHists(TFile* f)
{
  std::cout<<"Initializing Histograms "<<std::endl;
      f->mkdir("numu");
      f->mkdir("numu/1D");
      f->mkdir("numu/2D");

      f->mkdir("numubar");
      f->mkdir("numubar/1D");
      f->mkdir("numubar/2D");

      f->mkdir("nue");
      f->mkdir("nue/1D");
      f->mkdir("nue/2D");

      f->mkdir("nuebar");
      f->mkdir("nuebar/1D");
      f->mkdir("nuebar/2D");

  hPOT = new TH1D("hPOT", "POT", 1, 0, 1);
  hStatsMultiplier = new TH1D("hStatsMultiplier", "StatsMultiplier", 2000, 0, 2000);

  std::vector<std::string> nus;
  nus.push_back("numu"); nus.push_back("numubar");
  nus.push_back("nue");  nus.push_back("nuebar");
  std::vector<std::string>::iterator it;
  for (it = nus.begin(); it != nus.end(); ++it){
    std::string ns = *it;
    //std::cout << ns << std::endl;
    f->cd(ns.c_str());

    int index = -1;
    if (ns == "numu")       index = 0;
    else if (ns == "numubar") index = 1;
    else if (ns == "nue")     index = 2;
    else if (ns == "nuebar")  index = 3;

    static const std::string positions[] = {"Center", "WrongCenter", "Upstream", "Downstream", "Top", "Bottom", "Left", "Right"};
    std::vector<std::string> positions_vec (positions, positions + sizeof(positions) / sizeof(positions[0]) );

    f->cd(Form("%s/1D", ns.c_str()));

    //1D
    //single point in the center of the detector
    hFluxCenter[index]           = new TH1F("hFluxCenter",      Form("%s flux -- Detector center",ns.c_str()),    nbins, bins);

    //integrated over whole detector
    hFluxRndDetXYZ[index]        = new TH1F("hFluxRndDetXYZ",   Form("%s flux -- Random det XYZ",ns.c_str()),     nbins, bins);

    f->cd(Form("%s/2D", ns.c_str()));

    //2d Random Z
    // integrated over all energy bins
    hXYFluxRndZ[index]         = new TH2F("hXYFluxRndZ",       Form("%s flux all GeV, RndZ ",    ns.c_str()),  xbins ,  xmin,  xmax , ybins , ymin , ymax);
  
    // energy bins
    hXYFlux0_4GeVRndZ[index]   = new TH2F("hXYFlux0_4GeVRndZ", Form("%s flux 0-4GeV, RndZ ",     ns.c_str()),  xbins ,  xmin,  xmax , ybins , ymin , ymax);
    hXYFlux4_5GeVRndZ[index]   = new TH2F("hXYFlux4_5GeVRndZ", Form("%s flux 4-5GeV, RndZ ",     ns.c_str()),  xbins ,  xmin,  xmax , ybins , ymin , ymax);
    hXYFlux5_6GeVRndZ[index]   = new TH2F("hXYFlux5_6GeVRndZ", Form("%s flux 5-6GeV, RndZ ",     ns.c_str()),  xbins ,  xmin,  xmax , ybins , ymin , ymax);
    hXYFlux6_7GeVRndZ[index]   = new TH2F("hXYFlux6_7GeVRndZ", Form("%s flux 6-7GeV, RndZ ",     ns.c_str()),  xbins ,  xmin,  xmax , ybins , ymin , ymax);
    hXYFlux7_8GeVRndZ[index]   = new TH2F("hXYFlux7_8GeVRndZ", Form("%s flux 7-8GeV, RndZ ",     ns.c_str()),  xbins ,  xmin,  xmax , ybins , ymin , ymax);
    hXYFlux8_10GeVRndZ[index]  = new TH2F("hXYFlux8_10GeVRndZ", Form("%s flux 8-10GeV, RndZ ",   ns.c_str()),  xbins ,  xmin,  xmax , ybins , ymin , ymax);
    hXYFlux10_12GeVRndZ[index] = new TH2F("hXYFlux10_12GeVRndZ", Form("%s flux 10-12GeV, RndZ ", ns.c_str()),  xbins ,  xmin,  xmax , ybins , ymin , ymax);
    hXYFlux12_14GeVRndZ[index] = new TH2F("hXYFlux12_14GeVRndZ", Form("%s flux 12-14GeV, RndZ ", ns.c_str()),  xbins ,  xmin,  xmax , ybins , ymin , ymax);
    hXYFlux14_20GeVRndZ[index] = new TH2F("hXYFlux14_20GeVRndZ", Form("%s flux 14-20GeV, RndZ ", ns.c_str()),  xbins ,  xmin,  xmax , ybins , ymin , ymax);
    hXYFlux20_50GeVRndZ[index] = new TH2F("hXYFlux20_50GeVRndZ", Form("%s flux 20-50GeV, RndZ ", ns.c_str()),  xbins ,  xmin,  xmax , ybins , ymin , ymax);

  }
}

std::string GetOutDir(bool grid, bool debug)
{
  std::string fout_dir; 
  if(grid) fout_dir = OUTDIR_GRID;
  else{
    if(debug) fout_dir = OUTDIR_NOGRID;
    else fout_dir      = OUTDIR_NOGRID;
  }
  std::cout<<"CommonIncludes.h "<<fout_dir<<std::endl;
  return gSystem->ExpandPathName(fout_dir.c_str());
}

#endif
