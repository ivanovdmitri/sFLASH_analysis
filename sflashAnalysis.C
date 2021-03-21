#include "TObject.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TF1.h"
#include "TProfile.h"
#include "TLine.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <iostream>
#include <cstdio>
#include <cstring>
#include <map>

const Int_t NT0_MAX   = 10000;   // maximum number of time slices in any waveform readout
const Int_t N_MAX_OBJ = 10000;   // maximum number of objects for plotting
TCanvas *c1 = 0;                 // canvas for plotting


// PMT calibration formula for the TA photomultiplier tubes.
// Default values are the CRAYS calibration values.
// The formula as well as the CRAYS calibration values have been proved by B. K. Shin
// in an e-mail to D. Ivanov on September 19, 2019.
Double_t get_pmt_nph_per_nVs_CRAYS(Double_t hv, Double_t offset=267, Double_t alpha=-36.7009, Double_t beta=4.78707)
{
  return 1.0/(TMath::Exp(alpha) * TMath::Power(hv-offset,beta));
}



class Pulse
{
public:
  Double_t csignal;    // calibrated signal (signal rescaled by the appropriate calibration factor)
  Double_t d_csignal;  // calibrated signal uncertainty (signal uncertainty rescaled by the appropriate calibration factor)
  Double_t signal;     // pulse area minus the pedestal
  Double_t d_signal;   // uncertainty on the signal
  Double_t area;       // area of the pulse
  Double_t ped_mean_bp;   // mean pedestal per time slice before pulse
  Double_t ped_std_bp;    // standard deviation of the pedestal before pulse

  Double_t ped_mean_ap;   // mean pedestal per time slice after pulse
  Double_t ped_std_ap;    // standard deviation of the pedestal after pulse

  Int_t    width;      // width of the pulse in time slices
  Int_t    nt_ped_bp;  // number of first time slices used in determining the pedestal before pulse
  Int_t    nt_ped_ap;  // number of first time slices used in determining the pedestal after pulse
  Int_t    istart;     // signal start time slice
  Int_t    imax;       // time slice of the signal maximum
  Int_t    wf_i0;      // waveform start index
  Int_t    wf_nt;      // numbe of time slices in the entire waveform
  Pulse () { clear(); }
  virtual ~Pulse() { ; }
  void clear()
  {
    csignal          = 0;
    d_csignal        = 0;
    signal           = 0;
    d_signal         = 0;
    area             = 0;
    ped_mean_bp      = 0;
    ped_std_bp       = 0;
    ped_mean_ap      = 0;
    ped_std_ap       = 0;
    width            = 0;
    nt_ped_bp        = 0;
    nt_ped_ap        = 0;
    istart           = 0;
    imax             = 0;
    wf_nt            = 0;
    
  }
  Bool_t get(Int_t nt, const Double_t* wf)
  {
    // clear out the previous information
    clear();
    wf_i0 = 0;     // waveform start index
    wf_nt = nt;    // waveform number of time slices
    // Estimate pedestal with mean of first 100 bins in waveform
    Int_t n_ped_0 = nt / 5; // number of time slices (from start) for preliminary pedestal estimation
    const Double_t signal_stop_fract = 1e-2;  // if signal less than this fraction of the maximum, 
    // stop the signal window
    const Int_t signal_extend = 10;           // extend signal window to +/- this value
    // from the window determined by considering
    // signal behaviour
    Double_t ped_0 = 0.0;   // pedestal of the first n_ped_0 time bins
    for(Int_t i = 0; i < nt; i++)
      {
	if(i < n_ped_0)
	  ped_0 += wf[i];
	if(wf[i] < wf[imax])
	  imax = i;
      }
    ped_0 /= (Double_t) n_ped_0;
    // Find rising and falling edges of pulse
    for (istart = imax; istart >=0 && wf[istart] < ped_0 ; istart--);
    for (width  = imax - istart + 1;  istart+width <= nt && wf[istart+width-1] < signal_stop_fract*wf[imax]; width++);
    // Adjust the start and end points for the signal
    istart -= signal_extend;
    width  += 2 * signal_extend;
    if(istart < 0)
      istart = 0;
    if(istart + width > nt)
      width = nt - istart;
    // Re-compute pedestal (up to 10 bins before rising edge); pedestal before pulse
    ped_mean_bp = 0; // mean pedestal per time slice
    ped_std_bp  = 0; // standard deviation of the pedestal per time slice
    nt_ped_bp = (istart > 10 ? istart : 100);
    for (Int_t i = 0; i < nt_ped_bp; i++)
      {
	ped_mean_bp     += wf[i];
	ped_std_bp      += wf[i]*wf[i];
      }
    ped_mean_bp /= (Double_t) nt_ped_bp;
    ped_std_bp  /= (Double_t) nt_ped_bp;
    ped_std_bp  -= ped_mean_bp * ped_mean_bp;
    ped_std_bp  = TMath::Sqrt(ped_std_bp);
    ped_mean_ap = 0; // mean pedestal per time slice after pulse
    ped_std_ap  = 0; // standard deviation of the pedestal per time slice after pulse
    // get the width of the window for pedestal after pulse calculation
    nt_ped_ap = (istart + width < nt - 10 ? nt - width - istart : 10);
    for (Int_t i = nt - nt_ped_ap; i < nt; i++)
      {
	ped_mean_ap     += wf[i];
	ped_std_ap      += wf[i]*wf[i];
      }
    ped_mean_ap /= (Double_t) nt_ped_ap;
    ped_std_ap  /= (Double_t) nt_ped_ap;
    ped_std_ap  -= ped_mean_ap * ped_mean_ap;
    ped_std_ap  = TMath::Sqrt(ped_std_ap);
    // calculate the pulse area
    area = 0;
    double eff_ped = 0.5 * (ped_mean_bp + ped_mean_ap);
    for (Int_t i=istart; i < istart + width; i++)
      area += wf[i];
    // subtract the pedestal from the area to get the signal and convert the result to a positive value
    signal = TMath::Abs(area - eff_ped * ((Double_t)width));
    // determine the uncertainty on signal = area - pedestal
    d_signal = TMath::Sqrt(((Double_t)width) * ped_std_bp * ped_std_bp * (1.0 + 1.0 / (Double_t)nt_ped_bp));
    // return success
    return true;
  }
  
  void calibrate(Double_t calibration_factor)
  { 
    csignal    = calibration_factor * signal;
    d_csignal  = calibration_factor * d_signal;
  }
  
};

Int_t sflashAnalysis(const char* infile, const char* outfile)
{

  // open the input ROOT tree file and locate the raw data
  // ROOT tree in the file
  TFile *ifl = new TFile(infile,"r");
  if(ifl->IsZombie())
    return 2;
  TTree* tsFLASHwf = (TTree *)ifl->Get("tsFLASHwf");
  if(!tsFLASHwf)
    {
      fprintf(stderr,"error: failed to get tsFLASHwf from %s!\n", ifl->GetTitle());
      return 2;
    }

  // Variables from the raw waveform tree that are required for analyzing the coild and PMTs.
  Int_t nt0; // number of time slices for scope 0
  Double_t wf_coil[NT0_MAX];  // coil waveform
  Double_t wf_pmt_1[NT0_MAX]; // PMT_1 waveform
  Double_t C_per_Vs;          // coil calibration contant
  Double_t PMT_1_HV;          // PMT_1 high voltage setting [ V ]
  Double_t PMT_1_OFFSET;      // PMT_1 calibration information (presently from CRAYS)
  Double_t PMT_1_ALPHA;
  Double_t PMT_1_BETA;
  Int_t    run_type;          // 1 - Beam, 2 = UVLED
  Int_t    status;            // 1 - good, 0 = not
  Int_t    blind;             // 1 - blind open, 0 - blind closed
  Int_t    shutter;           // 1 - shutter open, 0 - shutter closed

  // point the tree branches to the variables that we need
  tsFLASHwf->SetBranchAddress("run_type",    &run_type);
  tsFLASHwf->SetBranchAddress("status",      &status);
  tsFLASHwf->SetBranchAddress("blind",       &blind);
  tsFLASHwf->SetBranchAddress("shutter",     &shutter);
  tsFLASHwf->SetBranchAddress("nt0",         &nt0);
  tsFLASHwf->SetBranchAddress("wf_coil",     wf_coil);
  tsFLASHwf->SetBranchAddress("wf_pmt_1",    wf_pmt_1);
  tsFLASHwf->SetBranchAddress("C_per_Vs",    &C_per_Vs);
  tsFLASHwf->SetBranchAddress("PMT_1_HV",    &PMT_1_HV); 
  tsFLASHwf->SetBranchAddress("PMT_1_OFFSET",&PMT_1_OFFSET);
  tsFLASHwf->SetBranchAddress("PMT_1_ALPHA", &PMT_1_ALPHA);
  tsFLASHwf->SetBranchAddress("PMT_1_BETA",  &PMT_1_BETA);

  
   
  // Start a new ROOT tree that is a copy of the original tree
  // plus additional branches that are the results of the analysis
  // in this routine. Link the new ROOT tree to a new (output)
  // ROOT file.
  Pulse*    coil  = new Pulse();      // singal analysis results of the coil
  Pulse*    pmt_1 = new Pulse();      // signal analysis of the pmt_1
  Int_t     beam_open_open      = 0;  // 1 = good open beam run, 0 not
  Int_t     beam_closed_open    = 0;  // 1 = good beam run, blind closed, shutter open, 0 not
  Int_t     beam_closed_closed  = 0;  // 1 = blind closed shutter closed, 0 not
  // starting a new ROOT output file
  TFile *ofl = new TFile(outfile,"recreate");
  if(ofl->IsZombie())
    return 2;
  // starting a new ROOT tree
  TTree *tsFLASHwfAnalyzed =  tsFLASHwf->CloneTree(0);
  // additional branches for the output ROOT tree
  tsFLASHwfAnalyzed->Branch("coil",                "Pulse",              &coil);
  tsFLASHwfAnalyzed->Branch("pmt_1",               "Pulse",              &pmt_1);
  tsFLASHwfAnalyzed->Branch("beam_open_open",      &beam_open_open,      "beam_open_open/I");
  tsFLASHwfAnalyzed->Branch("beam_closed_open",    &beam_closed_open,    "beam_closed_open/I");
  tsFLASHwfAnalyzed->Branch("beam_closed_closed",  &beam_closed_closed,  "beam_closed_closed/I");
  
  // analysis loop: going over all events in the original ROOT tree, analyzing the events,
  // and filling the original data + analysis results to a new ROOT tree.
  for (Int_t ievent =0; ievent < tsFLASHwf->GetEntries(); ievent++)
    {
      // load data from the tree that's being analyzed
      tsFLASHwf->GetEntry(ievent);

      
      // analyze the coil waveform
      coil->get(nt0,wf_coil);
      
      // determine beam charge in pC units from the coil signal and the calibration constant
      coil->calibrate(1e12 * C_per_Vs * 1e-9);
      
      // analyze the PMT_1 waveform
      pmt_1->get(nt0,wf_pmt_1);
      
      // determine the number of photons for the PMT_1 from the PMT_1 signal and the CRAYS calibration
      pmt_1->calibrate(get_pmt_nph_per_nVs_CRAYS(PMT_1_HV,PMT_1_OFFSET,PMT_1_ALPHA,PMT_1_BETA));
      
      
      
      // blind open shutter open ? 
      beam_open_open = 1; // initialize
      beam_open_open = (run_type == 1 ? beam_open_open : 0);                    // beam run
      beam_open_open = (status   == 1 ? beam_open_open : 0);                    // good run
      beam_open_open = (blind    == 1 ? beam_open_open : 0);                    // blind open
      beam_open_open = (shutter  == 1 ? beam_open_open : 0);                    // shutter open
      
      // blind closed door closed ?
      beam_closed_open = 1; // initialize
      beam_closed_open = (run_type == 1 ? beam_closed_open : 0);    // beam run
      beam_closed_open = (status   == 1 ? beam_closed_open : 0);    // good run
      beam_closed_open = (blind    == 0 ? beam_closed_open : 0);    // blind closed
      beam_closed_open = (shutter  == 1 ? beam_closed_open : 0);    // shutter open

      // blind open door open ?
      beam_closed_closed = 1; // initialize
      beam_closed_closed = (run_type == 1 ? beam_closed_closed : 0);  // beam run
      beam_closed_closed = (status   == 1 ? beam_closed_closed : 0);  // good run
      beam_closed_closed = (blind    == 0 ? beam_closed_closed : 0);  // blind closed
      beam_closed_closed = (shutter  == 0 ? beam_closed_closed : 0);  // shutter closed
      
      // Fill the tree with the data of the previous tree plus newly added analysis branches
      tsFLASHwfAnalyzed->Fill();
      
    }
  // Write out the new tree that has the original data plus the analysis results
  tsFLASHwfAnalyzed->Write();

  // close the input and output files
  ofl->Close();
  ifl->Close();
  
  cout << "\nDone\n";
  return 0;
}

void init_canv()
{
  gStyle->SetOptFit();
  if(!c1)
    {
      c1 = new TCanvas("c1","c1",900,500);
      c1->SetLeftMargin(0.125);
      c1->SetBottomMargin(0.125);
    }
  c1->SetGridx();
  c1->SetGridy();
  c1->SetTickx();
  c1->SetTicky();
  c1->cd();  
}



TGraph *get_p_vs_q(TTree *tsFLASHwf, 
			 const char* what="pmt_1", 
			 const char* cuts = "beam_open_open==1&&rl==5")
{  
  TString s_what = what;
  s_what.ToLower();
  if(!tsFLASHwf->GetBranch(s_what))
    {
      fprintf(stderr,"error: branch '%s' is absent in tree '%s'; analyze the tree for '%s' first!\n",
	      s_what.Data(), tsFLASHwf->GetName(),s_what.Data());
      return 0;
    }
  if(!tsFLASHwf->GetBranch("coil"))
    {
      fprintf(stderr,"error: branch 'coil' is absent in tree '%s'; analyze the tree first!\n",
	      tsFLASHwf->GetName());
      return 0;
    }  
  TString var_exp = "";
  TString var_cut = "";
  var_exp.Form("coil.csignal:%s.csignal",s_what.Data());
  if(strlen(cuts))
    var_cut.Form("coil.csignal/coil.d_csignal > 5.0 && %s.csignal/%s.d_csignal > 5.0 && %s",s_what.Data(),s_what.Data(),cuts);
  else
    var_cut.Form("coil.csignal/coil.d_csignal > 5.0 && %s.csignal/%s.d_csignal > 5.0",s_what.Data(),s_what.Data());
  tsFLASHwf->Draw(var_exp,var_cut,"goff");  
  TGraph *g = new TGraphErrors(tsFLASHwf->GetSelectedRows(),
			       tsFLASHwf->GetV1(),
			       tsFLASHwf->GetV2());
  g->SetMarkerStyle(6);
  g->SetLineColor(kBlack);
  g->SetMarkerColor(kBlack);
  TString s_title = "";
  TString s_unit="";
  if(s_what.Contains("pmt"))
    s_unit = "[ photons ]";
  else if (s_what.Contains("coil"))
    s_unit = "[ pC ]";
  else
    s_unit = " [ ? ]";
  s_what.ToUpper();
  s_title.Form(";COIL [ pC ];%s %s",s_what.Data(), s_unit.Data());
  g->SetTitle(s_title);
  g->GetXaxis()->SetTitleSize(0.055);
  g->GetYaxis()->SetTitleSize(0.055);
  return g;
}


class p_vs_q_fits_class
{
public:
  TGraph*  g_open_open;
  TGraph*  g_closed_open;
  TGraph*  g_closed_closed;
  Double_t slope_beam_blind;     // corrected by subtracting blind closed data
  Double_t d_slope_beam_blind;
  Double_t slope_beam_shutter;   // corrected by subtracting shutter closed data
  Double_t d_slope_beam_shutter;

  Double_t ratio_beam_blind;
  Double_t d_ratio_beam_blind;

  p_vs_q_fits_class()
  {
    slope_beam_blind     = 0;
    d_slope_beam_blind   = 0;
    slope_beam_shutter   = 0;
    d_slope_beam_shutter = 0;
    g_open_open          = 0;
    g_closed_open        = 0;
    g_closed_closed      = 0;
    ratio_beam_blind     = 0;
    d_ratio_beam_blind   = 0;
  }
  ~p_vs_q_fits_class()
  {
    if(g_open_open)
      delete g_open_open;
    if(g_closed_open)
      delete g_closed_open;
    if(g_closed_closed)
      delete g_closed_closed;
  }
  void get_slope()
  {
    TF1 *f_open_open     = g_open_open->GetFunction("pol1");
    TF1 *f_closed_open   = g_closed_open->GetFunction("pol1");
    TF1 *f_closed_closed = g_closed_closed->GetFunction("pol1");
    if(!f_open_open || !f_closed_open || !f_closed_closed)
      {
	cerr << "error: open and/or closed graphs haven't been fitted yet!" << endl; 
	return;
      }
    
    slope_beam_blind       = f_open_open->GetParameter(1) - ( f_closed_open->GetParameter(1) > 0 ? 
							      f_closed_open->GetParameter(1) :  0);
    d_slope_beam_blind     = TMath::Sqrt(f_open_open->GetParError(1)*f_open_open->GetParError(1) + 
					 f_closed_open->GetParError(1)*f_closed_open->GetParError(1));
    slope_beam_shutter     = f_open_open->GetParameter(1) - f_closed_closed->GetParameter(1);
    d_slope_beam_shutter   = TMath::Sqrt(f_open_open->GetParError(1)*f_open_open->GetParError(1) + 
					 f_closed_closed->GetParError(1)*f_closed_closed->GetParError(1));
    
    ratio_beam_blind       = f_open_open->GetParameter(1) / TMath::Abs(f_closed_open->GetParameter(1));

    d_ratio_beam_blind     = ratio_beam_blind * TMath::Sqrt(f_open_open->GetParError(1)*f_open_open->GetParError(1) /
							    f_open_open->GetParameter(1)/f_open_open->GetParameter(1)+ 
							    f_closed_open->GetParError(1)*f_closed_open->GetParError(1)/
							    f_closed_open->GetParameter(1)/f_closed_open->GetParameter(1));

  }
  void Draw(TCanvas *canv)
  {
    if(!g_open_open || !g_closed_open || !g_closed_closed)
      {
	cerr << "p_vs_q_fits_class::Draw: graphs have not been filled yet" << endl;
	return;
      }
    if(!canv)
      {
	cerr << "p_vs_q_fits_class::Draw: canvas has not been allocated yet" << endl;
	return;
      }
      g_open_open->SetLineColor(kRed);
      g_open_open->Draw("a,p");
      g_closed_open->Draw("p,sames");
      g_closed_open->SetLineColor(kBlue);
      g_closed_closed->Draw("p,sames");
      canv->Modified();
      canv->Update();
      TPaveStats *stats = (TPaveStats *)g_open_open->FindObject("stats");
      if(stats)
	{
	  stats->SetOptStat(0);
	  stats->SetOptFit(10011);
	  stats->SetX1NDC(0.14);
	  stats->SetX2NDC(0.50);
	  stats->SetY1NDC(0.80);
	  stats->SetY2NDC(0.90);
	}
      stats = (TPaveStats *)g_closed_open->FindObject("stats");
      if(stats)
	{
	  stats->SetOptStat(0);
	  stats->SetOptFit(10011);
	  stats->SetX1NDC(0.14);
	  stats->SetX2NDC(0.50);
	  stats->SetY1NDC(0.70);
	  stats->SetY2NDC(0.80);
	}
      stats = (TPaveStats *)g_closed_closed->FindObject("stats");
      if(stats)
	{
	  stats->SetOptStat(0);
	  stats->SetOptFit(10011);
	  stats->SetX1NDC(0.14);
	  stats->SetX2NDC(0.50);
	  stats->SetY1NDC(0.60);
	  stats->SetY2NDC(0.70);
	}
      canv->Modified();
      canv->Update();
  }
};


p_vs_q_fits_class* get_p_vs_q_fits(TTree *tsFLASHwf, 
				   const char* what="pmt_1",
				   Int_t rl = 5,
				   const char* other_cuts = "")
{
  TString cuts_open_open = "beam_open_open==1&&rl==";
  cuts_open_open += rl;
  if(strlen(other_cuts))
    {
      cuts_open_open += "&&";
      cuts_open_open += other_cuts;
    }
  TString cuts_closed_open = "beam_closed_open==1&&rl==";
  cuts_closed_open += rl;
  if(strlen(other_cuts))
    {
      cuts_closed_open += "&&";
      cuts_closed_open += other_cuts;
    }
  TString cuts_closed_closed = "beam_closed_closed==1&&rl==";
  cuts_closed_closed += rl;
  if(strlen(other_cuts))
    {
      cuts_closed_closed += "&&";
      cuts_closed_closed += other_cuts;
    }
  p_vs_q_fits_class *fits = new p_vs_q_fits_class;
  fits->g_open_open     = get_p_vs_q(tsFLASHwf,what,cuts_open_open);
  fits->g_closed_open   = get_p_vs_q(tsFLASHwf,what,cuts_closed_open);
  fits->g_closed_closed = get_p_vs_q(tsFLASHwf,what,cuts_closed_closed);
  fits->g_open_open->Fit("pol1","W,Q");
  fits->g_closed_open->Fit("pol1","W,Q");
  fits->g_closed_closed->Fit("pol1","W,Q");
  fits->get_slope();
  return fits;
}


void plot_p_vs_q(TTree *tsFLASHwf, 
		 const char* what="pmt_1", 
		 const char* cuts = "beam_open_open==1&&rl==5")
{
  TGraph *g = get_p_vs_q(tsFLASHwf,what,cuts);
  init_canv();
  g->Fit("pol1","W,Q");
  g->Draw("a,p");
  c1->Modified();
  c1->Update();
  TPaveStats *stats = (TPaveStats *)g->FindObject("stats");
  if(stats)
    {
      stats->SetOptStat(0);
      stats->SetOptFit(10011);
    }
  c1->Modified();
  c1->Update();
}


void plot_p_vs_q_rl(TTree *tsFLASHwf, 
		    const char* what="pmt_1",
		    Int_t rl = 5,
		    const char* other_cuts = "")
{
  p_vs_q_fits_class *fits = get_p_vs_q_fits(tsFLASHwf,what,rl,other_cuts);
  init_canv();
  TString s_title = "";
  s_title += rl;
  s_title += " R.L.";
  fits->g_open_open->SetTitle(s_title);
  fits->Draw(c1);
}


class rl_result
{
public:
  TGraphErrors* g_beam_blind;   // noise subtracted using closed blind open shutter
  TGraphErrors* g_beam_blind_to_noise_ratio; // noise divided
  TGraphErrors* g_beam_shutter; // noise subtracted using closed blind closed shutter
  TGraphErrors* g_ratio;        // g_beam_shutter / g_beam_blind
  rl_result(Int_t nrl)
  {
    g_beam_blind   = new TGraphErrors(nrl);
    g_beam_shutter = new TGraphErrors(nrl);
    g_ratio        = new TGraphErrors(nrl);
    g_beam_blind_to_noise_ratio = new TGraphErrors(nrl);
  }

  void ScaleBeamBlindByInvOf(TGraph *g)
  {
    if(!g_beam_blind)
      return;
    for (Int_t i=0; i < g_ratio->GetN(); i++)
      {
	Double_t x, y, ey, sf;
	g_beam_blind->GetPoint(i,x,y);
	ey = g_beam_blind->GetErrorY(i);
	sf = 1.0 / g->Eval(x);
	g_beam_blind->SetPoint(i,x,sf*y);
	g_beam_blind->SetPointError(i,0,sf*ey);
      }
  }

  void ScaleBeamBlindBy(TGraph *g)
  {
    if(!g_beam_blind)
      return;
    for (Int_t i=0; i < g_ratio->GetN(); i++)
      {
	Double_t x, y, ey, sf;
	g_beam_blind->GetPoint(i,x,y);
	ey = g_beam_blind->GetErrorY(i);
	sf = g->Eval(x);
	g_beam_blind->SetPoint(i,x,sf*y);
	g_beam_blind->SetPointError(i,0,sf*ey);
      }
  }


  void calc_ratio()
  {
    for (Int_t i=0; i < g_ratio->GetN(); i++)
      {
	Double_t x, y1, ey1, y2, ey2, r, er;
	g_beam_blind->GetPoint(i,x,y1);
	ey1 = g_beam_blind->GetErrorY(i);
	g_beam_shutter->GetPoint(i,x,y2);
	ey2 = g_beam_shutter->GetErrorY(i);
	r = y2 / y1;
	er = r * TMath::Sqrt(ey2*ey2/y2/y2+ey1*ey1/y1/y1);
	g_ratio->SetPoint(i,x,r);
	g_ratio->SetPointError(i,0,er);
      }
  }
  ~rl_result()
  {
    if(g_beam_blind)
      delete g_beam_blind;
    if(g_beam_shutter)
      delete g_beam_shutter;
    if(g_ratio)
      delete g_ratio;
    if(g_beam_blind_to_noise_ratio)
      delete g_beam_blind_to_noise_ratio;
  }
  void SetTitle(const char* title1="", const char* title2="")
  {

    if(!g_beam_blind || !g_beam_shutter || !g_ratio)
      return;

    g_beam_blind->SetTitle(title1);
    g_beam_blind->SetMarkerStyle(20);
    g_beam_blind->SetMarkerSize(1.0);
    g_beam_blind->GetXaxis()->SetTitleSize(0.055);
    g_beam_blind->GetYaxis()->SetTitleSize(0.055);
    g_beam_blind->SetLineColor(kBlack);
    g_beam_blind->SetMarkerColor(kBlack);
    g_beam_blind->SetLineStyle(1);

    g_beam_blind_to_noise_ratio->SetTitle(title1);
    g_beam_blind_to_noise_ratio->SetMarkerStyle(20);
    g_beam_blind_to_noise_ratio->SetMarkerSize(1.0);
    g_beam_blind_to_noise_ratio->GetXaxis()->SetTitleSize(0.055);
    g_beam_blind_to_noise_ratio->GetYaxis()->SetTitleSize(0.055);
    g_beam_blind_to_noise_ratio->SetLineColor(kBlack);
    g_beam_blind_to_noise_ratio->SetMarkerColor(kBlack);
    g_beam_blind_to_noise_ratio->SetLineStyle(1);
    
    g_beam_shutter->SetTitle(title2);
    g_beam_shutter->SetMarkerStyle(25);
    g_beam_shutter->SetMarkerSize(1.0);
    g_beam_shutter->GetXaxis()->SetTitleSize(0.055);
    g_beam_shutter->GetYaxis()->SetTitleSize(0.055);
    g_beam_shutter->SetLineColor(kBlack);
    g_beam_shutter->SetMarkerColor(kBlack);
    g_beam_shutter->SetLineStyle(2);


    g_ratio->SetTitle(title1);
    g_ratio->SetMarkerStyle(20);
    g_ratio->SetMarkerSize(1.0);
    g_ratio->GetXaxis()->SetTitleSize(0.055);
    g_ratio->GetYaxis()->SetTitleSize(0.055);
    g_ratio->SetLineColor(kRed);
    g_ratio->SetMarkerColor(kRed);
    g_ratio->SetLineStyle(1);


  }

  void Draw(const char* opt = "e1p,C,both")
  {
    init_canv();
    
    // If true, draw both beam-blind (main) and beam-shutter 
    // (of secondary importance) results.
    // Otherwise, draw only the main result.
    Bool_t draw_both = false;
    TString g_opt = opt;
    g_opt.ToLower();
    if (g_opt.Contains("both"))
      {
	draw_both = true;
	Ssiz_t i = g_opt.Index(",both");
	Ssiz_t l = (Ssiz_t) strlen(",both");
	if (i == -1)
	  {
	    i  = g_opt.Index("both");
	    l  = (Ssiz_t) strlen("both");
	  }
	g_opt.Remove(i,l);
      }
    
    TString g_opt_a = g_opt;
    g_opt_a.Prepend("a,");
    
    if(draw_both)
      {
	g_beam_shutter->Draw(g_opt_a);
	g_beam_blind->Draw(g_opt+",same");
	return;
      }
    
    g_beam_blind->Draw(g_opt_a);
    
  }

  void DrawRatio()
  {
    g_ratio->Draw("a,e1p");
  }
  
};

rl_result* rl_study(TTree *tsFLASHwf, const char* what="pmt_1",const char* other_cuts = "")
{
  const Int_t nrl = 6;
  Int_t rl[nrl] = {0,2,3,5,8,10};
  rl_result *result = new rl_result(nrl);
  result->SetTitle(";j;M_{1j} [ photons / pC ]", ";j;M_{1j} [ photons / pC ]");
  for (Int_t i=0; i<nrl; i++)
    {
      p_vs_q_fits_class *fits = get_p_vs_q_fits(tsFLASHwf,what,rl[i],other_cuts);
      result->g_beam_blind->SetPoint(i, (Double_t)rl[i], fits->slope_beam_blind);
      result->g_beam_blind->SetPointError(i, 0, fits->d_slope_beam_blind);
      result->g_beam_blind_to_noise_ratio->SetPoint(i, (Double_t)rl[i], fits->ratio_beam_blind);
      result->g_beam_blind_to_noise_ratio->SetPointError(i, 0, fits->d_ratio_beam_blind);
      result->g_beam_shutter->SetPoint(i, (Double_t)rl[i], fits->slope_beam_shutter);
      result->g_beam_shutter->SetPointError(i, 0, fits->d_slope_beam_shutter);
      delete fits;
    }
  result->calc_ratio();
  return result;
}


rl_result* fy_study(TTree *tsFLASHwf, const char* what="pmt_1",const char* other_cuts = "")
{
  const Int_t nrl = 6;
  Int_t rl[nrl] = {0,2,3,5,8,10};
  rl_result *result = new rl_result(nrl);
  result->SetTitle(";j;Y [ photons / MeV ]", ";j;Y [ photons / MeV ]");
 
  
  TString s_PMT_TB_FK = what;
  s_PMT_TB_FK.ToUpper();
  s_PMT_TB_FK += "_TB_FK";
  TString s_PMT_TS_FK = what;
  s_PMT_TS_FK.ToUpper();
  s_PMT_TS_FK += "_TS_FK";
 
  if(!tsFLASHwf->GetBranch(s_PMT_TB_FK) || !tsFLASHwf->GetBranch(s_PMT_TS_FK))
    {
      fprintf(stderr,"error: energy deposition tracing branches %s, %s not found in tree %s!\n",
	      s_PMT_TB_FK.Data(),s_PMT_TS_FK.Data(),tsFLASHwf->GetName());
      return 0;
    }
  
  for (Int_t i=0; i<nrl; i++)
    {


      // obtain the energy deposition tracing factors
      TString cut_exp = "beam_open_open&&rl==";
      cut_exp += rl[i];
      if(strlen(other_cuts))
	{
	  cut_exp += "&&";
	  cut_exp += other_cuts;
	}
      TString var_exp = s_PMT_TB_FK;
      var_exp +=":";
      var_exp +=  s_PMT_TS_FK;
      tsFLASHwf->Draw(var_exp,cut_exp,"goff");
      if(tsFLASHwf->GetSelectedRows() <  1)
	{
	  fprintf(stderr,"error: failed to get %s %s values for rl = %d from tree %s\n", 
		  s_PMT_TB_FK.Data(), s_PMT_TS_FK.Data(), rl[i], tsFLASHwf->GetName());
	  return 0;
	}
      
      // energy tracing factor for the PMT for active region constrained by the blind 
      // calculated by FLUKA simulation
      Double_t PMT_TB_FK = tsFLASHwf->GetV1()[0];
      
      // energy tracking factor for the PMT for active region constrained by the shutter
      // calculated by FLUKA simulation
      Double_t PMT_TS_FK = tsFLASHwf->GetV2()[0];

      // do the linear fitts to NPH vs Charge to determine the noise-subtracted slope values
      // for beam to blind and beam to shutter regions.  Rescale results by the inverse
      // of the corresponding energy tracing factors to get the fluorescence yield values.
      p_vs_q_fits_class *fits = get_p_vs_q_fits(tsFLASHwf,what,rl[i],other_cuts);
      result->g_beam_blind->SetPoint(i, (Double_t)rl[i], fits->slope_beam_blind/PMT_TB_FK);
      result->g_beam_blind->SetPointError(i, 0, fits->d_slope_beam_blind/PMT_TB_FK);
      result->g_beam_shutter->SetPoint(i, (Double_t)rl[i], fits->slope_beam_shutter/PMT_TS_FK);
      result->g_beam_shutter->SetPointError(i, 0, fits->d_slope_beam_shutter/PMT_TS_FK);
      
      // clean up
      delete fits;
    }
  
  // calculate the ratio of beam to shutter and beam to blind results
  result->calc_ratio();

  // return the answers
  return result;
}


TString get_unique_object_name(const char* basename = "obj")
{
  for (Int_t i=0; i<N_MAX_OBJ; i++)
    {
      TString name = basename;
      name += "_";
      name += i;
      // skip if already present
      if(gROOT->FindObject(name))
        continue;
      return name;
    }
  fprintf(stderr,"warning: ran out of object for basename %s\n",basename);
  return "";
}


Bool_t plotWF(TTree *tsFLASHwf,
	      const char* what   = "pmt_1",
	      Int_t ievent       = 0,
	      Int_t    nbins     = 500,
	      Double_t xlo_ns    = 0.0, 
	      Double_t xup_ns    = 500.0,
	      Bool_t draw_pulse  = true,
	      const char* xtitle = "t [ nS ]",
	      const char* ytitle = "V_{OUT} [ V ]")
{

  if(!tsFLASHwf)
    {
      fprintf(stderr,"error: tsFLASHwf has not been initialized!\n");
      return false;
    }


  Int_t nt0;
  Int_t nt1;
  Double_t wf[NT0_MAX];
  Double_t dt;  // width of the time slice as read out [ s ]

  
  TString wf_what = "wf_";
  wf_what += what;
  wf_what.ToLower();
  
  TString dt_what = "dt_";
  dt_what += what;
  dt_what.ToLower();

  if(!tsFLASHwf->GetBranch("nt0") || !tsFLASHwf->GetBranch("nt1"))
    {
      fprintf(stderr,"error: important branches nt0 nt1 missing in tree %s\n",tsFLASHwf->GetName());
      return false;
    }
  if(!tsFLASHwf->GetBranch(wf_what))
    {
      fprintf(stderr,"error: branch '%s' that corresponds to waveform of '%s' not present in the tree %s\n",
	      wf_what.Data(), what, tsFLASHwf->GetName());
      return false;
    }
  if(!tsFLASHwf->GetBranch(dt_what))
    {
      fprintf(stderr,"error: branch '%s' that corresponds to waveform of '%s' is not in tree %s\n",
	      dt_what.Data(), what, tsFLASHwf->GetName());
      return false;
    }
  tsFLASHwf->SetBranchAddress("nt0",     &nt0);
  tsFLASHwf->SetBranchAddress("nt1",     &nt1);
  tsFLASHwf->SetBranchAddress(wf_what,   wf);
  tsFLASHwf->SetBranchAddress(dt_what,   &dt);

  
  Pulse *pulse = 0; // non zero if there is the pulse analysis branch for the thing we are looking at
  TString pulse_branch_name = what;
  pulse_branch_name.ToLower();
  if(tsFLASHwf->GetBranch(pulse_branch_name))
    {
      pulse = new Pulse();
      tsFLASHwf->SetBranchAddress(pulse_branch_name,&pulse);
    }
  if(!tsFLASHwf->GetEntry(ievent))
    {
      fprintf(stderr,"failed to get ievent =  %d from tree %s\n", ievent, tsFLASHwf->GetName());
      return false;
    }
  // important so that it is OK to eventually get rid of objects to which the branches have been pointed
  tsFLASHwf->ResetBranchAddresses();
  
  // make sure that the two scopes have same settings
  if(nt0 != nt1)
    {
      fprintf(stderr,"strange event: number of time slices in scope 0 (%d) different from that of scope 1 (%d); ",
	      nt0,nt1);
      fprintf(stderr,"such events are not supported by this routine, plot their waveforms manually\n");
      return false;
    }
  Int_t nt = nt0;  // assuming same number of time slices for both scopes


  TString s_title = what;
  s_title.ToUpper();
  s_title += ";";
  s_title += xtitle;
  s_title += ";";
  s_title += ytitle;
  TProfile* p = new TProfile(get_unique_object_name("pWF"),s_title,nbins,xlo_ns,xup_ns);
  p->GetXaxis()->SetTitleSize(0.055);
  p->GetYaxis()->SetTitleSize(0.055);
  p->GetYaxis()->SetTitleOffset(1.1);
  p->SetLineColor(kBlack);

  for (Int_t i=0; i < nt; i++)
    p->Fill(1e9*((Double_t)i+0.5)*dt,wf[i]);


  init_canv();
  
  p->SetStats(0);
  p->Draw("hist");
  
  if(pulse && draw_pulse)
    {
      TLine *lpedCalc = new TLine(pulse->wf_i0,pulse->ped_mean_bp,(Double_t)pulse->nt_ped_bp,pulse->ped_mean_bp);
      lpedCalc->SetLineColor(kRed);
      lpedCalc->SetLineWidth(3);
      lpedCalc->Draw();


      TLine *lped = new TLine(pulse->wf_i0,pulse->ped_mean_bp,(Double_t)pulse->wf_nt,pulse->ped_mean_bp);
      lped->SetLineColor(kRed);
      lped->SetLineWidth(2);
      lped->SetLineStyle(2);
      lped->Draw();


      TLine *lped1 = new TLine(pulse->wf_i0,pulse->ped_mean_bp-pulse->ped_std_bp,
			       (Double_t)pulse->wf_nt,
			       pulse->ped_mean_bp-pulse->ped_std_bp);
      lped1->SetLineColor(kRed);
      lped1->SetLineWidth(2);
      lped1->SetLineStyle(2);
      lped1->Draw();

      TLine *lped2 = new TLine(pulse->wf_i0,pulse->ped_mean_bp+pulse->ped_std_bp,
			       (Double_t)pulse->wf_nt,
			       pulse->ped_mean_bp+pulse->ped_std_bp);
      lped2->SetLineColor(kRed);
      lped2->SetLineWidth(2);
      lped2->SetLineStyle(2);
      lped2->Draw();
      
    
      
      TLine *l_signal_i1 = new TLine((Double_t)pulse->istart,p->GetMinimum(),pulse->istart,p->GetMaximum());
      l_signal_i1->SetLineColor(kRed);
      l_signal_i1->SetLineWidth(2);
      l_signal_i1->Draw();
      
      TLine *l_signal_i2 = new TLine((Double_t)pulse->istart+(Double_t)pulse->width,
				     p->GetMinimum(),
				     (Double_t)pulse->istart+(Double_t)pulse->width,
				     p->GetMaximum());
      l_signal_i2->SetLineColor(kRed);
      l_signal_i2->SetLineWidth(2);
      l_signal_i2->Draw();
      
      
      TString title_and_s_s2n;
      TString tmp_what = what;
      tmp_what.ToUpper();
      title_and_s_s2n.Form("%s: SIGNAL / #sigma_{SIGNAL} = %.2E / %.2E = %.1f", 
			   tmp_what.Data(),
			   pulse->signal,pulse->d_signal,
			   pulse->signal/pulse->d_signal);
      p->SetTitle(title_and_s_s2n);
      
      delete pulse;
    }


  return true;
}


TGraphErrors *scan_sig_to_noise(TTree *tsFLASHwf, const char* cuts = "beam_open_open==1&&rl==5",
				Double_t step_size = 1.0,
				Double_t min  = 5.0,
				Double_t max  = 100.0)
{  
  gROOT->cd();
  if(gROOT->FindObject("hSigToChargeRatio"))
    delete (gROOT->FindObject("hSigToChargeRatio"));
  TH1D *hSigToChargeRatio = new TH1D("hSigToChargeRatio","",150,0.0,4500.0);
  TGraphErrors *g_sigma_to_mean_ratio = new TGraphErrors(0);
  Int_t npts = 0;
  for (Double_t cut_val = min; cut_val < max; cut_val += step_size)
    {
      Double_t sigma_to_mean_ratio = 0;
      Double_t sigma_to_mean_ratio_error = 0;
      hSigToChargeRatio->Reset();
      TString cut_exp;
      cut_exp.Form("%s&&pmt_1.signal/pmt_1.d_signal>%f&&coil.signal/coil.d_signal>%f",
		   cuts,cut_val,cut_val);
      tsFLASHwf->Draw("pmt_1.csignal/coil.csignal>>hSigToChargeRatio",cut_exp,"goff");
      hSigToChargeRatio->Fit("gaus","Q,0");
      TF1 *f = hSigToChargeRatio->GetFunction("gaus");
      if(f)
	{
	  Double_t mean = f->GetParameter(1);
	  Double_t mean_error = f->GetParError(1);
	  Double_t sigma = f->GetParameter(2);
	  Double_t sigma_error = f->GetParError(2);
	  sigma_to_mean_ratio = sigma / mean;
	  sigma_to_mean_ratio_error = sigma_to_mean_ratio * 
	    TMath::Sqrt(mean_error*mean_error/mean/mean+sigma_error*sigma_error/sigma/sigma);
	}
      g_sigma_to_mean_ratio->SetPoint(npts,cut_val,sigma_to_mean_ratio);
      g_sigma_to_mean_ratio->SetPointError(npts,0.5,sigma_to_mean_ratio_error);
      npts++;
    }
  g_sigma_to_mean_ratio->SetTitle(";Signal-to-noise ratio cut; Sigma / Mean of (PMT_1 Signal / Coil Charge)");
  g_sigma_to_mean_ratio->SetMarkerStyle(20);
  delete hSigToChargeRatio;
  return g_sigma_to_mean_ratio;
}

Int_t sflashAnalysis()
{
  fprintf(stdout,"have:\n");
  fprintf(stdout,"sflashAnalysis(infile,outfile) - produces output tree with analyzed signals\n");
  fprintf(stdout,"scan_sig_to_noise(tsFLASHwf)   - returns a scan of sigma/mean ratio vs signal to noise ratio cut\n");
  fprintf(stdout,"plotWF(TTree* tsFLASHwf,const char* what, Int_t ievent) - plots a waveform for what for event ievent\n");
  return 0;
}
