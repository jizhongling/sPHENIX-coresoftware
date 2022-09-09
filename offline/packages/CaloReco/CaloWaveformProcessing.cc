#include "CaloWaveformProcessing.h"

#include <string>
#include <iostream>

#include <TF1.h>
#include <TFile.h>
#include <TProfile.h>
// #include <torch/script.h> // One-stop header.
#include "Math/WrappedTF1.h"
#include "Math/WrappedMultiTF1.h"
#include "Fit/BinData.h"
#include "Fit/UnBinData.h"
#include "HFitInterface.h"
#include "Fit/Fitter.h"
#include "Fit/Chi2FCN.h"
#include "TThread.h"
#include <pthread.h>
#include <ROOT/TThreadExecutor.hxx>
#include <ROOT/TThreadedObject.hxx>

TProfile* CaloWaveformProcessing::h_template = nullptr;


double CaloWaveformProcessing::template_function(double *x, double *par)
{
  Double_t v1 = par[0]*h_template->Interpolate(x[0]-par[1])+par[2];
  return v1;
}





void CaloWaveformProcessing::initialize_processing()
{


  if (m_modeltype == 1)
    {
      TFile* fin = TFile::Open(m_template_input_file.c_str());
      assert(fin);
      assert(fin->IsOpen());
      h_template = static_cast<TProfile*>(fin->Get("hp_electrons_fine_emcal_36_8GeV"));
    }
}



std::vector<std::vector<float>> CaloWaveformProcessing::process_waveform(std::vector<std::vector<float>> waveformvector)
{

int size1 = waveformvector.size();
 std::vector<std::vector<float>> fitresults;
if (m_modeltype == 1)
  {
    for (int i = 0; i < size1;i++)
      {
	waveformvector.at(i).push_back(i);
      }
    fitresults =  CaloWaveformProcessing::calo_processing_templatefit(waveformvector);
  }
  return fitresults;
}





std::vector<std::vector<float>> CaloWaveformProcessing::calo_processing_templatefit(std::vector<std::vector<float>> chnlvector)
{
  ROOT::TThreadExecutor t(_nthreads);
  auto func = [&](std::vector<float> &v) {
    int size1 = v.size()-1;
    auto h = new TH1F(Form("h_%d",(int)round(v.at(size1))),"",size1,0,size1);
    float maxheight = 0;
    int maxbin = 0;
    for (int i = 0; i < size1;i++)
      {
	h->SetBinContent(i+1,v.at(i));
	if (v.at(i)>maxheight)
	  {
	    maxheight = v.at(i);
	    maxbin = i;
	  }
      }
    float pedestal = 1500;
    if (maxbin > 4)
      {
	pedestal=  0.5* ( v.at(maxbin - 4) + v.at(maxbin-5));
      }
    else if (maxbin > 3)
      {
	pedestal=( v.at(maxbin - 4));
      }
    else 
      {
	pedestal = 0.5* ( v.at(size1-3) + v.at(size1-2)); 
      }
    auto f = new TF1(Form("f_%d",(int)round(v.at(size1))),template_function,0,31,3);
    ROOT::Math::WrappedMultiTF1 * fitFunction = new ROOT::Math::WrappedMultiTF1(*f, 3 );
    ROOT::Fit::BinData data(v.size()-1,1);
    ROOT::Fit::FillData(data,h);
    ROOT::Fit::Chi2Function *EPChi2 = new ROOT::Fit::Chi2Function(data, *fitFunction);
    ROOT::Fit::Fitter *fitter = new ROOT::Fit::Fitter();
    fitter->Config().MinimizerOptions().SetMinimizerType("GSLMultiFit");
    double params[] = {static_cast<double>(maxheight),static_cast<double>(maxbin-5),static_cast<double>(pedestal)};
    fitter->Config().SetParamsSettings(3,params);
    fitter->FitFCN(*EPChi2,0,data.Size(),true);
    for (int i =0;i<3;i++)
      {
	v.push_back(f->GetParameter(i));
      }
    h->Delete();
    f->Delete();
    delete fitFunction;
    delete fitter;
    delete EPChi2;
  };
  
  t.Foreach(func, chnlvector);
  
  int size3 = chnlvector.size();
  std::vector<std::vector<float>> fit_params;
  std::vector<float> fit_params_tmp;
  for (int i = 0; i < size3;i++)
    {
      std::vector<float> tv = chnlvector.at(i);
      int size2 = tv.size();
      for (int q = 3; q > 0;q--)
	{
	  fit_params_tmp.push_back(tv.at(size2-q));
	}
      fit_params.push_back(fit_params_tmp);
      fit_params_tmp.clear();
    }


  chnlvector.clear();
  return fit_params;
  fit_params.clear();
}



