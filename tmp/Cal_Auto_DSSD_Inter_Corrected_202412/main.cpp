#include <iostream>
#include <fstream>
#include <vector>
#include "TSystem.h"
#include "TApplication.h"
#include "TString.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include "setup.h"
#include "AutoCalibration.h"

using namespace std;

int main(int argc, char* argv[])
{
	TApplication* app=new TApplication("app",0,0);
	TFile *file_in;
	int debug_ch = -1;
	TString hist_filename, energy_filename, str_tmp;

	// open the calibration file and get tree
#ifdef DEBUG_MSG_OFFLINE
	if(argc == 4)
	{
		hist_filename = argv[1];
		energy_filename = argv[2];
		debug_ch = atoi(argv[3]);
	}
	else
	{
		cout<<"input hist filename"<<endl;
		cin>>hist_filename;
		cout<<"input energy filename"<<endl;
		cin>>energy_filename;
		cout<<"debug ch:"<<endl;
		cin>>debug_ch;
	}
	if(debug_ch < 0 || debug_ch > TOTAL_CH)
	{
		cout<<"debug ch:"<<debug_ch<<" error!"<<endl;
		return -1;
	}

#else
	if(argc == 3)
	{
		hist_filename = argv[1];
		energy_filename = argv[2];
	}
	else
	{
		cout<<"input hist filename"<<endl;
		cin>>hist_filename;
		cout<<"input Energy filename"<<endl;
		cin>>energy_filename;
	}
#endif


	///////////////////
	file_in =new TFile(hist_filename.Data());
	if(file_in->IsZombie())
	{
		cout<<"open "<<hist_filename<<" error, plz check!"<<endl;
		return -2;
	}

	int ii,jj;
	TH1F *h1[TOTAL_CH];
	for(ii=0;ii<TOTAL_CH;ii++)
	{
		str_tmp=TString::Format("h1_E%03d",ii);
		h1[ii]=(TH1F*)file_in->Get(str_tmp.Data());
	}


	double Cal_Par[4];
	int ret;


	AutoCalibration* ptr_Auto=new AutoCalibration();
	//ptr_Auto->Set_StdAlpha_Energy();
	if(ptr_Auto->Set_PeakEnergy(energy_filename.Data()) <= 0)
	{
		cout<<"read peak energy file "<<energy_filename<<" error, plz check!"<<endl;
		return -3;
	}
	ptr_Auto->Set_SearchPar(8,0.1,500);
	ptr_Auto->Set_FitRange(0.003,0.005);
#ifdef DEBUG_MSG_OFFLINE
	char test='y';
	double sigma = 3.0;
	double threshold=0.14;

	h1[debug_ch]->Rebin(4);
	if(debug_ch <= DSSD_Y_CH_HIGH)
	{
		sigma = 3.0;
		threshold =0.6;
		ptr_Auto->Set_SearchPar(sigma,threshold,500);
		ptr_Auto->Set_FitRange(0.003,0.005);
	}
	else
	{
		sigma = 3.0;
		threshold =0.6;
		ptr_Auto->Set_SearchPar(sigma,threshold,500);
		ptr_Auto->Set_FitRange(0.003,0.005);
	}

	while(true)
	{
		ptr_Auto->GetFit(h1[debug_ch],Cal_Par);
		str_tmp=str_tmp.Format("%-d\t%.2f\t%.6f\t%.2f\t%-g", debug_ch, Cal_Par[0],Cal_Par[1],Cal_Par[2],Cal_Par[3]);
		cout<<str_tmp<<endl;
		cout<<"Debug continue(Y/N)"<<endl;
		cin>>test;
		if(test != 'Y' && test != 'y')
			break;
		else
		{
			cout<<"debug ch num:"<<endl;
			cin>>debug_ch;
			if(debug_ch <0 && debug_ch >= TOTAL_SI_CH)
				break;
			cout<<"sigma value("<<sigma<<"):"<<endl;
			cin>>sigma;
			cout<<"threshold value("<<threshold<<"):"<<endl;
			cin>>threshold;
			ptr_Auto->Set_SearchPar(sigma,threshold,500);
		}
	}
#else
	str_tmp=hist_filename.Remove(hist_filename.Sizeof()-6,5);
	TString filename_output="ener_cal_"+str_tmp+".root";
	str_tmp="ener_cal_"+str_tmp+".dat";
	ofstream ofs(str_tmp.Data());
	///////////////////
	TFile *file_output =new TFile(filename_output.Data(),"RECREATE");
	if(file_output->IsZombie())
	{
		cout<<"open "<<hist_filename<<" error, plz check!"<<endl;
	}
	file_output->cd();

	TH1F *h1_output[TOTAL_CH];
	TH1F *h1_total = new TH1F("h1_total","h1_total", 30000, 0, 30000);
	TH1F *h1_Resolution=new TH1F("h1_Resolution","h1_Resolution",303,0,303);


	for(ii=0;ii<TOTAL_CH;ii++)
	{
		str_tmp=TString::Format("h1_E%03d",ii);
		h1_output[ii]= new TH1F(str_tmp.Data(),str_tmp.Data(),30000,0,30000);
	}


	for( ii=0; ii <= DSSD_Y_CH_HIGH; ++ii)
	{
		//DSSDX
		if(ii < DSSD_Y_CH_LOW)
		{
			h1[ii]->Rebin(4);
			ptr_Auto->Set_SearchPar(3.0,0.10,500);
			ptr_Auto->Set_FitRange(0.003,0.005);
		}
		//DSSD_Y
		else if(ii <= DSSD_Y_CH_HIGH)
		{
			h1[ii]->Rebin(4);
			ptr_Auto->Set_SearchPar(3.0,0.10,500);
			ptr_Auto->Set_FitRange(0.003, 0.005);
		}

		ret=ptr_Auto->GetFit(h1[ii],Cal_Par);
		if(ret < 0)
		{
			cout<<"autocal for ch"<< ii<<"fail!! ret:"<<ret<<endl;
			Cal_Par[0]=0.0;
			Cal_Par[1]=1.0;
			Cal_Par[2]=0.0;
			Cal_Par[3]=0.0;
		}
		h1_Resolution->SetBinContent(ii+1, Cal_Par[2]);
		str_tmp=str_tmp.Format("%-d\t%10.2f\t%10.6f\t%6.2f\t%-g",ii,Cal_Par[0],Cal_Par[1],Cal_Par[2],Cal_Par[3]);
		ofs<<str_tmp<<endl;
		for(jj=0; jj<h1[ii]->GetNbinsX(); ++jj)
		{
				h1_output[ii]->Fill(h1[ii]->GetBinCenter(jj+1)*Cal_Par[1]+Cal_Par[0], h1[ii]->GetBinContent(jj+1));
		}
		(*h1_total) = (*h1_total)+(*h1_output[ii]);
	}

//	h1_total->DrawCopy("hist");
	if(file_output != NULL)
	{
		file_output->cd();
		for(ii=0;ii<TOTAL_CH;ii++)
			h1_output[ii]->Write("",TObject::kOverwrite);
		h1_total->Write("",TObject::kOverwrite);
		h1_Resolution->Write("",TObject::kOverwrite);
	}
	file_output->Close();
#endif
//	cout<<"CTRL+C quit"<<endl;
//	app->Run();
	return 0;
}

