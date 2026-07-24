#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include <unistd.h>
#include "AutoCalibration.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TString.h"
#include "TGraph.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "setup.h"

using namespace std;

int AutoCalibration::Set_PeakEnergy(double* E,int num)
{
	peak_num=num;
	map_Predefine_Peak_Energy.clear();
	for(int i=0;i<num;i++)
		map_Predefine_Peak_Energy.insert(pair<double,double>(E[i],1.0));

	threshold=0.3;
	sigma=5.0;
	CHI_LIMIT=500.0;
	Range_L=0.01;
	Range_R=0.01;

	return peak_num;
}

int AutoCalibration::Set_PeakEnergy(const char* filename)
{
	FILE* fp_file=fopen(filename,"r");
	if(fp_file == NULL)
		return -1;
	double E;
	double bool_fit;
	peak_num=0;
	fscanf(fp_file,"%d",&pre_fit_num);
	fscanf(fp_file,"%lf",&EDiff_Select);
	std::map<double,double> map_fit_E;

#ifdef DEBUG_MSG_OFFLINE
	cout<<"prefitnum:"<<pre_fit_num<<",E_Select:"<<EDiff_Select<<endl;
#endif
	while(fscanf(fp_file,"%lf %lf",&E, &bool_fit)!=EOF)
	{
#ifdef DEBUG_MSG_OFFLINE
		cout<<"read from file E:"<<E<<",boolfit:"<<bool_fit<<endl;
#endif
		if(peak_num < pre_fit_num )
		{
			map_fit_E.emplace(E,1.0);
		}
		map_Predefine_Peak_Energy.emplace(E,bool_fit);
		peak_num++;
	}
	auto iter_tmp = map_fit_E.rbegin();
	int idx = 0;
	while(iter_tmp != map_fit_E.rend() && idx < MAX_PEAK_NUM)
	{
		Array_Predefine_Peak_Fit[idx] = iter_tmp->first;
		idx++;
		iter_tmp++;
	}

	stand_source = false;
	threshold=0.1;
	sigma=3.0;
	CHI_LIMIT=500.0;
	Range_L=0.005;
	Range_R=0.005;

	f1_gaus=new TF1("fit","gaus(0)",100,32100);
//	peak_num *=3;
	return peak_num;
}


void AutoCalibration::Set_StdAlpha_Energy()
{
	const int PEAK_NUM=3;
	double STD_Alpha_E[PEAK_NUM]={5156.59,5485.56,5804.77};
	Set_PeakEnergy(STD_Alpha_E,PEAK_NUM);

	Array_Predefine_Peak_Fit[0] = 5804.77;
	Array_Predefine_Peak_Fit[1] = 5485.56;
	pre_fit_num = 2;

	threshold=0.5;
	sigma=5.0;
	CHI_LIMIT=500.0;
	Range_L=0.003;
	Range_R=0.01;
	stand_source = true;

	f1_gaus=new TF1("fit","gaus(0)",100,32100);

}

void AutoCalibration::Set_StdGamma_Energy()
{
	int i;
	const Int_t   PEAK_NUM_152Eu=11;
	double STD_152Eu_E[PEAK_NUM_152Eu]={121.8, 244.7, 344.3, 411.1, 444.0, 778.9, 867.4, 964.1, 1085.8, 1112.0, 1408.0};
	//double STD_152Eu_Int[PEAK_NUM_152Eu]={28.37,7.53,26.57,2.23,3.12,12.97,4.214,14.63,10.13,1112.1,20.85};

	const Int_t   PEAK_NUM_133Ba=5;
	double STD_133Ba_E[PEAK_NUM_133Ba]={80.9, 276.4,302.8,356.0,383.8};
	//double STD_133Ba_Int[PEAK_NUM_133Ba]={34.11,7.147,18.30,61.94,8.905};

	//Int_t PEAK_NUM_60Co=2;
	//double STD_60Co_E[2]={1173.2,1332.5};
	//double STD_60Co_Int[2]={100,100};

	map_Predefine_Peak_Energy.clear();

	for(i=0;i<PEAK_NUM_152Eu;i++)
		map_Predefine_Peak_Energy.insert(pair<double,double>(STD_152Eu_E[i],1.0));
	for(i=0;i<PEAK_NUM_133Ba;i++)
		map_Predefine_Peak_Energy.insert(pair<double,double>(STD_133Ba_E[i],1.0));
//	for(i=0;i<PEAK_NUM_60Co;i++)
//		map_Predefine_Peak_Energy.insert(pair<double,double>(STD_60Co_E[i],1.0));

	peak_num=PEAK_NUM_152Eu+PEAK_NUM_133Ba;
	threshold=0.1;
	sigma=20.0;
	CHI_LIMIT=100.0;
	Range_L=0.01;
	Range_R=0.01;
	stand_source = true;

	Array_Predefine_Peak_Fit[0] = 1408.0;
	Array_Predefine_Peak_Fit[1] = 1112.0;
	pre_fit_num = 2;

	f1_gaus=new TF1("fit","gaus(0)+pol0(3)",100,32100);
}


int AutoCalibration::Get_PreFit_Par()
{
	int ret = -1;

	if(pre_fit_num == 0)
	{
		f1_linear->SetParameter(0,-60.0);
		f1_linear->SetParameter(1,1.0);
		ret = 0;
	}
	else if (pre_fit_num == 1)
	{
		std::cout<<"only one peak for pre fit....."<<std::endl;
		f1_linear->SetParameter(0,-60.0);
		f1_linear->SetParameter(1,1.0);
		ret = -1;
	}
	else
	{
		if (pre_fit_num == 2 && stand_source == true)
		{
			int i=0;
			auto iter_peak  = map_Spec_Peak_Energy.rbegin();
			while(iter_peak != map_Spec_Peak_Energy.rend())
			{
				Array_Spec_Peak_Fit[i] = iter_peak->first;
				++iter_peak;
				++i;
				if(i >= pre_fit_num)
					break;
			}
#ifdef DEBUG_MSG_OFFLINE
			cout<<"x1:"<< Array_Spec_Peak_Fit[1]<<endl;
			cout<<"x0:"<< Array_Spec_Peak_Fit[0]<<endl;

			cout<<"y1:"<< Array_Predefine_Peak_Fit[1]<<endl;
			cout<<"y0:"<< Array_Predefine_Peak_Fit[0]<<endl;
#endif


			TGraph* gr=new TGraph(pre_fit_num, Array_Spec_Peak_Fit, Array_Predefine_Peak_Fit);
			gr->Fit(f1_linear,"Q");
			ret = 0;
#ifdef DEBUG_MSG_OFFLINE
			cout<<"b:"<<f1_linear->GetParameter(0)<<",k:"<<f1_linear->GetParameter(1)<<endl;
#endif
		}
		else if (pre_fit_num == 2 && stand_source == false)
		{
			int i=0;
			auto iter_peak  = map_Spec_Peak_Energy.rbegin();
			while(iter_peak != map_Spec_Peak_Energy.rend())
			{
				Array_Spec_Peak_Fit[i] = iter_peak->first;
				++iter_peak;
				++i;
				if(i >= pre_fit_num)
					break;
			}
#ifdef DEBUG_MSG_OFFLINE
			cout<<"x1:"<< Array_Spec_Peak_Fit[1]<<endl;
			cout<<"x0:"<< Array_Spec_Peak_Fit[0]<<endl;

			cout<<"y1:"<< Array_Predefine_Peak_Fit[1]<<endl;
			cout<<"y0:"<< Array_Predefine_Peak_Fit[0]<<endl;
#endif

			TGraph* gr=new TGraph(pre_fit_num, Array_Spec_Peak_Fit, Array_Predefine_Peak_Fit);
			gr->Fit(f1_linear,"Q");
			ret = 0;
		}
		else if(pre_fit_num > 2)
		{
			int spec_peak_num  = 0;
			auto iter_peak  = map_Spec_Peak_Energy.rbegin();
			while(iter_peak != map_Spec_Peak_Energy.rend() && spec_peak_num < MAX_PEAK_NUM)
			{
				Array_Spec_Peak_Fit[spec_peak_num] = iter_peak->first;
				++iter_peak;
				++spec_peak_num;
			}
			/// search for the correct peak for pre-calibration
			int ii = 0, jj = 0, kk = 0;
			int ll = 0, mm = 0, nn = 0;
			double Peak_Y[3] = {0};
			double Peak_X[3] = {0};
			for(ii = 0; ii < pre_fit_num; ii++)
				for(jj = ii + 1; jj < pre_fit_num; jj++)
					for(kk = jj + 1; kk < pre_fit_num; kk++)
					{
						for(ll = 0; ll < spec_peak_num; ll++)
							for(mm = ll + 1; mm < spec_peak_num; mm++)
								for(nn = mm + 1; nn < spec_peak_num; nn++)
								{
#ifdef DEBUG_MSG_OFFLINE
									double delta_y2 = Array_Predefine_Peak_Fit[ii] - Array_Predefine_Peak_Fit[jj];
									double delta_y1 = Array_Predefine_Peak_Fit[jj] - Array_Predefine_Peak_Fit[kk];
									double delta_x2 = Array_Spec_Peak_Fit[ll] - Array_Spec_Peak_Fit[mm];
									double delta_x1 = Array_Spec_Peak_Fit[mm] - Array_Spec_Peak_Fit[nn];
									cout<<"x2:"<< Array_Spec_Peak_Fit[ll]<<endl;
									cout<<"x1:"<< Array_Spec_Peak_Fit[mm]<<endl;
									cout<<"x0:"<< Array_Spec_Peak_Fit[nn]<<endl;

									cout<<"y2:"<< Array_Predefine_Peak_Fit[ii]<<endl;
									cout<<"y1:"<< Array_Predefine_Peak_Fit[jj]<<endl;
									cout<<"y0:"<< Array_Predefine_Peak_Fit[kk]<<endl;
									cout<<"pre fit phase1:"<<(delta_y2/delta_x2 - delta_y1/delta_x1)/(delta_y2/delta_x2)<< endl;
#endif
									{
										Peak_X[0] = Array_Spec_Peak_Fit[ll];
										Peak_X[1] = Array_Spec_Peak_Fit[mm];
										Peak_X[2] = Array_Spec_Peak_Fit[nn];
										Peak_Y[0] = Array_Predefine_Peak_Fit[ii];
										Peak_Y[1] = Array_Predefine_Peak_Fit[jj];
										Peak_Y[2] = Array_Predefine_Peak_Fit[kk];
										TGraph* gr=new TGraph(3, Peak_X, Peak_Y);
										gr->Fit(f1_linear,"Q");
#ifdef DEBUG_MSG_OFFLINE
										cout<<"pre fit phase2: k:"<<f1_linear->GetParameter(1)<<", b:"<<f1_linear->GetParameter(0)<< endl;
#endif
										if(f1_linear->GetParameter(1) > 0.7 && f1_linear->GetParameter(1) < 1.1 && f1_linear->GetChisquare() < 5000)
											return 0;
									}
								}

					}
			for(ii = 0; ii < pre_fit_num; ii++)
				for(jj = ii + 1; jj < pre_fit_num; jj++)
					for(ll = 0; ll < spec_peak_num; ll++)
						for(mm = ll + 1; mm < spec_peak_num; mm++)
						{
							Peak_X[0] = Array_Spec_Peak_Fit[ll];
							Peak_X[1] = Array_Spec_Peak_Fit[mm];
							Peak_Y[0] = Array_Predefine_Peak_Fit[ii];
							Peak_Y[1] = Array_Predefine_Peak_Fit[jj];
							TGraph* gr=new TGraph(2, Peak_X, Peak_Y);
							gr->Fit(f1_linear,"Q");
							if(f1_linear->GetParameter(1) > 0.7 && f1_linear->GetParameter(1) < 1.1)
								return 0;
						}
		}
	}
	return ret;
}


int AutoCalibration::GetFit(TH1F* h, double* par)
{
	int i=0;
	ch_count++;
	if(!h)
	{
		cout<<"get hist error, check!"<<endl;
		return -1;
	}

	map_Spec_Peak_Energy.clear();

	//Use TSpectrum to find the peak candidates
#ifdef DEBUG_MSG_OFFLINE
	cout<<"check debug"<<endl;
	TCanvas* c1=new TCanvas("c1","c1");
	c1->cd();
	h->Draw();
	c1->Modified();
#endif
	TSpectrum *s = new TSpectrum(2*peak_num);
	Int_t nfound = s->Search(h,sigma,"",threshold);
	Double_t *xpos=s->GetPositionX();
	Double_t *ypos=s->GetPositionY();

#ifdef DEBUG_MSG_OFFLINE
	s->Print();
#endif

//	if(nfound > map_Predefine_Peak_Energy.size()*5)
//	{
//		cout<<"peak num found ("<< nfound <<") too much than peaks num" <<map_Predefine_Peak_Energy.size()<<", plz check the search parameters!"<<endl;
//		return -2;
//	}
	if(nfound<2)
	{
		cout<<"peak found too less, plz check the search parameters!"<<endl;
		return -3;
	}

	for(i = 0; i < nfound; ++i)
	{
		f1_gaus->SetParameter(2,sigma);
		f1_gaus->SetParameter(1,xpos[i]);
		f1_gaus->SetParameter(0,ypos[i]);
		h->Fit(f1_gaus,"Q","",(1-Range_L)*xpos[i],(1+Range_R)*xpos[i]);
		xpos[i]=f1_gaus->GetParameter(1); //
		map_Spec_Peak_Energy.emplace(xpos[i],f1_gaus->GetParameter(2));
#ifdef DEBUG_MSG_OFFLINE
		cout<<"fit peak "<<i<<" finish"<<endl;
		cout<<"Peak pos:"<<f1_gaus->GetParameter(1)<<",sigma:"<<f1_gaus->GetParameter(2)<<endl;
		c1->Modified();
	//	c1->WaitPrimitive();
#endif
	}
	///cal

	if(Get_PreFit_Par() < 0)
	{
		par[0] = 0.0;
		par[1] = 1.0 ;
		par[2] = 0.0;
		par[3] = 0xFFF;
		cout<<"pre fit failed"<<endl;
		return -1;
	}

	cout<<h->GetName()<<", pre fit parameters: k:"<<f1_linear->GetParameter(1)<<", b:"<<f1_linear->GetParameter(0)<< endl;

#ifdef DEBUG_MSG_OFFLINE
	c1->WaitPrimitive();
#endif

	//////
	auto iter_peak   = map_Spec_Peak_Energy.rbegin();
	par[2]=iter_peak->second;
	i = 0;
	std::vector<double> source_E;
	for(auto iter_source = map_Predefine_Peak_Energy.rbegin(); iter_source != map_Predefine_Peak_Energy.rend(); ++iter_source)
	{
		if(iter_source->second != 0)
			source_E.push_back(iter_source->first);
	}
	std::vector<bool> source_used(source_E.size(), false);
	while(iter_peak != map_Spec_Peak_Energy.rend() && i < MAX_PEAK_NUM)
	{
		double eval_E = f1_linear->Eval(iter_peak->first);
		int best_idx = -1;
		double best_diff = EDiff_Select;
		for(size_t idx = 0; idx < source_E.size(); ++idx)
		{
			if(source_used[idx])
				continue;
			double diff = std::fabs(eval_E - source_E[idx]);
			if(diff <= best_diff)
			{
				best_diff = diff;
				best_idx = idx;
			}
		}
		if(best_idx >= 0)
		{
			cal_x[i] = iter_peak->first;
			cal_y[i] = source_E[best_idx];
			source_used[best_idx] = true;
#ifdef DEBUG_MSG_OFFLINE
			cout<<"use prefit match x:"<<cal_x[i]<<", y:"<<cal_y[i]<<", eval:"<<eval_E<<", diff:"<<best_diff<<endl;
#endif
			i++;
		}
		iter_peak++;
	}


	if(i < 2)
	{
		par[0] = f1_linear->GetParameter(0);
		par[1] = f1_linear->GetParameter(1);
		par[2] = 0.0;
		par[3] = 0xFFE;
		cout<<"fit peak matched too less: "<<i<<endl;
		return -4;
	}
#ifdef DEBUG_MSG_OFFLINE
	for(int j = 0;j < i; j++)
		cout<<"use:"<<j<<":"<<cal_x[j]<<" "<<cal_y[j]<<endl;
#endif
	TGraph* gr2=new TGraph(i,cal_x,cal_y);
	gr2->Fit(f1_linear,"Q");
	par[0]=f1_linear->GetParameter(0);
	par[1]=f1_linear->GetParameter(1);
	par[2]=2.355*par[1]*par[2];
	par[3]=f1_linear->GetChisquare();
	cout<<h->GetName()<<": par value:"<<par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<par[3]<<endl;
#ifdef DEBUG_MSG_OFFLINE
	cout<<"curr pars: sigma:"<<sigma<<", threshold:"<<threshold<<", EDiff_Select:"<<EDiff_Select<<endl;
	cout<<"double click on canvas , continue..."<<endl;
	c1->WaitPrimitive();
#endif
	return 0;
}

