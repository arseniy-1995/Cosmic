//#include <iostream.h>
#include "TFile.h"
#include "TTree.h"


Double_t S1 = 0.0;
Double_t S2 = 0.0;
Double_t X1 = 0.0;
Double_t X2 = 0.0;


const Int_t Smooth_index = 2;	// количество прогонов сглаживания гистограммы

const Double_t L_scint = 18; // размер сцинтиллятора, см


// П-образная функция для фита, 5 параматров
Double_t fit_ln(Double_t* x, Double_t* par) {

	Double_t arg1 = par[1] * x[0] + par[2];
	Double_t arg2 = -par[3] * x[0] + par[4];

	Double_t fitval = par[0] / ((1.0 + TMath::Exp(arg1)) * (1.0 + TMath::Exp(arg2)));

	return fitval;
}


TF1* func_1 = new TF1("my_fit_1", fit_ln, -5, 5, 5);


// расстояние между точками для логарифма отношения, на полувысоте
Double_t delta_L05(Double_t S1, Double_t X1, Double_t S2, Double_t X2) {

	Double_t temp = fabs(-X1 / S1) + fabs(X2 / S2);

	printf("-X1/S1=%.3f   -X2/S2=%.3f\n", -X1 / S1, -X2 / S2);

	printf("L0.5= %.3f\n", temp);

	return temp;

}

// расстояние между точками для логарифма отношения, на полувысоте
Double_t delta_L(Double_t S1, Double_t X1, Double_t S2, Double_t X2) {

	Double_t u1 = -((X1 + log(5.0 - 2.0 * sqrt(6.0))) / S1);

	Double_t u2 = ((X2 - log(5.0 + 2.0 * sqrt(6.0))) / S2);

	printf("u1=%.3f   u2=%.3f\n", u1, u2);


	Double_t temp = fabs(u1) + fabs(u2);

	printf("L= %.3f\n", temp);

	return temp;

}



Double_t delta_t01 (Double_t b01){
 
return 1000.0/(b01-50.0)+400.0;
  
}



Double_t delta_t02 (Double_t b02){
  
return 1000.0/(b02-50.0)+400.0;
  
}





Double_t bt01_left (Double_t *x, Double_t *par){
 
  Double_t xx=x[0];
  
  return 1000.0/(xx-260.0)+30.0;
  
}

Double_t bt01_right (Double_t *x, Double_t *par){
 
  Double_t xx=x[0];
  return 1000.0/(xx-400.0)+50.0;
  
}

Double_t bt02_left (Double_t *x, Double_t *par){
 
  Double_t xx=x[0];
  return 1000.0/(xx-260.0)+30.0;
  
}

Double_t bt02_right (Double_t *x, Double_t *par){
 
  Double_t xx=x[0];
  return 1000.0/(xx-400.0)+50.0;
  
}

const Double_t b01_pedestal = 0.0;	   // пьедисталы для вычета из амплитуд
const Double_t b02_pedestal = 0.0;

int main() {

 
  
//TFile *f=new TFile("stend_12.02.20.root");

TFile* f = new TFile("../macros/stend_09.12.20_cosmic.root");



TTree* ntpl = (TTree*)f->Get("ntpl");
  
TCanvas *MyC = new TCanvas("MyC", "MyC", 0, 0, 2000, 1000);

TF1 *f01_left = new TF1 ("f01_left", bt01_left,500, 2000, 1);
TF1 *f01_right = new TF1 ("f01_right",bt01_right,850, 2000, 1);

TF1 *f02_left = new TF1 ("f02_left", bt02_left,520, 2000, 1);
TF1 *f02_right = new TF1 ("f02_right", bt02_right,820, 2000, 1);



TH1F *h1_time = new TH1F("h1_time","t01-t02",200,-300,300);

TH1F *h1_time_delta = new TH1F("h1_time_delta","(t01-delta_t01)-(t02-delta_t02)",200,-300,300);



TH1F *h1_time_1 = new TH1F("h1_time_1","t01+t02",200,-300,300);

TH1F *h1_time_delta_1 = new TH1F("h1_time_delta_1","(t01-delta_t01)+(t02-delta_t02)",200,-300,300);


//TH1F *h1_time_delta_1 = new TH1F("h1_time_delta_1","time resolution_delta",500,-500,500);

//TH1F *h1_time_delta_2 = new TH1F("h1_time_delta_2","time resolution_delta",500,-500,500);

//TH1F *h1_time_delta_3 = new TH1F("h1_time_delta_3","time resolution_delta",500,-500,500);

//TH1F *h1_time_delta_4 = new TH1F("h1_time_delta_4","time resolution_delta",500,-500,500);



TH1F *h1_time_SiPM01_delta_1 = new TH1F("h1_time_SiPM01_delta_1","time resolution_delta",200,-300,300);

TH1F *h1_time_SiPM01_delta_2 = new TH1F("h1_time_SiPM01_delta_2","time resolution_delta",200,-300,300);

TH1F *h1_time_SiPM01_delta_3 = new TH1F("h1_time_SiPM01_delta_3","time resolution_delta",200,-300,300);

TH1F *h1_time_SiPM01_delta_4 = new TH1F("h1_time_SiPM01_delta_4","time resolution_delta",200,-300,300);


TH1F *h1_time_SiPM02_delta_1 = new TH1F("h1_time_SiPM02_delta_1","time resolution_delta",200,-300,300);

TH1F *h1_time_SiPM02_delta_2 = new TH1F("h1_time_SiPM02_delta_2","time resolution_delta",200,-300,300);

TH1F *h1_time_SiPM02_delta_3 = new TH1F("h1_time_SiPM02_delta_3","time resolution_delta",200,-300,300);

TH1F *h1_time_SiPM02_delta_4 = new TH1F("h1_time_SiPM02_delta_4","time resolution_delta",200,-300,300);


TH1F *h1_b1b2 = new TH1F("h1_b1b2","A1/A2",200,0,3.0);

TH1F* h1_b1b2_log = new TH1F("h1_b1b2_log", "0.5*log(A1/A2)", 200, -2.0, 2.0);


//cout<<"!!!"<<endl;

TH1F *h1_1 = new TH1F("h1_1","t01-delta_t01",200,-300,300);
TH1F *h1_2 = new TH1F("h1_2","t02-delta_t02",200,-300,300);



TH2 *h2_1 = new TH2D("h2_1","b01:t01+Cut",500,0,2100,500,0,2100);
TH2 *h3_1 = new TH2D("h3_1","b01:t01+Cut",500,0,2100,1000,0,2100);

TH2 *h2_2 = new TH2D("h2_2","b02:t02+Cut",500,0,2100,500,0,2100);
TH2 *h3_2 = new TH2D("h3_2","b02:t02+Cut",500,0,2100,500,0,2100);

TH2 *h4 = new TH2D("h4","b01:b02+Cut",500,0,2100,500,0,2100);
TH2 *h4_ = new TH2D("h4_","b01:b02+Cut",500,0,2100,500,0,2100);

TH2 *h4_2 = new TH2D("h4_2","t01:t02+Cut",500,0,2100,500,0,2100);
TH2 *h4_3 = new TH2D("h4_3","t01:t02+Cut",500,0,2100,500,0,2100);




//h2 = new TProfile2D("h2","h2",40,-4,4,40,-4,4,0,20);

//TCut cut_bb = "(b01-100.0)*(b02-100.0)>2500.0&&b01>100.0&&b02>100.0";

TCut cut_bb = "";


//TCut cut_time = "t01>250.0&&t02>250.0&&t01<400.0&&t02<400.0";

TCut cut_time = "";
//TCut cut_aplitude = "b01>80.0&&b02>80.0&&b01<1130.0&&b02<1130.0&&b05>80.0&&b05<900.0";


TCut cut_aplitude = "b01>200.0&&b02>200.0&&b01<450.0&&b02<450.0&&b05>80.0&&b05<900.0";


//TCut cut_aplitude = "b01>200.0&&b02>200.0&&b01<1100.0&&b02<1100.0&&b05>80.0&&b05<900.0";

//TCut cut_aplitude = "b01>40.0&&b02>40.0&&b01<1800.0&&b02<1800.0&&b05>30.0&&b05<1000.0&&(((b01>800.0&&b01<1000.0)||(b02>800.0&&b02<1000.0))&&b01<1000.0&&b02<1000.0)";



//TCut cut_aplitude = "b04<1980.0&&b05<1980.0&&b06<1980.0&&b07<1980.0&&b08>80.0&&b08<1800.0";

TCut cut_bt01 = "(t01-260.0)*(b01-30.0)>1000.0&&(t01-400.0)*(b01-50.0)<1000.0";
TCut cut_bt02 = "(t02-260.0)*(b02-30.0)>1000.0&&(t02-400.0)*(b02-50.0)<1000.0";





//TCut cut_time = "";
//TCut cut_aplitude = "";
//TCut cut_bt01 = "";
//TCut cut_bt02 = "";



TCut cut_b_interval_1 = "b01<500.0&&b02<500.0";

TCut cut_b_interval_2 = "b01>500.0&&b02>500.0&&b01<1000.0&&b02<1000.0";

TCut cut_b_interval_3 = "b01>1000.0&&b02>1000.0&&b01<1500.0&&b02<1500.0";

TCut cut_b_interval_4 = "b01>1500.0&&b02>1500.0&&b01<2000.0&&b02<2000.0";


TCut cut_b_interval_1_SiPM01 = "b01<200.0";

TCut cut_b_interval_2_SiPM01 = "b01>200.0&&b01<400.0";

TCut cut_b_interval_3_SiPM01 = "b01>400.0&&b01<600.0";

TCut cut_b_interval_4_SiPM01 = "b01>600.0&&b01<800.0";


TCut cut_b_interval_1_SiPM02 = "b02<200.0";

TCut cut_b_interval_2_SiPM02 = "b02>200.0&&b02<400.0";

TCut cut_b_interval_3_SiPM02 = "b02>600.0&&b02<800.0";

TCut cut_b_interval_4_SiPM02 = "b02>800.0&&b02<1000.0";






ntpl->Draw("t01-t02>>h1_time",cut_time&&cut_aplitude&&cut_bt01&&cut_bt02&&cut_bb);

ntpl->Draw("(t01-delta_t01(b01))-(t02-delta_t02(b02))>>h1_time_delta",cut_time&&cut_aplitude&&cut_bt01&&cut_bt02&&cut_bb);


ntpl->Draw("t01+t02-600.0>>h1_time_1",cut_time&&cut_aplitude&&cut_bt01&&cut_bt02);

ntpl->Draw("(t01-delta_t01(b01))+(t02-delta_t02(b02))+120>>h1_time_delta_1",cut_time&&cut_aplitude&&cut_bt01&&cut_bt02&&cut_bb);



//ntpl->Draw("(t04-delta_t04(b04))+(t05-delta_t05(b05))-((t06-delta_t06(b06))+(t07-delta_t07(b07)))+275.0>>h1_time_delta_1",cut_time&&cut_aplitude&&cut_bt04&&cut_bt05&&cut_bt06&&cut_bt07&&cut_b_interval_1);

//ntpl->Draw("(t04-delta_t04(b04))+(t05-delta_t05(b05))-((t06-delta_t06(b06))+(t07-delta_t07(b07)))+275.0>>h1_time_delta_2",cut_time&&cut_aplitude&&cut_bt04&&cut_bt05&&cut_bt06&&cut_bt07&&cut_b_interval_2);

//ntpl->Draw("(t04-delta_t04(b04))+(t05-delta_t05(b05))-((t06-delta_t06(b06))+(t07-delta_t07(b07)))+275.0>>h1_time_delta_3",cut_time&&cut_aplitude&&cut_bt04&&cut_bt05&&cut_bt06&&cut_bt07&&cut_b_interval_3);

//ntpl->Draw("(t04-delta_t04(b04))+(t05-delta_t05(b05))-((t06-delta_t06(b06))+(t07-delta_t07(b07)))+275.0>>h1_time_delta_4",cut_time&&cut_aplitude&&cut_bt04&&cut_bt05&&cut_bt06&&cut_bt07&&cut_b_interval_4);


ntpl->Draw("t01-delta_t01(b01)>>h1_1",cut_time&&cut_aplitude&&cut_bt01&&cut_bt02&&cut_bb);
ntpl->Draw("t02-delta_t02(b02)>>h1_2",cut_time&&cut_aplitude&&cut_bt01&&cut_bt02&&cut_bb);



ntpl->Draw("t01-delta_t01(b01)>>h1_time_SiPM01_delta_1",cut_time&&cut_aplitude&&cut_bt01&&cut_bt02&&cut_bb&&cut_b_interval_1_SiPM01);
ntpl->Draw("t01-delta_t01(b01)>>h1_time_SiPM01_delta_2",cut_time&&cut_aplitude&&cut_bt01&&cut_bt02&&cut_bb&&cut_b_interval_2_SiPM01);
ntpl->Draw("t01-delta_t01(b01)>>h1_time_SiPM01_delta_3",cut_time&&cut_aplitude&&cut_bt01&&cut_bt02&&cut_bb&&cut_b_interval_3_SiPM01);
ntpl->Draw("t01-delta_t01(b01)>>h1_time_SiPM01_delta_4",cut_time&&cut_aplitude&&cut_bt01&&cut_bt02&&cut_bb&&cut_b_interval_4_SiPM01);


ntpl->Draw("t02-delta_t02(b02)>>h1_time_SiPM02_delta_1",cut_time&&cut_aplitude&&cut_bt01&&cut_bt02&&cut_bb&&cut_b_interval_1_SiPM02);
ntpl->Draw("t02-delta_t02(b02)>>h1_time_SiPM02_delta_2",cut_time&&cut_aplitude&&cut_bt01&&cut_bt02&&cut_bb&&cut_b_interval_2_SiPM02);
ntpl->Draw("t02-delta_t02(b02)>>h1_time_SiPM02_delta_3",cut_time&&cut_aplitude&&cut_bt01&&cut_bt02&&cut_bb&&cut_b_interval_3_SiPM02);
ntpl->Draw("t02-delta_t02(b02)>>h1_time_SiPM02_delta_4",cut_time&&cut_aplitude&&cut_bt01&&cut_bt02&&cut_bb&&cut_b_interval_4_SiPM02);






ntpl->Draw("b01:b02>>h4", cut_bt01&&cut_bt02&&cut_time&&cut_aplitude&&cut_bb);
ntpl->Draw("b01:b02>>h4_");


ntpl->Draw("(b01-10)/(b02-20)>>h1_b1b2", cut_bt01&&cut_bt02&&cut_time&&cut_aplitude&&cut_bb);


ntpl->Draw("t01:t02>>h4_2", cut_bt01&&cut_bt02&&cut_time&&cut_aplitude&&cut_bb);
ntpl->Draw("t01:t02>>h4_3");


MyC->Divide(4,3);


MyC->cd(1);

gStyle->SetOptFit(1111);
h1_time->Fit("gaus", "","", -200,200);

h1_time->GetXaxis()->SetTitle("Time, channel number");
h1_time->GetYaxis()->SetTitle("Events number");
h1_time->GetXaxis()->CenterTitle();
h1_time->GetYaxis()->CenterTitle();

h1_time->Draw();


MyC->cd(2);

ntpl->Draw("b01:t01>>h2_1");
ntpl->Draw("b01:t01>>h3_1",cut_bt01&&cut_bt02&&cut_time&&cut_aplitude&&cut_bb);
//ntpl->Draw("b04:t04>>h3",cut_bt04);

h2_1->GetXaxis()->SetTitle("Time, channel number");
h2_1->GetYaxis()->SetTitle("Amplitude, channel number");
h2_1->GetXaxis()->CenterTitle();
h2_1->GetYaxis()->CenterTitle();

h2_1->Draw();
//h3->SetMarkerStyle(5);
h3_1->SetMarkerColor(kRed);
h3_1->Draw("same");

//f04_left->SetMarkerColor(kGreen);
//f04_right->SetMarkerColor(kGreen);

//f01_left->Draw("same");
//f01_right->Draw("same");

MyC->cd(3);

ntpl->Draw("b02:t02>>h2_2");
ntpl->Draw("b02:t02>>h3_2",cut_bt01&&cut_bt02&&cut_time&&cut_aplitude&&cut_bb);
//ntpl->Draw("b04:t04>>h3",cut_bt04);

h2_2->GetXaxis()->SetTitle("Time, channel number");
h2_2->GetYaxis()->SetTitle("Amplitude, channel number");
h2_2->GetXaxis()->CenterTitle();
h2_2->GetYaxis()->CenterTitle();

h2_2->Draw();
//h3->SetMarkerStyle(5);
h3_2->SetMarkerColor(kRed);
h3_2->Draw("same");

//f02_left->Draw("same");
//f02_right->Draw("same");


MyC->cd(4);


gStyle->SetOptFit(1111);
h1_1->Fit("gaus", "","", -200,200);

h1_1->GetXaxis()->SetTitle("Time, channel number");
h1_1->GetYaxis()->SetTitle("Events number");
h1_1->GetXaxis()->CenterTitle();
h1_1->GetYaxis()->CenterTitle();

h1_1->Draw();



MyC->cd(5);


gStyle->SetOptFit(1111);
h1_2->Fit("gaus", "","", -200,200);

h1_2->GetXaxis()->SetTitle("Time, channel number");
h1_2->GetYaxis()->SetTitle("Events number");
h1_2->GetXaxis()->CenterTitle();
h1_2->GetYaxis()->CenterTitle();

h1_2->Draw();

  

MyC->cd(6);


h1_time_delta->GetXaxis()->SetTitle("Time, channel number");
h1_time_delta->GetYaxis()->SetTitle("Events number");
h1_time_delta->GetXaxis()->CenterTitle();
h1_time_delta->GetYaxis()->CenterTitle();

gStyle->SetOptFit(1111);
h1_time_delta->Fit("gaus", "","", -200,200);
h1_time_delta->Draw();

MyC->cd(7);

h4_3->GetXaxis()->SetTitle("Time, channel number");
h4_3->GetYaxis()->SetTitle("Time, channel number");
h4_3->GetXaxis()->CenterTitle();
h4_3->GetYaxis()->CenterTitle();

h4_3->Draw();
//h3->SetMarkerStyle(5);
h4_2->SetMarkerColor(kRed);
h4_2->Draw("same");

MyC->cd(8);


TString var7 = Form("0.5*log((b01 - %f)/(b02 - %f)) >> h1_b1b2_log", b01_pedestal, b02_pedestal);

h1_b1b2_log->Smooth(Smooth_index);	 // Сглаживание


ntpl->Draw(var7,"sqrt(b01*b02)>100");

gStyle->SetOptFit(1111);
h1_b1b2_log->Fit("my_fit_1", "", "", -5, 5);
h1_b1b2_log->Draw();

S1 = func_1->GetParameter(1);
X1 = func_1->GetParameter(2);
S2 = func_1->GetParameter(3);
X2 = func_1->GetParameter(4);

printf("lambda_scint= %.3f cm\n", L_scint / delta_L(S1, X1, S2, X2));



MyC->cd(9);

h4_->GetXaxis()->SetTitle("Amplitude, channel number");
h4_->GetYaxis()->SetTitle("Amplitude, channel number");
h4_->GetXaxis()->CenterTitle();
h4_->GetYaxis()->CenterTitle();

h4_->Draw();
//h3->SetMarkerStyle(5);
h4->SetMarkerColor(kRed);
h4->Draw("same");

MyC->cd(10);


h1_b1b2->GetXaxis()->SetTitle("A1/A2");
h1_b1b2->GetYaxis()->SetTitle("Events number");
h1_b1b2->GetXaxis()->CenterTitle();
h1_b1b2->GetYaxis()->CenterTitle();

gStyle->SetOptFit(1111);
h1_b1b2->Fit("gaus", "","", 0.4,2.5);
h1_b1b2->Draw();




MyC->cd(11);

h1_time_1->GetXaxis()->SetTitle("Time, channel number");
h1_time_1->GetYaxis()->SetTitle("Events number");
h1_time_1->GetXaxis()->CenterTitle();
h1_time_1->GetYaxis()->CenterTitle();

gStyle->SetOptFit(1111);
h1_time_1->Fit("gaus", "","", -200,200);
h1_time_1->Draw();

MyC->cd(12);

h1_time_delta_1->GetXaxis()->SetTitle("Time, channel number");
h1_time_delta_1->GetYaxis()->SetTitle("Events number");
h1_time_delta_1->GetXaxis()->CenterTitle();
h1_time_delta_1->GetYaxis()->CenterTitle();

gStyle->SetOptFit(1111);
h1_time_delta_1->Fit("gaus", "","", -200,200);
h1_time_delta_1->Draw();










/*

MyC->cd(13);

gStyle->SetOptFit(1111);
h1_time_delta_3->Fit("gaus", "","", -200,200);
h1_time_delta_3->Draw();

MyC->cd(14);

gStyle->SetOptFit(1111);
h1_time_delta_4->Fit("gaus", "","", -200,200);
h1_time_delta_4->Draw();

*/


/*

MyC->cd(13);

gStyle->SetOptFit(1111);
h1_time_SiPM01_delta_1->Fit("gaus", "","", -200,500);
h1_time_SiPM01_delta_1->Draw();

MyC->cd(14);
1
gStyle->SetOptFit(1111);
h1_time_SiPM01_delta_2->Fit("gaus", "","", -200,500);
h1_time_SiPM01_delta_2->Draw();


MyC->cd(15);

gStyle->SetOptFit(1111);
h1_time_SiPM01_delta_3->Fit("gaus", "","", -200,500);
h1_time_SiPM01_delta_3->Draw();

MyC->cd(16);

gStyle->SetOptFit(1111);
h1_time_SiPM01_delta_4->Fit("gaus", "","", -200,500);
h1_time_SiPM01_delta_4->Draw();


*/

MyC->Update();

//f.Close();

return 0;
}
