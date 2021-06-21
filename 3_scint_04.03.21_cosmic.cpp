// This is an independent project of an individual developer. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com
//-V2008_CYCLOMATIC_COMPLEXITY=100


//#include <iostream.h>
#include "TFile.h"											
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <fstream>
#include <TObject.h>
//#include <TThread.h>

#include "ROOT/TFuture.hxx"

#include <stdio.h>
#include <omp.h>





const Int_t N_number = 15600000; // 4884983 число событий из файла

const Int_t Smooth_index = 2;	// количество прогонов сглаживания гистограммы

//const Double_t L_scint = 100.0 + 2.0 * (6.0 + 12.0); // размер сцинтиллятора + световоды с двух сторон, см, счетчики 120x120 мм

const Double_t L_scint = 106.0 + 2.0 * (13.0 + 12.0); // размер сцинтиллятора + световоды с двух сторон, см, счетчики 200x200 мм


const Double_t k_TDC = 0.1; // 1 канал ЗЦП=0,1 нс

const Double_t b1_1_pedestal = 44.0;	   // пьедисталы для вычета из амплитуд
const Double_t b1_2_pedestal = 34.0;
const Double_t b2_1_pedestal = 54.0;
const Double_t b2_2_pedestal = 68.0;
const Double_t b3_1_pedestal = 64.0;
const Double_t b3_2_pedestal = 58.0;

// диапозон по времени для фита

const Double_t t1_1_left = 1300.0;
const Double_t t1_1_right = 1600.0;

const Double_t t1_2_left = 1300.0;
const Double_t t1_2_right = 1600.0;

const Double_t t2_1_left = 1450.0;
const Double_t t2_1_right = 1650.0;

const Double_t t2_2_left = 1350.0;
const Double_t t2_2_right = 1500.0;

const Double_t t3_1_left = 1300.0;
const Double_t t3_1_right = 1500.0;

const Double_t t3_2_left = 1450.0;
const Double_t t3_2_right = 1700.0;



//// время-амплитудная поправка

///
const Double_t par0_tb1_1_left = 15000.0;
const Double_t par1_tb1_1_left = 1400.0 + 90.0 * 0;
const Double_t par2_tb1_1_left = b1_1_pedestal;

const Double_t par0_tb1_1_right = 15000.0;
const Double_t par1_tb1_1_right = 1400.0 + 90.0 * 0;
const Double_t par2_tb1_1_right = b1_1_pedestal;


const Double_t par0_tb1_2_left = 15000.0;
const Double_t par1_tb1_2_left = 1400.0 + 50.0 * 0;
const Double_t par2_tb1_2_left = b1_2_pedestal;

const Double_t par0_tb1_2_right = 15000.0;
const Double_t par1_tb1_2_right = 1400.0 + 50.0 * 0;
const Double_t par2_tb1_2_right = b1_2_pedestal;

///
const Double_t par0_tb2_1_left = 15000.0;
const Double_t par1_tb2_1_left = 1400.0 + 90.0 * 0;
const Double_t par2_tb2_1_left = b2_1_pedestal;

const Double_t par0_tb2_1_right = 15000.0;
const Double_t par1_tb2_1_right = 1400.0 + 90.0 * 0;
const Double_t par2_tb2_1_right = b2_1_pedestal;


const Double_t par0_tb2_2_left = 15000.0;
const Double_t par1_tb2_2_left = 1400.0 + 90.0 * 0;
const Double_t par2_tb2_2_left = b2_2_pedestal;

const Double_t par0_tb2_2_right = 15000.0;
const Double_t par1_tb2_2_right = 1400.0 + 90.0 * 0;
const Double_t par2_tb2_2_right = b2_2_pedestal;

///
const Double_t par0_tb3_1_left = 15000.0;
const Double_t par1_tb3_1_left = 1400.0 + 50.0 * 0;
const Double_t par2_tb3_1_left = b3_1_pedestal;

const Double_t par0_tb3_1_right = 15000.0;
const Double_t par1_tb3_1_right = 1400.0 + 50.0 * 0;
const Double_t par2_tb3_1_right = b3_1_pedestal;


const Double_t par0_tb3_2_left = 15000.0;
const Double_t par1_tb3_2_left = 1400.0 + 20.0 * 0;
const Double_t par2_tb3_2_left = b3_2_pedestal;

const Double_t par0_tb3_2_right = 15000.0;
const Double_t par1_tb3_2_right = 1400.0 + 20.0 * 0;
const Double_t par2_tb3_2_right = b3_2_pedestal;


/////
/////
///// Каты
/////
/////



const TCut cut_bb_1 = Form("sqrt((b1_1-%f)*(b1_2-%f))>300.0&&sqrt((b1_1-%f)*(b1_2-%f))<1900.0", 1.0 * b1_1_pedestal, 1.0 * b1_2_pedestal, 1.0 * b1_1_pedestal, 1.0 * b1_2_pedestal);
const TCut cut_bb_2 = Form("sqrt((b2_1-%f)*(b2_2-%f))>440.0&&sqrt((b2_1-%f)*(b2_2-%f))<1900.0", 1.0 * b2_1_pedestal, 1.0 * b2_2_pedestal, 1.0 * b2_1_pedestal, 1.0 * b2_2_pedestal);
const TCut cut_bb_3 = Form("sqrt((b3_1-%f)*(b3_2-%f))>390.0&&sqrt((b3_1-%f)*(b3_2-%f))<1900.0", 1.0 * b3_1_pedestal, 1.0 * b3_2_pedestal, 1.0 * b3_1_pedestal, 1.0 * b3_2_pedestal);

//const TCut cut_bb_1 = Form("sqrt((b1_1-%f)*(b1_2-%f))>250.0&&sqrt((b1_1-%f)*(b1_2-%f))<400.0", 1.0 * b1_1_pedestal, 1.0 * b1_2_pedestal, 1.0 * b1_1_pedestal, 1.0 * b1_2_pedestal);
//const TCut cut_bb_2 = Form("sqrt((b2_1-%f)*(b2_2-%f))>500.0&&sqrt((b2_1-%f)*(b2_2-%f))<750.0", 1.0 * b2_1_pedestal, 1.0 * b2_2_pedestal, 1.0 * b2_1_pedestal, 1.0 * b2_2_pedestal);
//const TCut cut_bb_3 = Form("sqrt((b3_1-%f)*(b3_2-%f))>250.0&&sqrt((b3_1-%f)*(b3_2-%f))<400.0", 1.0 * b3_1_pedestal, 1.0 * b3_2_pedestal, 1.0 * b3_1_pedestal, 1.0 * b3_2_pedestal);

//TCut cut_bb_1 = "";
//TCut cut_bb_2 = "";
//TCut cut_bb_3 = "";

//TCut cut_time_1 = "t1_1>800.0&&t1_2>800.0&&t1_1<1200.0&&t1_2<1200.0";
//TCut cut_time_2 = "t2_1>800.0&&t2_2>800.0&&t2_1<1200.0&&t2_2<1200.0";
//TCut cut_time_3 = "t3_1>800.0&&t3_2>800.0&&t3_1<1200.0&&t3_2<1200.0";

const TCut cut_time_1 = "t1_1>900.0&&t1_2>900.0&&t1_1<1400.0&&t1_2<1400.0";
const TCut cut_time_2 = "t2_1>1000.0&&t2_2>1000.0&&t2_1<1500.0&&t2_2<1500.0";
//const TCut cut_time_3 = "t3_1>1200.0&&t3_2>1550.0&&t3_1<2000.0&&t3_2<2000.0";

const TCut cut_time_3 = "t3_2<1480";

//TCut cut_time_1 = "t1_1>1040.0&&t1_2>1070.0&&t1_1<1070.0&&t1_2<1100.0";
//TCut cut_time_2 = "t2_1>1070.0&&t2_2>1070.0&&t2_1<1120.0&&t2_2<1130.0";
//TCut cut_time_3 = "t3_1>1010.0&&t3_2>1020.0&&t3_1<1050.0&&t3_2<1060.0";


// const TCut cut_time_3 = "";


const TCut cut_aplitude_1 = "b1_1>10.0&&b1_2>10.0&&b1_1<2000.0&&b1_2<2000.0";
const TCut cut_aplitude_2 = "b2_1>10.0&&b2_2>10.0&&b2_1<2000.0&&b2_2<2000.0";
const TCut cut_aplitude_3 = "b3_1>10.0&&b3_2>10.0&&b3_1<2000.0&&b3_2<2000.0";
const TCut cut_bb_pedestal = Form("b1_1>%f && b1_2>%f && b2_1>%f && b2_2>%f && b3_1>%f && b3_1>%f", b1_1_pedestal, b1_2_pedestal, b2_1_pedestal, b2_2_pedestal, b3_1_pedestal, b3_2_pedestal);



// const TCut cut_full = cut_aplitude_1 && cut_aplitude_2 && cut_aplitude_3 && cut_time_1 && cut_time_2 && cut_time_3 && cut_bb_1 && cut_bb_2 && cut_bb_3 && cut_bb_pedestal;

//	TCut cut_full = cut_aplitude_1 && cut_aplitude_2 && cut_aplitude_3 && cut_bb_1 && cut_bb_2 && cut_bb_3 && cut_bb_pedestal;

const TCut cut_full = cut_bb_1 && cut_bb_2 && cut_bb_3 && cut_time_3;

//	TCut cut_full = cut_bb_1 && cut_bb_2;

//	TCut cut_full = cut_bb_3;

//const TCut cut_b1b2 = "b2_1>80.0&&b2_2>110.0";

//const TCut cut_full = cut_bb_1&&cut_bb_3&&cut_b1b2;



Double_t S1 = 0.0;
Double_t S2 = 0.0;
Double_t X1 = 0.0;
Double_t X2 = 0.0;





// П-образная функция для фита, 5 параматров
Double_t fit_ln(Double_t* x, Double_t* par) {

	Double_t arg1 = par[1] * x[0] + par[2];
	Double_t arg2 = -par[3] * x[0] + par[4];

	Double_t fitval = par[0] / ((1.0 + TMath::Exp(arg1)) * (1.0 + TMath::Exp(arg2)));

	return fitval;
}


// гипербола для фита, 3 параметра
Double_t fit_giperbola(Double_t* x, Double_t* par) {

	
	Double_t fitval = par[0] / (x[0] - par[1]) + par[2] ;

	return fitval;
}




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


// временное разрешение
Double_t sigma_t(Double_t sigma_31, Double_t sigma_32, Double_t sigma_tr) {

	Double_t temp2 = (pow(sigma_31, 2.0) + pow(sigma_32, 2.0) - pow(sigma_tr, 2.0)) / 2.0;

	Double_t temp = sqrt(temp2);

	printf("Sigma_31 =%.3f\n", sigma_31);
	printf("Sigma_32 =%.3f\n", sigma_32);
	printf("Sigma_tr =%.3f\n", sigma_tr);
	
	if (temp2 > 0) {
		printf("Sigma t=%.3f\n", temp);
	}
	else {

	  printf("Sigma t=i%.3f\n", sqrt(fabs(temp2)));
	}
	
	

	return temp;

}


Double_t delta_t(Double_t x, Double_t par0, Double_t par1, Double_t par2) {

	return par0 / (x - par2) + par1;

}





Double_t bt_left_1_1(Double_t* x, Double_t* par) {

	Double_t xx = x[0];

	return par[0] / (xx - par[1]) + par[2];

}

Double_t bt_right_1_1(Double_t* x, Double_t* par) {

	Double_t xx = x[0];

	return par[0] / (xx - par[1]) + par[2];

}

Double_t bt_left_1_2(Double_t* x, Double_t* par) {

	Double_t xx = x[0];

	return par[0] / (xx - par[1]) + par[2];

}

Double_t bt_right_1_2(Double_t* x, Double_t* par) {

	Double_t xx = x[0];

	return par[0] / (xx - par[1]) + par[2];

}

///
/// 
/// 

Double_t bt_left_2_1(Double_t* x, Double_t* par) {

	Double_t xx = x[0];

	return par[0] / (xx - par[1]) + par[2];

}

Double_t bt_right_2_1(Double_t* x, Double_t* par) {

	Double_t xx = x[0];

	return par[0] / (xx - par[1]) + par[2];

}

Double_t bt_left_2_2(Double_t* x, Double_t* par) {

	Double_t xx = x[0];

	return par[0] / (xx - par[1]) + par[2];

}

Double_t bt_right_2_2(Double_t* x, Double_t* par) {

	Double_t xx = x[0];

	return par[0] / (xx - par[1]) + par[2];

}

///
/// 
/// 

Double_t bt_left_3_1(Double_t* x, Double_t* par) {

	Double_t xx = x[0];

	return par[0] / (xx - par[1]) + par[2];

}

Double_t bt_right_3_1(Double_t* x, Double_t* par) {

	Double_t xx = x[0];

	return par[0] / (xx - par[1]) + par[2];

}

Double_t bt_left_3_2(Double_t* x, Double_t* par) {

	Double_t xx = x[0];

	return par[0] / (xx - par[1]) + par[2];

}

Double_t bt_right_3_2(Double_t* x, Double_t* par) {

	Double_t xx = x[0];

	return par[0] / (xx - par[1]) + par[2];

}

void workItem0()
{
	printf("Running workItem0...\n");
}



TF1* func_1 = new TF1("my_fit_1", fit_ln, -5, 5, 5);
TF1* func_2 = new TF1("my_fit_2", fit_ln, -250, 250, 5);
TF1* func_3 = new TF1("my_fit_3", fit_giperbola, 800, 1500, 3);

// ТЕСТИРОВАНИЕ 8 СЧЕТЧИКОВ 120x120 мм

// TFile* f = new TFile("../macros/stend_04.03.21_3xscintil.root");

// TFile *f = new TFile("../macros/stend_10.03.21_.root");

// TFile* f = TFile::Open("../macros/stend_10.03.21_.root");

// TFile* f = TFile::Open("../macros/stend_12.03.21.root");  // без смазки, счетчик 3 в центре	   

// TFile* f = TFile::Open("../macros/stend_17.03.21.root"); // со смазкой (мало событий)

// TFile* f = TFile::Open("../macros/stend_18.03.21.root"); // со смазкой, счетчик 3 в центре	 1 - сверху, 2 - снизу

// TFile* f = TFile::Open("../macros/stend_19.03.21.root"); // со смазкой, счетчик 1 в центре	 3 - сверху, 2 - снизу

// TFile* f = TFile::Open("../macros/stend_23.03.21.root"); // со смазкой, счетчик 2 в центре	3 - сверху, 1 - снизу

// TFile* f = TFile::Open("../macros/stend_25.03.21.root"); // со смазкой, счетчик 4 в центре	    3 - сверху, 1 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ

// TFile* f = TFile::Open("../macros/stend_26.03.21.root"); // схема с маленьким запускающим счетчиком (около правого ФЭУ, зона 1). со смазкой, счетчик 4 в центре	    3 - сверху, 1 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ

// TFile* f = TFile::Open("../macros/stend_29.03.21.root"); // схема с маленьким запускающим счетчиком (по центру, зона 2). со смазкой, счетчик 4 в центре	    3 - сверху, 1 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ

// TFile* f = TFile::Open("../macros/stend_31.03.21.root"); // схема с маленьким запускающим счетчиком (по центру, зона 2). со смазкой, счетчик 4 в центре	    3 - сверху, 1 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс)

// TFile* f = TFile::Open("../macros/stend_31.03.21_.root"); // схема с маленьким запускающим счетчиком (около правого ФЭУ, зона 1). со смазкой, счетчик 4 в центре	    3 - сверху, 1 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс)

// TFile* f = TFile::Open("../macros/stend_01.04.21.root"); // со смазкой, счетчик 4 в центре	    3 - сверху, 1 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс), вариант триггера №2

// TFile* f = TFile::Open("../macros/stend_02.04.21.root"); // со смазкой, счетчик 5 в центре	    3 - сверху, 1 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс), вариант триггера №3


// TFile* f = TFile::Open("../macros/stend_02.04.21__.root"); // со смазкой, счетчик 5 в центре	    3 - сверху, 1 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс), запуск от ФЭУ 3_1

// TFile* f = TFile::Open("../macros/stend_06.04.21.root"); // со смазкой, счетчик 7 в центре	    3 - сверху, 1 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс), вариант триггера №2

TFile* f = TFile::Open("../macros/stend_07.04.21.root"); // со смазкой, счетчик 6 в центре	    3 - сверху, 1 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс), вариант триггера №2

// TFile* f = TFile::Open("../macros/stend_08.04.21.root"); // со смазкой, счетчик 2 в центре	    3 - сверху, 1 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс), вариант триггера №2

// TFile* f = TFile::Open("../macros/stend_09.04.21.root"); // со смазкой, счетчик 3 в центре	    2 - сверху, 1 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс), вариант триггера №2

// TFile* f = TFile::Open("../macros/stend_12.04.21.root"); // со смазкой, счетчик 8 в центре	    3 - сверху, 1 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс), вариант триггера №2

// TFile* f = TFile::Open("../macros/stend_13.04.21.root"); // со смазкой, счетчик 1 в центре	    3 - сверху, 8 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс), вариант триггера №2


// ТЕСТИРОВАНИЕ 6 СЧЕТЧИКОВ 200x200 мм

// TFile* f = TFile::Open("../macros/stend_20.04.21.root"); // (каналы по амплитуде 1_2 и 3_1 перепутаны местами при наборе, ВЫЗВАТЬ Rename()) со смазкой, счетчик 4 в центре	    5 - сверху, 2 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс), вариант триггера №2


// TFile* f = TFile::Open("../macros/stend_21.04.21.root"); // со смазкой, счетчик 5 в центре	    4 - сверху, 2 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс), вариант триггера №2

// TFile* f = TFile::Open("../macros/stend_22.04.21.root"); // со смазкой, счетчик 2 в центре	    4 - сверху, 5 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс), вариант триггера №2

// TFile* f = TFile::Open("../macros/stend_26.04.21.root"); // со смазкой, счетчик 3 в центре	    1 - сверху, 6 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс), вариант триггера №2

// TFile* f = TFile::Open("../macros/stend_27.04.21.root"); // со смазкой, счетчик 1 в центре	    3 - сверху, 6 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс), вариант триггера №2

// TFile* f = TFile::Open("../macros/stend_29.04.21.root"); // со смазкой, счетчик 6 в центре	    3 - сверху, 1 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс), вариант триггера №2


// ТЕСТИРОВАНИЕ 6 СЧЕТЧИКОВ 200x200 мм + калориметр для измерения по зонам

//TFile* f = TFile::Open("../macros/tqmu_18.05.21.root"); // со смазкой, счетчик №1 счетчик 1, №2 счетчик 3, №3 счетчик 6   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс)





TTree* ntpl = (TTree*)f->Get("ntpl");

//TTree* ntpl = (TTree*)f->Get("T");


//TString[]



TTree* ntpl_select = ntpl->CopyTree(cut_full, "", N_number);	  // создает новое дерево с условиями отбора

//TTree* ntpl_select = 0 ;	  // создает новое дерево с условиями отбора

//	TTree* ntpl_select = (TTree*)ntpl_select->Get("ntpl_select");


/*
  void AddBranch() {

	Float_t a1_1, a1_2, a2_1, a2_2, a3_1, a3_2;
	Float_t t1_1, t1_2, t2_1, t2_2, t3_1, t3_2;


	auto a1_1 = ntpl->Branch("a1_1", &a1_1, "a1_1/F");
	auto a1_2 = ntpl->Branch("a1_2", &a1_2, "a1_2/F");
	auto a2_1 = ntpl->Branch("a2_1", &a2_1, "a2_1/F");
	auto a2_2 = ntpl->Branch("a2_2", &a2_2, "a2_2/F");
	auto a3_1 = ntpl->Branch("a3_1", &a3_1, "a3_1/F");
	auto a3_2 = ntpl->Branch("a3_2", &a3_2, "a3_2/F");

	auto t1_1 = ntpl->Branch("t1_1", &t1_1, "t1_1/F");
	auto t1_2 = ntpl->Branch("t1_2", &t1_2, "t1_2/F");
	auto t2_1 = ntpl->Branch("t2_1", &t2_1, "t2_1/F");
	auto t2_2 = ntpl->Branch("t2_2", &t2_2, "t2_2/F");
	auto t3_1 = ntpl->Branch("t3_1", &t3_1, "t3_1/F");
	auto t3_2 = ntpl->Branch("t3_2", &t3_2, "t3_2/F");

	Long64_t nentries = ntpl->GetEntries();    // Read the number of entries in the ntpl.
	for (Long64_t i = 0; i < nentries; i++) {
		
		a1_1 = Copy(aF[0]);
		
		a1_1->Fill();
		a1_2->Fill();
		a2_1->Fill();
		a2_2->Fill();
		a3_1->Fill();
		a3_2->Fill();

		t1_1->Fill();
		t1_2->Fill();
		t2_1->Fill();
		t2_2->Fill();
		t3_1->Fill();
		t3_2->Fill();
	}
	ntpl->Write("", TObject::kOverwrite);       // Save only the new version of the tree.
}

*/


void Rename() {

	ntpl->GetBranch("b1_2")->SetTitle("b1_2_");
	ntpl->GetBranch("b1_2")->SetName("b1_2_");
	ntpl->GetLeaf("b1_2")->SetTitle("b1_2_");
	ntpl->GetLeaf("b1_2")->SetName("b1_2_");

	ntpl->GetBranch("b3_1")->SetTitle("b1_2");
	ntpl->GetBranch("b3_1")->SetName("b1_2");
	ntpl->GetLeaf("b3_1")->SetTitle("b1_2");
	ntpl->GetLeaf("b3_1")->SetName("b1_2");


	ntpl->GetBranch("b1_2_")->SetTitle("b3_1");
	ntpl->GetBranch("b1_2_")->SetName("b3_1");
	ntpl->GetLeaf("b1_2_")->SetTitle("b3_1");
	ntpl->GetLeaf("b1_2_")->SetName("b3_1");

	//ntpl->Write();

	ntpl_select = ntpl->CopyTree(cut_full, "", N_number);	  // создает новое дерево с условиями отбора

	/*
	
	

	ntpl_select->GetBranch("b1_2")->SetTitle("b1_2_");
	ntpl_select->GetBranch("b1_2")->SetName("b1_2_");
	ntpl_select->GetLeaf("b1_2")->SetTitle("b1_2_");
	ntpl_select->GetLeaf("b1_2")->SetName("b1_2_");

	ntpl_select->GetBranch("b3_1")->SetTitle("b1_2");
	ntpl_select->GetBranch("b3_1")->SetName("b1_2");
	ntpl_select->GetLeaf("b3_1")->SetTitle("b1_2");
	ntpl_select->GetLeaf("b3_1")->SetName("b1_2");


	ntpl_select->GetBranch("b1_2_")->SetTitle("b3_1");
	ntpl_select->GetBranch("b1_2_")->SetName("b3_1");
	ntpl_select->GetLeaf("b1_2_")->SetTitle("b3_1");
	ntpl_select->GetLeaf("b1_2_")->SetName("b3_1");
	*/
}

// временное разрешение
void MyC_1() {

	
	TH1F* h1_time = new TH1F("h1_time", "((t3_1-t3_2)-0.5*((t1_1-t1_2)-(t2_1-t2_2)))/sqrt(3) WITH AMPLITUDE CORRECTION", 200, -300, 300);

	TH2* h_time_amplitude = new TH2D("h_time_amplitude", "0.5*(t3_1-t3_2):0.5*ln(A3_1/A3_2)+Cut WITH AMPLITUDE CORRECTION", 100, -5.0, 5.0, 1000, -500.0, 500.0);
	
	TH1F* h1_t1t2_1 = new TH1F("h1_t1t2_1", "0.5*(t1_1+t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	TH1F* h2_t1t2_1 = new TH1F("h2_t1t2_1", "0.5*(t2_1+t2_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	TH1F* h3_t1t2_1 = new TH1F("h3_t1t2_1", "0.5*(t3_1+t3_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);

	TH1F* h1_t1t2_2 = new TH1F("h1_t1t2_2", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	TH1F* h2_t1t2_2 = new TH1F("h2_t1t2_2", "0.5*(t2_1-t2_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	TH1F* h3_t1t2_2 = new TH1F("h3_t1t2_2", "0.5*(t3_1-t3_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);

	
	
	TCanvas* MyC_1 = new TCanvas("MyC_1", "MyC_1", 0, 0, 2000, 1000);	  // для определения числа ф.э.

	MyC_1->Divide(3, 3);


	MyC_1->cd(1);

	TString vartime = Form("((t3_1-delta_t(b3_1, %f, %f, %f))-(t3_2-delta_t(b3_2, %f, %f, %f))-0.5*(((t1_1-delta_t(b1_1, %f, %f, %f))-(t1_2-delta_t(b1_2, %f, %f, %f)))-((t2_1-delta_t(b2_1, %f, %f, %f))-(t2_2-delta_t(b2_2, %f, %f, %f)))))/(%f)>>h1_time", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right, par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right, sqrt(3));


	ntpl_select->Draw(vartime);

	h1_time->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h1_time->Fit("gaus", "", "", -200, 200);



	h1_time->GetXaxis()->SetTitle("Time, channel number");
	h1_time->GetYaxis()->SetTitle("Events number");
	h1_time->GetXaxis()->CenterTitle();
	h1_time->GetYaxis()->CenterTitle();

	h1_time->Draw();


	MyC_1->cd(2);



	TString vartb = Form("0.5*((t3_1-delta_t(b3_1, %f, %f, %f))-(t3_2-delta_t(b3_2, %f, %f, %f))):0.5*log((b3_1-%f)/(b3_2-%f))>>h_time_amplitude", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, b3_1_pedestal, b3_2_pedestal);

	ntpl_select->Draw(vartb);


	h_time_amplitude->GetXaxis()->SetTitle("Ln(A3_1/A3_2)/2");
	h_time_amplitude->GetYaxis()->SetTitle("Time (T3_1-T3_2)/2, channel number");
	h_time_amplitude->GetXaxis()->CenterTitle();
	h_time_amplitude->GetYaxis()->CenterTitle();

	TF1* f1 = new TF1("f1", "[0]+[1]*x", -5, 5);
	f1->SetParameters(0., 1.);
	f1->SetLineColor(kRed);
	h_time_amplitude->Fit(f1);
	h_time_amplitude->Draw();
	f1->Draw("same");


	MyC_1->cd(3);




	MyC_1->cd(4);



	TString vart1_1 = Form("0.5*((t1_1-delta_t(b1_1, %f, %f, %f))+(t1_2-delta_t(b1_2, %f, %f, %f)))>>h1_t1t2_1", par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right);


	ntpl_select->Draw(vart1_1);

	h1_t1t2_1->Smooth(Smooth_index);	 // Сглаживание

	h1_t1t2_1->Fit("gaus", "", "", -200, 200);

	h1_t1t2_1->GetXaxis()->SetTitle("Time, channel number");
	h1_t1t2_1->GetYaxis()->SetTitle("Events number");
	h1_t1t2_1->GetXaxis()->CenterTitle();
	h1_t1t2_1->GetYaxis()->CenterTitle();


	h1_t1t2_1->Draw();

	MyC_1->cd(5);



	TString vart2_1 = Form("0.5*((t2_1-delta_t(b2_1, %f, %f, %f))+(t2_2-delta_t(b2_2, %f, %f, %f)))>>h2_t1t2_1", par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);


	ntpl_select->Draw(vart2_1);

	h2_t1t2_1->Smooth(Smooth_index);	 // Сглаживание

	h2_t1t2_1->Fit("gaus", "", "", -200, 200);

	h2_t1t2_1->GetXaxis()->SetTitle("Time, channel number");
	h2_t1t2_1->GetYaxis()->SetTitle("Events number");
	h2_t1t2_1->GetXaxis()->CenterTitle();
	h2_t1t2_1->GetYaxis()->CenterTitle();


	h2_t1t2_1->Draw();

	MyC_1->cd(6);



	TString vart3_1 = Form("0.5*((t3_1-delta_t(b3_1, %f, %f, %f))+(t3_2-delta_t(b3_2, %f, %f, %f)))>>h3_t1t2_1", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right);


	ntpl_select->Draw(vart3_1);

	h3_t1t2_1->Smooth(Smooth_index);	 // Сглаживание
	
	h3_t1t2_1->Fit("gaus", "", "", -200, 200);

	h3_t1t2_1->GetXaxis()->SetTitle("Time, channel number");
	h3_t1t2_1->GetYaxis()->SetTitle("Events number");
	h3_t1t2_1->GetXaxis()->CenterTitle();
	h3_t1t2_1->GetYaxis()->CenterTitle();


	h3_t1t2_1->Draw();

	MyC_1->cd(7);


	TString vart1_2 = Form("0.5*((t1_1-delta_t(b1_1, %f, %f, %f))-(t1_2-delta_t(b1_2, %f, %f, %f)))>>h1_t1t2_2", par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right);

	h1_t1t2_2->Smooth(Smooth_index);	 // Сглаживание
	
	ntpl_select->Draw(vart1_2);



	h1_t1t2_2->GetXaxis()->SetTitle("Time, channel number");
	h1_t1t2_2->GetYaxis()->SetTitle("Events number");
	h1_t1t2_2->GetXaxis()->CenterTitle();
	h1_t1t2_2->GetYaxis()->CenterTitle();

	gStyle->SetOptFit(1111);
	h1_t1t2_2->Fit("gaus", "", "", -150, 150);
	h1_t1t2_2->Fit("my_fit_2", "+", "", -150, 150);


	h1_t1t2_2->Draw();

	S1 = func_2->GetParameter(1);
	X1 = func_2->GetParameter(2);
	S2 = func_2->GetParameter(3);
	X2 = func_2->GetParameter(4);


	printf("v_scint= %.3f ns/m\n", k_TDC * delta_L(S1, X1, S2, X2) / (0.01 * L_scint));

	MyC_1->cd(8);



	TString vart2_2 = Form("0.5*((t2_1-delta_t(b2_1, %f, %f, %f))-(t2_2-delta_t(b2_2, %f, %f, %f)))>>h2_t1t2_2", par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);

	h2_t1t2_2->Smooth(Smooth_index);	 // Сглаживание

	ntpl_select->Draw(vart2_2);



	h2_t1t2_2->GetXaxis()->SetTitle("Time, channel number");
	h2_t1t2_2->GetYaxis()->SetTitle("Events number");
	h2_t1t2_2->GetXaxis()->CenterTitle();
	h2_t1t2_2->GetYaxis()->CenterTitle();

	gStyle->SetOptFit(1111);
	h2_t1t2_2->Fit("gaus", "", "", -150, 150);
	h2_t1t2_2->Fit("my_fit_2", "+", "", -150, 150);
	h2_t1t2_2->Draw();

	S1 = func_2->GetParameter(1);
	X1 = func_2->GetParameter(2);
	S2 = func_2->GetParameter(3);
	X2 = func_2->GetParameter(4);

	
	printf("v_scint= %.3f ns/m\n", k_TDC * delta_L(S1, X1, S2, X2) / (0.01 * L_scint));

	MyC_1->cd(9);



	TString vart3_2 = Form("0.5*((t3_1-delta_t(b3_1, %f, %f, %f))-(t3_2-delta_t(b3_2,  %f, %f, %f)))>>h3_t1t2_2", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right);

	h3_t1t2_2->Smooth(Smooth_index);	 // Сглаживание

	ntpl_select->Draw(vart3_2);

	h3_t1t2_2->GetXaxis()->SetTitle("Time, channel number");
	h3_t1t2_2->GetYaxis()->SetTitle("Events number");
	h3_t1t2_2->GetXaxis()->CenterTitle();
	h3_t1t2_2->GetYaxis()->CenterTitle();

	gStyle->SetOptFit(1111);
	h3_t1t2_2->Fit("gaus", "", "", -150, 150);
	h3_t1t2_2->Fit("my_fit_2", "+", "", -150, 150);
	h3_t1t2_2->Draw();

	S1 = func_2->GetParameter(1);
	X1 = func_2->GetParameter(2);
	S2 = func_2->GetParameter(3);
	X2 = func_2->GetParameter(4);

	printf("v_scint= %.3f ns/m\n", k_TDC * delta_L(S1, X1, S2, X2) / (0.01 * L_scint));

	MyC_1->Update();

}

// амплитуда-амплитуда
void MyC_2() {

	TH1F* h1_b1b2_1 = new TH1F("h1_b1b2_1", "0.5*ln(A1_1/A1_2) -pedestal", 200, -2.0, 2.0);
	TH1F* h2_b1b2_1 = new TH1F("h2_b1b2_1", "0.5*ln(A2_1/A2_2) -pedestal", 200, -2.0, 2.0);
	TH1F* h3_b1b2_1 = new TH1F("h3_b1b2_1", "0.5*ln(A3_1/A3_2) -pedestal", 200, -2.0, 2.0);

	TH1F* h1_b1b2_2 = new TH1F("h1_b1b2_2", "sqrt(A1_1*A1_2) -pedestal", 1500, 0, 1500.0);
	TH1F* h2_b1b2_2 = new TH1F("h2_b1b2_2", "sqrt(A2_1*A2_2) -pedestal", 1500, 0, 1500.0);
	TH1F* h3_b1b2_2 = new TH1F("h3_b1b2_2", "sqrt(A3_1*A3_2) -pedestal", 1500, 0, 1500.0);
	
	
	TH2* h4_1 = new TH2D("h4_1", "b1_1:b1_2 -pedestal", 2100, 0, 2100, 2100, 0, 2100);
	TH2* h5_1 = new TH2D("h5_1", "b2_1:b2_2 -pedestal", 2100, 0, 2100, 2100, 0, 2100);
	TH2* h6_1 = new TH2D("h6_1", "b3_1:b3_2 -pedestal", 2100, 0, 2100, 2100, 0, 2100);

	TH2* h4_1_same = new TH2D("h4_1_same", "b1_1:b1_2+Cut -pedestal", 2100, 0, 2100, 2100, 0, 2100);
	TH2* h5_1_same = new TH2D("h5_1_same", "b2_1:b2_2+Cut -pedestal", 2100, 0, 2100, 2100, 0, 2100);
	TH2* h6_1_same = new TH2D("h6_1_same", "b3_1:b3_2+Cut -pedestal", 2100, 0, 2100, 2100, 0, 2100);
	
	
	// set the parameters to the mean and RMS of the histogram
	func_1->SetParameters(10.0, 5.0, -5.0, 5.0, -5.0);
	func_2->SetParameters(10.0, 1.0, -100.0, 1.0, -100.0);
	// give the parameters meaningful names
	func_1->SetParNames("A", "S1", "X1", "S2", "X2");
	func_2->SetParNames("A", "S1", "X1", "S2", "X2");


	func_3->SetParameters(10000.0, 1000.0, 10.0);
	func_3->SetParNames("k", "t0", "b0");

	
	TCanvas* MyC_2 = new TCanvas("MyC_2", "MyC_2", 0, 0, 2000, 1000);     


	MyC_2->Divide(3, 3);

	MyC_2->cd(1);


	TString var1 = Form("(b1_1 - %f) : (b1_2 - %f) >> h4_1", b1_1_pedestal, b1_2_pedestal);
	TString var1_same = Form("(b1_1 - %f) : (b1_2 - %f) >> h4_1_same", b1_1_pedestal, b1_2_pedestal);

	ntpl->Draw(var1, "", "", N_number);
	ntpl_select->Draw(var1_same);

	h4_1->GetXaxis()->SetTitle("Amplitude, channel number");
	h4_1->GetYaxis()->SetTitle("Amplitude, channel number");
	h4_1->GetXaxis()->CenterTitle();
	h4_1->GetYaxis()->CenterTitle();

	
	h4_1->Draw();
	h4_1_same->SetMarkerColor(kRed);
	h4_1_same->Draw("same");

	MyC_2->cd(2);

	TString var2 = Form("(b2_1 - %f) : (b2_2 - %f) >> h5_1", b2_1_pedestal, b2_2_pedestal);
	TString var2_same = Form("(b2_1 - %f) : (b2_2 - %f) >> h5_1_same", b2_1_pedestal, b2_2_pedestal);

	
	
	ntpl->Draw(var2, "", "", N_number);
	ntpl_select->Draw(var2_same);

	h5_1->GetXaxis()->SetTitle("Amplitude, channel number");
	h5_1->GetYaxis()->SetTitle("Amplitude, channel number");
	h5_1->GetXaxis()->CenterTitle();
	h5_1->GetYaxis()->CenterTitle();

	h5_1->Draw();
	h5_1_same->SetMarkerColor(kRed);
	h5_1_same->Draw("same");


	MyC_2->cd(3);

	TString var3 = Form("(b3_1 - %f) : (b3_2 - %f) >> h6_1", b3_1_pedestal, b3_2_pedestal);
	TString var3_same = Form("(b3_1 - %f) : (b3_2 - %f) >> h6_1_same", b3_1_pedestal, b3_2_pedestal);


	ntpl->Draw(var3, "", "", N_number);
	ntpl_select->Draw(var3_same);


	h6_1->GetXaxis()->SetTitle("Amplitude, channel number");
	h6_1->GetYaxis()->SetTitle("Amplitude, channel number");
	h6_1->GetXaxis()->CenterTitle();
	h6_1->GetYaxis()->CenterTitle();

	h6_1->Draw();
	h6_1_same->SetMarkerColor(kRed);
	h6_1_same->SetMarkerStyle(6);
	h6_1_same->Draw("same");

	MyC_2->cd(4);

	TString var4 = Form("sqrt((b1_1 - %f)*(b1_2 - %f)) >> h1_b1b2_2", b1_1_pedestal, b1_2_pedestal);

	h1_b1b2_2->Smooth(Smooth_index);	 // Сглаживание

	
	ntpl_select->Draw(var4);

	gStyle->SetOptFit(1111);
	h1_b1b2_2->Fit("landau", "", "", 100.0, 800.0);

	h1_b1b2_2->SetFillColor(kGreen);
	h1_b1b2_2->Draw();

	MyC_2->cd(5);

	TString var5 = Form("sqrt((b2_1 - %f)*(b2_2 - %f)) >> h2_b1b2_2", b2_1_pedestal, b2_2_pedestal);

	h2_b1b2_2->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(var5);
	gStyle->SetOptFit(1111);
	h2_b1b2_2->Fit("landau", "", "", 100.0, 800.0);

	h2_b1b2_2->SetFillColor(kGreen);
	h2_b1b2_2->Draw();

	MyC_2->cd(6);

	TString var6 = Form("sqrt((b3_1 - %f)*(b3_2 - %f)) >> h3_b1b2_2", b3_1_pedestal, b3_2_pedestal);

	h3_b1b2_2->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(var6);
	gStyle->SetOptFit(1111);
	h3_b1b2_2->Fit("landau", "", "", 100.0, 800.0);

	h3_b1b2_2->SetFillColor(kGreen);
	h3_b1b2_2->Draw();

	MyC_2->cd(7);

	TString var7 = Form("0.5*log((b1_1 - %f)/(b1_2 - %f)) >> h1_b1b2_1", b1_1_pedestal, b1_2_pedestal);

	h1_b1b2_1->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(var7);

	gStyle->SetOptFit(1111);
	h1_b1b2_1->Fit("my_fit_1", "", "", -2, 2);
	h1_b1b2_1->Draw();

	S1 = func_1->GetParameter(1);
	X1 = func_1->GetParameter(2);
	S2 = func_1->GetParameter(3);
	X2 = func_1->GetParameter(4);
	 	
	printf("lambda_scint= %.3f cm\n", L_scint / delta_L(S1, X1, S2, X2));


	MyC_2->cd(8);

	TString var8 = Form("0.5*log((b2_1 - %f)/(b2_2 - %f)) >> h2_b1b2_1", b2_1_pedestal, b2_2_pedestal);

	h2_b1b2_1->Smooth(Smooth_index);	 // Сглаживание

	
	ntpl_select->Draw(var8);

	gStyle->SetOptFit(1111);
	h2_b1b2_1->Fit("my_fit_1", "", "", -2, 2);
	h2_b1b2_1->Draw();

	S1 = func_1->GetParameter(1);
	X1 = func_1->GetParameter(2);
	S2 = func_1->GetParameter(3);
	X2 = func_1->GetParameter(4);
	
	printf("lambda_scint= %.3f cm\n", L_scint / delta_L(S1, X1, S2, X2));

	MyC_2->cd(9);

	TString var9 = Form("0.5*log((b3_1 - %f)/(b3_2 - %f)) >> h3_b1b2_1", b3_1_pedestal, b3_2_pedestal);

	h3_b1b2_1->Smooth(Smooth_index);	 // Сглаживание
 	
	ntpl_select->Draw(var9);

	gStyle->SetOptFit(1111);
	h3_b1b2_1->Fit("my_fit_1", "", "", -2, 2);
	h3_b1b2_1->Draw();

	S1 = func_1->GetParameter(1);
	X1 = func_1->GetParameter(2);
	S2 = func_1->GetParameter(3);
	X2 = func_1->GetParameter(4);
	
	printf("lambda_scint= %.3f cm\n", L_scint / delta_L(S1, X1, S2, X2));



	MyC_2->Update();

}

//время, исходный файл
void MyC_3() {

	TH1F* ht1_1 = new TH1F("ht1_1", "t1_1 without cuts", 2500, 0, 2500);
	TH1F* ht1_2 = new TH1F("ht1_2", "t1_2 without cuts", 2500, 0, 2500);
	TH1F* ht2_1 = new TH1F("ht2_1", "t2_1 without cuts", 2500, 0, 2500);
	TH1F* ht2_2 = new TH1F("ht2_2", "t2_2 without cuts", 2500, 0, 2500);
	TH1F* ht3_1 = new TH1F("ht3_1", "t3_1 without cuts", 2500, 0, 2500);
	TH1F* ht3_2 = new TH1F("ht3_2", "t3_2 without cuts", 2500, 0, 2500);


	TCanvas* MyC_3 = new TCanvas("MyC_3", "MyC_3", 0, 0, 2000, 1000);	  

	MyC_3->Divide(2, 3);

	MyC_3->cd(1);

	ntpl->Draw("t1_1>>ht1_1", "", "", N_number);

	ht1_1->GetXaxis()->SetTitle("Time, channel number");
	ht1_1->GetYaxis()->SetTitle("Events");
	ht1_1->GetXaxis()->CenterTitle();
	ht1_1->GetYaxis()->CenterTitle();

	MyC_3->cd(2);

	ntpl->Draw("t1_2>>ht1_2", "", "", N_number);

	ht1_2->GetXaxis()->SetTitle("Time, channel number");
	ht1_2->GetYaxis()->SetTitle("Events");
	ht1_2->GetXaxis()->CenterTitle();
	ht1_2->GetYaxis()->CenterTitle();

	MyC_3->cd(3);

	ntpl->Draw("t2_1>>ht2_1", "", "", N_number);

	ht2_1->GetXaxis()->SetTitle("Time, channel number");
	ht2_1->GetYaxis()->SetTitle("Events");
	ht2_1->GetXaxis()->CenterTitle();
	ht2_1->GetYaxis()->CenterTitle();


	MyC_3->cd(4);

	ntpl->Draw("t2_2>>ht2_2", "", "", N_number);

	ht2_2->GetXaxis()->SetTitle("Time, channel number");
	ht2_2->GetYaxis()->SetTitle("Events");
	ht2_2->GetXaxis()->CenterTitle();
	ht2_2->GetYaxis()->CenterTitle();

	MyC_3->cd(5);

	ntpl->Draw("t3_1>>ht3_1", "", "", N_number);

	ht3_1->GetXaxis()->SetTitle("Time, channel number");
	ht3_1->GetYaxis()->SetTitle("Events");
	ht3_1->GetXaxis()->CenterTitle();
	ht3_1->GetYaxis()->CenterTitle();


	MyC_3->cd(6);

	ntpl->Draw("t3_2>>ht3_2", "", "", N_number);

	ht3_2->GetXaxis()->SetTitle("Time, channel number");
	ht3_2->GetYaxis()->SetTitle("Events");
	ht3_2->GetXaxis()->CenterTitle();
	ht3_2->GetYaxis()->CenterTitle();


	MyC_3->Update();

}

// амплитуда,исходный файл
void MyC_4() {



	TH1F* hb1_1 = new TH1F("hb1_1", "b1_1 without cuts", 1500, 0, 1500);
	TH1F* hb1_2 = new TH1F("hb1_2", "b1_2 without cuts", 1500, 0, 1500);
	TH1F* hb2_1 = new TH1F("hb2_1", "b2_1 without cuts", 1500, 0, 1500);
	TH1F* hb2_2 = new TH1F("hb2_2", "b2_2 without cuts", 1500, 0, 1500);
	TH1F* hb3_1 = new TH1F("hb3_1", "b3_1 without cuts", 1500, 0, 1500);
	TH1F* hb3_2 = new TH1F("hb3_2", "b3_2 without cuts", 1500, 0, 1500);



	TCanvas* MyC_4 = new TCanvas("MyC_4", "MyC_4", 0, 0, 2000, 1000);	  



	MyC_4->Divide(2, 3);

	MyC_4->cd(1);

	ntpl->Draw("b1_1>>hb1_1", "", "", N_number);

	hb1_1->GetXaxis()->SetTitle("Amplitude, channel number");
	hb1_1->GetYaxis()->SetTitle("Events");
	hb1_1->GetXaxis()->CenterTitle();
	hb1_1->GetYaxis()->CenterTitle();

	MyC_4->cd(2);

	ntpl->Draw("b1_2>>hb1_2", "", "", N_number);

	hb1_2->GetXaxis()->SetTitle("Amplitude, channel number");
	hb1_2->GetYaxis()->SetTitle("Events");
	hb1_2->GetXaxis()->CenterTitle();
	hb1_2->GetYaxis()->CenterTitle();


	MyC_4->cd(3);

	ntpl->Draw("b2_1>>hb2_1", "", "", N_number);

	hb2_1->GetXaxis()->SetTitle("Amplitude, channel number");
	hb2_1->GetYaxis()->SetTitle("Events");
	hb2_1->GetXaxis()->CenterTitle();
	hb2_1->GetYaxis()->CenterTitle();


	MyC_4->cd(4);


	ntpl->Draw("b2_2>>hb2_2", "", "", N_number);

	hb2_2->GetXaxis()->SetTitle("Amplitude, channel number");
	hb2_2->GetYaxis()->SetTitle("Events");
	hb2_2->GetXaxis()->CenterTitle();
	hb2_2->GetYaxis()->CenterTitle();

	MyC_4->cd(5);

	ntpl->Draw("b3_1>>hb3_1", "", "", N_number);

	hb3_1->GetXaxis()->SetTitle("Amplitude, channel number");
	hb3_1->GetYaxis()->SetTitle("Events");
	hb3_1->GetXaxis()->CenterTitle();
	hb3_1->GetYaxis()->CenterTitle();


	MyC_4->cd(6);

	ntpl->Draw("b3_2>>hb3_2", "", "", N_number);

	hb3_2->GetXaxis()->SetTitle("Amplitude, channel number");
	hb3_2->GetYaxis()->SetTitle("Events");
	hb3_2->GetXaxis()->CenterTitle();
	hb3_2->GetYaxis()->CenterTitle();

	MyC_4->Update();

}

// время-амплитуда
void MyC_5() {


	TH2* h1_1 = new TH2D("h1_1", "b1_1:t1_1", 500, 0, 2100, 1500, 0, 1500);
	TH2* h1_2 = new TH2D("h1_2", "b1_2:t1_2", 500, 0, 2100, 1500, 0, 1500);
	TH2* h2_1 = new TH2D("h2_1", "b2_1:t2_1", 500, 0, 2100, 1500, 0, 1500);
	TH2* h2_2 = new TH2D("h2_2", "b2_2:t2_2", 500, 0, 2100, 1500, 0, 1500);
	TH2* h3_1 = new TH2D("h3_1", "b3_1:t3_1", 500, 0, 2100, 1500, 0, 1500);
	TH2* h3_2 = new TH2D("h3_2", "b3_2:t3_2", 500, 0, 2100, 1500, 0, 1500);




	TH2* h1_1_same = new TH2D("h1_1_same", "b1_1:t1_1+Cut", 500, 0, 2100, 1500, 0, 1500);
	TH2* h1_2_same = new TH2D("h1_2_same", "b1_2:t1_2+Cut", 500, 0, 2100, 1500, 0, 1500);
	TH2* h2_1_same = new TH2D("h2_1_same", "b2_1:t2_1+Cut", 500, 0, 2100, 1500, 0, 1500);
	TH2* h2_2_same = new TH2D("h2_2_same", "b2_2:t2_2+Cut", 500, 0, 2100, 1500, 0, 1500);
	TH2* h3_1_same = new TH2D("h3_1_same", "b3_1:t3_1+Cut", 500, 0, 2100, 1500, 0, 1500);
	TH2* h3_2_same = new TH2D("h3_2_same", "b3_2:t3_2+Cut", 500, 0, 2100, 1500, 0, 1500);



	TF1* f_bt_left_1_1 = new TF1("f_bt_left_1_1", bt_left_1_1, 1000, 1500, 3);
	TF1* f_bt_right_1_1 = new TF1("f_bt_right_1_1", bt_right_1_1, 1000, 1500, 3);
	TF1* f_bt_left_1_2 = new TF1("f_bt_left_1_2", bt_left_1_2, 1000, 1500, 3);
	TF1* f_bt_right_1_2 = new TF1("f_bt_right_1_2", bt_right_1_2, 1000, 1500, 3);

	TF1* f_bt_left_2_1 = new TF1("f_bt_left_2_1", bt_left_2_1, 1000, 1500, 3);
	TF1* f_bt_right_2_1 = new TF1("f_bt_right_2_1", bt_right_2_1, 1000, 1500, 3);
	TF1* f_bt_left_2_2 = new TF1("f_bt_left_2_2", bt_left_2_2, 1000, 1500, 3);
	TF1* f_bt_right_2_2 = new TF1("f_bt_right_2_2", bt_right_2_2, 1000, 1500, 3);

	TF1* f_bt_left_3_1 = new TF1("f_bt_left_3_1", bt_left_3_1, 1000, 1500, 3);
	TF1* f_bt_right_3_1 = new TF1("f_bt_right_3_1", bt_right_3_1, 1000, 1500, 3);
	TF1* f_bt_left_3_2 = new TF1("f_bt_left_3_2", bt_left_3_2, 1000, 1500, 3);
	TF1* f_bt_right_3_2 = new TF1("f_bt_right_3_2", bt_right_3_2, 1000, 1500, 3);


	f_bt_left_1_1->SetParNames("coefficient", "mean", "constant");
	f_bt_right_1_1->SetParNames("coefficient", "mean", "constant");
	f_bt_left_1_2->SetParNames("coefficient", "mean", "constant");
	f_bt_right_1_2->SetParNames("coefficient", "mean", "constant");

	f_bt_left_2_1->SetParNames("coefficient", "mean", "constant");
	f_bt_right_2_1->SetParNames("coefficient", "mean", "constant");
	f_bt_left_2_2->SetParNames("coefficient", "mean", "constant");
	f_bt_right_2_2->SetParNames("coefficient", "mean", "constant");

	f_bt_left_3_1->SetParNames("coefficient", "mean", "constant");
	f_bt_right_3_1->SetParNames("coefficient", "mean", "constant");
	f_bt_left_3_2->SetParNames("coefficient", "mean", "constant");
	f_bt_right_3_2->SetParNames("coefficient", "mean", "constant");





	TCanvas* MyC_5 = new TCanvas("MyC_5", "MyC_5", 0, 0, 2000, 1000);	  



	MyC_5->Divide(2, 3);

	MyC_5->cd(1);

	ntpl->Draw("b1_1:t1_1>>h1_1", "", "", N_number);
	ntpl_select->Draw("b1_1:t1_1>>h1_1_same");

	h1_1->GetXaxis()->SetTitle("Time, channel number");
	h1_1->GetYaxis()->SetTitle("Amplitude, channel numbe");
	h1_1->GetXaxis()->CenterTitle();
	h1_1->GetYaxis()->CenterTitle();

	h1_1->Draw();
	h1_1_same->SetMarkerColor(kRed);
	h1_1_same->Draw("same");

	gStyle->SetOptFit(1111);
	h1_1_same->Fit("my_fit_3", "", "", t1_1_left, t1_1_right);


	//f_bt_left_1_1->SetMarkerColor(kGreen);
	f_bt_right_1_1->SetMarkerColor(kGreen);

	//f_bt_left_1_1->SetParameters(par0_tb1_1_left, par1_tb1_1_left, par2_tb1_1_left);
	f_bt_right_1_1->SetParameters(par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right);

	//f_bt_left_1_1->Draw("same");
	f_bt_right_1_1->Draw("same");





	MyC_5->cd(2);

	ntpl->Draw("b1_2:t1_2>>h1_2", "", "", N_number);
	ntpl_select->Draw("b1_2:t1_2>>h1_2_same");

	h1_2->GetXaxis()->SetTitle("Time, channel number");
	h1_2->GetYaxis()->SetTitle("Amplitude, channel numbe");
	h1_2->GetXaxis()->CenterTitle();
	h1_2->GetYaxis()->CenterTitle();

	h1_2->Draw();
	h1_2_same->SetMarkerColor(kRed);
	h1_2_same->Draw("same");

	gStyle->SetOptFit(1111);
	h1_2_same->Fit("my_fit_3", "", "", t1_2_left, t1_2_right);

	//f_bt_left_1_2->SetMarkerColor(kGreen);
	f_bt_right_1_2->SetMarkerColor(kGreen);

	//f_bt_left_1_2->SetParameters(par0_tb1_2_left, par1_tb1_2_left, par2_tb1_2_left);
	f_bt_right_1_2->SetParameters(par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right);

	//f_bt_left_1_2->Draw("same");
	f_bt_right_1_2->Draw("same");

	MyC_5->cd(3);

	ntpl->Draw("b2_1:t2_1>>h2_1", "", "", N_number);
	ntpl_select->Draw("b2_1:t2_1>>h2_1_same");

	h2_1->GetXaxis()->SetTitle("Time, channel number");
	h2_1->GetYaxis()->SetTitle("Amplitude, channel numbe");
	h2_1->GetXaxis()->CenterTitle();
	h2_1->GetYaxis()->CenterTitle();

	h2_1->Draw();
	h2_1_same->SetMarkerColor(kRed);
	h2_1_same->Draw("same");

	gStyle->SetOptFit(1111);
	h2_1_same->Fit("my_fit_3", "", "", t2_1_left, t2_1_right);

	//f_bt_left_2_1->SetMarkerColor(kGreen);
	f_bt_right_2_1->SetMarkerColor(kGreen);

	//f_bt_left_2_1->SetParameters(par0_tb2_1_left, par1_tb2_1_left, par2_tb2_1_left);
	f_bt_right_2_1->SetParameters(par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right);

	//f_bt_left_2_1->Draw("same");
	f_bt_right_2_1->Draw("same");


	MyC_5->cd(4);

	ntpl->Draw("b2_2:t2_2>>h2_2", "", "", N_number);
	ntpl_select->Draw("b2_2:t2_2>>h2_2_same");

	h2_2->GetXaxis()->SetTitle("Time, channel number");
	h2_2->GetYaxis()->SetTitle("Amplitude, channel numbe");
	h2_2->GetXaxis()->CenterTitle();
	h2_2->GetYaxis()->CenterTitle();

	h2_2->Draw();
	h2_2_same->SetMarkerColor(kRed);
	h2_2_same->Draw("same");

	gStyle->SetOptFit(1111);
	h2_2_same->Fit("my_fit_3", "", "", t2_2_left, t2_2_right);

	//f_bt_left_2_2->SetMarkerColor(kGreen);
	f_bt_right_2_2->SetMarkerColor(kGreen);

	//f_bt_left_2_2->SetParameters(par0_tb2_2_left, par1_tb2_2_left, par2_tb2_2_left);
	f_bt_right_2_2->SetParameters(par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);

	//f_bt_left_2_2->Draw("same");
	f_bt_right_2_2->Draw("same");


	MyC_5->cd(5);

	ntpl->Draw("b3_1:t3_1>>h3_1", "", "", N_number);
	ntpl_select->Draw("b3_1:t3_1>>h3_1_same");

	h3_1->GetXaxis()->SetTitle("Time, channel number");
	h3_1->GetYaxis()->SetTitle("Amplitude, channel numbe");
	h3_1->GetXaxis()->CenterTitle();
	h3_1->GetYaxis()->CenterTitle();

	h3_1->Draw();
	h3_1_same->SetMarkerColor(kRed);
	h3_1_same->Draw("same");

	gStyle->SetOptFit(1111);
	h3_1_same->Fit("my_fit_3", "", "", t3_1_left, t3_1_right);
	
	//f_bt_left_3_1->SetMarkerColor(kGreen);
	f_bt_right_3_1->SetMarkerColor(kGreen);

	//f_bt_left_3_1->SetParameters(par0_tb3_1_left, par1_tb3_1_left, par2_tb3_1_left);
	f_bt_right_3_1->SetParameters(par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right);

	//f_bt_left_3_1->Draw("same");
	f_bt_right_3_1->Draw("same");


	MyC_5->cd(6);


	ntpl->Draw("b3_2:t3_2>>h3_2", "", "", N_number);
	ntpl_select->Draw("b3_2:t3_2>>h3_2_same");

	h3_2->GetXaxis()->SetTitle("Time, channel number");
	h3_2->GetYaxis()->SetTitle("Amplitude, channel numbe");
	h3_2->GetXaxis()->CenterTitle();
	h3_2->GetYaxis()->CenterTitle();

	h3_2->Draw();
	h3_2_same->SetMarkerColor(kRed);
	h3_2_same->Draw("same");

	gStyle->SetOptFit(1111);
	h3_2_same->Fit("my_fit_3", "", "", t3_2_left, t3_2_right);
	
	//f_bt_left_3_2->SetMarkerColor(kGreen);
	f_bt_right_3_2->SetMarkerColor(kGreen);

	//f_bt_left_3_2->SetParameters(par0_tb3_2_left, par1_tb3_2_left, par2_tb3_2_left);
	f_bt_right_3_2->SetParameters(par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right);

	//f_bt_left_3_2->Draw("same");
	f_bt_right_3_2->Draw("same");


	MyC_5->Update();


}

// время-время
void MyC_6() {



	TH2* h4_2 = new TH2D("h4_2", "t1_1:t1_2", 2100, 0, 2100, 2100, 0, 2100);
	TH2* h5_2 = new TH2D("h5_2", "t2_1:t2_2", 2100, 0, 2100, 2100, 0, 2100);
	TH2* h6_2 = new TH2D("h6_2", "t3_1:t3_2", 2100, 0, 2100, 2100, 0, 2100);

	TH2* h4_2_same = new TH2D("h4_2_same", "t1_1:t1_2+Cut", 2100, 0, 2100, 2100, 0, 2100);
	TH2* h5_2_same = new TH2D("h5_2_same", "t2_1:t2_2+Cut", 2100, 0, 2100, 2100, 0, 2100);
	TH2* h6_2_same = new TH2D("h6_2_same", "t3_1:t3_2+Cut", 2100, 0, 2100, 2100, 0, 2100);


	TCanvas* MyC_6 = new TCanvas("MyC_6", "MyC_6", 0, 0, 2000, 1000);	  


	MyC_6->Divide(3, 1);

	MyC_6->cd(1);


	ntpl->Draw("t1_1:t1_2>>h4_2", "", "", N_number);
	ntpl_select->Draw("t1_1:t1_2>>h4_2_same");

	h4_2->GetXaxis()->SetTitle("Time, channel number");
	h4_2->GetYaxis()->SetTitle("Time, channel numbe");
	h4_2->GetXaxis()->CenterTitle();
	h4_2->GetYaxis()->CenterTitle();

	h4_2->Draw();
	h4_2_same->SetMarkerColor(kRed);
	h4_2_same->Draw("same");


	MyC_6->cd(2);


	ntpl->Draw("t2_1:t2_2>>h5_2", "", "", N_number);
	ntpl_select->Draw("t2_1:t2_2>>h5_2_same");

	h5_2->GetXaxis()->SetTitle("Time, channel number");
	h5_2->GetYaxis()->SetTitle("Time, channel numbe");
	h5_2->GetXaxis()->CenterTitle();
	h5_2->GetYaxis()->CenterTitle();

	h5_2->Draw();
	h5_2_same->SetMarkerColor(kRed);
	h5_2_same->Draw("same");


	MyC_6->cd(3);

	ntpl->Draw("t3_1:t3_2>>h6_2", "", "", N_number);
	ntpl_select->Draw("t3_1:t3_2>>h6_2_same");

	h6_2->GetXaxis()->SetTitle("Time, channel number");
	h6_2->GetYaxis()->SetTitle("Time, channel numbe");
	h6_2->GetXaxis()->CenterTitle();
	h6_2->GetYaxis()->CenterTitle();

	h6_2->Draw();
	h6_2_same->SetMarkerColor(kRed);
	h6_2_same->Draw("same");


	MyC_6->Update();


}

// полусумма/разница времен без коррекции
void MyC_7() {



	TH1F* h1_time_no_corr = new TH1F("h1_time_no_corr", "((t3_1-t3_2)-0.5*((t1_1-t1_2)-(t2_1-t2_2)))/sqrt(3) WITHOUT AMPLITUDE CORRECTION", 200, -300, 300);


	TH2* h_time_amplitude_no_corr = new TH2D("h_time_amplitude_no_corr", "0.5*(t3_1-t3_2):0.5*ln(A3_1/A3_2)+Cut WITHOUT AMPLITUDE CORRECTION", 100, -5.0, 5.0, 1000, -500.0, 500.0);


	TH1F* h1_t1t2_1_no_corr = new TH1F("h1_t1t2_1_no_corr", "0.5*(t1_1+t1_2) WITHOUT AMPLITUDE CORRECTION", 1000, 500.0, 1500.0);
	TH1F* h2_t1t2_1_no_corr = new TH1F("h2_t1t2_1_no_corr", "0.5*(t2_1+t2_2) WITHOUT AMPLITUDE CORRECTION", 1000, 500.0, 1500.0);
	TH1F* h3_t1t2_1_no_corr = new TH1F("h3_t1t2_1_no_corr", "0.5*(t3_1+t3_2) WITHOUT AMPLITUDE CORRECTION", 1000, 500.0, 1500.0);

	TH1F* h1_t1t2_2_no_corr = new TH1F("h1_t1t2_2_no_corr", "0.5*(t1_1-t1_2) WITHOUT AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	TH1F* h2_t1t2_2_no_corr = new TH1F("h2_t1t2_2_no_corr", "0.5*(t2_1-t2_2) WITHOUT AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	TH1F* h3_t1t2_2_no_corr = new TH1F("h3_t1t2_2_no_corr", "0.5*(t3_1-t3_2) WITHOUT AMPLITUDE CORRECTION", 800, -400.0, 400.0);


	TCanvas* MyC_7 = new TCanvas("MyC_7", "MyC_7", 0, 0, 2000, 1000);	  



	MyC_7->Divide(3, 3);


	MyC_7->cd(1);

	TString vartime_no_corr = Form("((t3_1-t3_2)-0.5*((t1_1-t1_2)-(t2_1-t2_2)))/(%f)>>h1_time_no_corr", sqrt(3));

	h1_time_no_corr->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(vartime_no_corr);

	gStyle->SetOptFit(1111);
	h1_time_no_corr->Fit("gaus", "", "", -200, 200);



	h1_time_no_corr->GetXaxis()->SetTitle("Time, channel number");
	h1_time_no_corr->GetYaxis()->SetTitle("Events number");
	h1_time_no_corr->GetXaxis()->CenterTitle();
	h1_time_no_corr->GetYaxis()->CenterTitle();

	h1_time_no_corr->Draw();


	MyC_7->cd(2);



	TString vartb_no_corr = Form("0.5*(t3_1-t3_2):0.5*log((b3_1-%f)/(b3_2-%f))>>h_time_amplitude_no_corr", b3_1_pedestal, b3_2_pedestal);

	h_time_amplitude_no_corr->Smooth(Smooth_index);	 // Сглаживание

	
	ntpl_select->Draw(vartb_no_corr);


	h_time_amplitude_no_corr->GetXaxis()->SetTitle("Ln(A3_1/A3_2)/2");
	h_time_amplitude_no_corr->GetYaxis()->SetTitle("Time (T3_1-T3_2)/2, channel number");
	h_time_amplitude_no_corr->GetXaxis()->CenterTitle();
	h_time_amplitude_no_corr->GetYaxis()->CenterTitle();

	TF1* f1_no_corr = new TF1("f1_no_corr", "[0]+[1]*x", -5, 5);
	f1_no_corr->SetParameters(0., 1.);
	f1_no_corr->SetLineColor(kRed);
	h_time_amplitude_no_corr->Fit(f1_no_corr);
	h_time_amplitude_no_corr->Draw();
	f1_no_corr->Draw("same");


	MyC_7->cd(3);




	MyC_7->cd(4);



	TString vart1_1_no_corr = Form("0.5*(t1_1+t1_2)>>h1_t1t2_1_no_corr");

	h1_t1t2_1_no_corr->Smooth(Smooth_index);	 // Сглаживание

	ntpl_select->Draw(vart1_1_no_corr);

	h1_t1t2_1_no_corr->Fit("gaus", "", "", 700, 1300);

	h1_t1t2_1_no_corr->GetXaxis()->SetTitle("Time, channel number");
	h1_t1t2_1_no_corr->GetYaxis()->SetTitle("Events number");
	h1_t1t2_1_no_corr->GetXaxis()->CenterTitle();
	h1_t1t2_1_no_corr->GetYaxis()->CenterTitle();

	h1_t1t2_1_no_corr->Draw();

	MyC_7->cd(5);



	TString vart2_1_no_corr = Form("0.5*(t2_1+t2_2)>>h2_t1t2_1_no_corr");

	h2_t1t2_1_no_corr->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(vart2_1_no_corr);

	h2_t1t2_1_no_corr->Fit("gaus", "", "", 700, 1300);

	h2_t1t2_1_no_corr->GetXaxis()->SetTitle("Time, channel number");
	h2_t1t2_1_no_corr->GetYaxis()->SetTitle("Events number");
	h2_t1t2_1_no_corr->GetXaxis()->CenterTitle();
	h2_t1t2_1_no_corr->GetYaxis()->CenterTitle();

	h2_t1t2_1_no_corr->Draw();

	MyC_7->cd(6);



	TString vart3_1_no_corr = Form("0.5*(t3_1+t3_2)>>h3_t1t2_1_no_corr");

	h3_t1t2_1_no_corr->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(vart3_1_no_corr);

	h3_t1t2_1_no_corr->Fit("gaus", "", "", 700, 1300);

	h3_t1t2_1_no_corr->GetXaxis()->SetTitle("Time, channel number");
	h3_t1t2_1_no_corr->GetYaxis()->SetTitle("Events number");
	h3_t1t2_1_no_corr->GetXaxis()->CenterTitle();
	h3_t1t2_1_no_corr->GetYaxis()->CenterTitle();

	h3_t1t2_1_no_corr->Draw();

	MyC_7->cd(7);


	TString vart1_2_no_corr = Form("0.5*(t1_1-t1_2)>>h1_t1t2_2_no_corr");

	h1_t1t2_2_no_corr->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(vart1_2_no_corr);


	h1_t1t2_2_no_corr->GetXaxis()->SetTitle("Time, channel number");
	h1_t1t2_2_no_corr->GetYaxis()->SetTitle("Events number");
	h1_t1t2_2_no_corr->GetXaxis()->CenterTitle();
	h1_t1t2_2_no_corr->GetYaxis()->CenterTitle();

	gStyle->SetOptFit(1111);
	h1_t1t2_2_no_corr->Fit("my_fit_2", "", "", -150, 150);
	h1_t1t2_2_no_corr->Fit("gaus", "+", "", -150, 150);
	h1_t1t2_2_no_corr->Draw();

	S1 = func_2->GetParameter(1);
	X1 = func_2->GetParameter(2);
	S2 = func_2->GetParameter(3);
	X2 = func_2->GetParameter(4);


	printf("v_scint= %.3f ns/m\n", k_TDC * delta_L(S1, X1, S2, X2) / (0.01 * L_scint));

	MyC_7->cd(8);



	TString vart2_2_no_corr = Form("0.5*(t2_1-t2_2)>>h2_t1t2_2_no_corr");

	h2_t1t2_2_no_corr->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(vart2_2_no_corr);



	h2_t1t2_2_no_corr->GetXaxis()->SetTitle("Time, channel number");
	h2_t1t2_2_no_corr->GetYaxis()->SetTitle("Events number");
	h2_t1t2_2_no_corr->GetXaxis()->CenterTitle();
	h2_t1t2_2_no_corr->GetYaxis()->CenterTitle();

	gStyle->SetOptFit(1111);
	h2_t1t2_2_no_corr->Fit("my_fit_2", "", "", -150, 150);
	h2_t1t2_2_no_corr->Fit("gaus", "+", "", -150, 150);
	h2_t1t2_2_no_corr->Draw();

	S1 = func_2->GetParameter(1);
	X1 = func_2->GetParameter(2);
	S2 = func_2->GetParameter(3);
	X2 = func_2->GetParameter(4);


	printf("v_scint= %.3f ns/m\n", k_TDC * delta_L(S1, X1, S2, X2) / (0.01 * L_scint));

	MyC_7->cd(9);



	TString vart3_2_no_corr = Form("0.5*(t3_1-t3_2)>>h3_t1t2_2_no_corr");

	h3_t1t2_2_no_corr->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(vart3_2_no_corr);

	h3_t1t2_2_no_corr->GetXaxis()->SetTitle("Time, channel number");
	h3_t1t2_2_no_corr->GetYaxis()->SetTitle("Events number");
	h3_t1t2_2_no_corr->GetXaxis()->CenterTitle();
	h3_t1t2_2_no_corr->GetYaxis()->CenterTitle();

	gStyle->SetOptFit(1111);
	h3_t1t2_2_no_corr->Fit("my_fit_2", "", "", -150, 150);
	h3_t1t2_2_no_corr->Fit("gaus", "+", "", -150, 150);
	h3_t1t2_2_no_corr->Draw();

	S1 = func_2->GetParameter(1);
	X1 = func_2->GetParameter(2);
	S2 = func_2->GetParameter(3);
	X2 = func_2->GetParameter(4);

	printf("T0.5= %.3f\n", delta_L05(S1, X1, S2, X2));
	printf("T= %.3f\n", delta_L(S1, X1, S2, X2));
	printf("v_scint= %.3f ns/m\n", k_TDC * delta_L(S1, X1, S2, X2) / (0.01 * L_scint));


	MyC_7->Update();


}

// для определения порога дискриминатора (в каналах амплитуды) // по среднему геометрическому амплитуд
void MyC_8() {

	TH1F* h1_b1b2_2_1 = new TH1F("h1_b1b2_2_1", "sqrt(A1_1*A1_2) -pedestal without cuts", 1450, 50.0, 1500.0);
	TH1F* h2_b1b2_2_1 = new TH1F("h2_b1b2_2_1", "sqrt(A2_1*A2_2) -pedestal without cuts", 1450, 50.0, 1500.0);
	TH1F* h3_b1b2_2_1 = new TH1F("h3_b1b2_2_1", "sqrt(A3_1*A3_2) -pedestal without cuts", 1450, 50.0, 1500.0);

	TH1F* h1_b1b2_2_2 = new TH1F("h1_b1b2_2_2", "sqrt(A1_1*A1_2) -pedestal with cuts", 1450, 50.0, 1500.0);
	TH1F* h2_b1b2_2_2 = new TH1F("h2_b1b2_2_2", "sqrt(A2_1*A2_2) -pedestal with cuts", 1450, 50.0, 1500.0);
	TH1F* h3_b1b2_2_2 = new TH1F("h3_b1b2_2_2", "sqrt(A3_1*A3_2) -pedestal with cuts", 1450, 50.0, 1500.0);

	TH1F* h1_b1b2_2_divide = new TH1F("h1_b1b2_2_divide", "Divide (h with cuts)/h, h=sqrt(A1_1*A1_2) -pedestal", 1450, 50.0, 1500.0);
	TH1F* h2_b1b2_2_divide = new TH1F("h2_b1b2_2_divide", "Divide (h with cuts)/h, h=sqrt(A2_1*A2_2) -pedestal", 1450, 50.0, 1500.0);
	TH1F* h3_b1b2_2_divide = new TH1F("h3_b1b2_2_divide", "Divide (h with cuts)/h, h=sqrt(A3_1*A3_2) -pedestal", 1450, 50.0, 1500.0);

	TCanvas* MyC_8 = new TCanvas("MyC_8", "MyC_8", 0, 0, 2000, 1000);	  


	MyC_8->Divide(3, 3);


	MyC_8->cd(1);

	TString var4_1 = Form("sqrt((b1_1 - %f)*(b1_2 - %f)) >> h1_b1b2_2_1", b1_1_pedestal, b1_2_pedestal);

	h1_b1b2_2_1->Smooth(Smooth_index);	 // Сглаживание

	
	ntpl->Draw(var4_1, "", "", N_number);

	//	gStyle->SetOptFit(1111);
	//	h1_b1b2_2_1->Fit("landau", "", "", 0, 600);
	//	h1_b1b2_2_1->Draw();

	MyC_8->cd(2);

	TString var5_1 = Form("sqrt((b2_1 - %f)*(b2_2 - %f)) >> h2_b1b2_2_1", b2_1_pedestal, b2_2_pedestal);

	h2_b1b2_2_1->Smooth(Smooth_index);	 // Сглаживание

	ntpl->Draw(var5_1, "", "", N_number);

	//	gStyle->SetOptFit(1111);
	//	h2_b1b2_2_1->Fit("landau", "", "", 0, 600);
	//	h2_b1b2_2_1->Draw();


	MyC_8->cd(3);

	TString var6_1 = Form("sqrt((b3_1 - %f)*(b3_2 - %f)) >> h3_b1b2_2_1", b3_1_pedestal, b3_2_pedestal);

	h3_b1b2_2_1->Smooth(Smooth_index);	 // Сглаживание


	ntpl->Draw(var6_1, "", "", N_number);

	//	gStyle->SetOptFit(1111);
	//	h3_b1b2_2_1->Fit("landau", "", "", 0, 600);
	//	h3_b1b2_2_1->Draw();


	MyC_8->cd(4);

	TString var4_2 = Form("sqrt((b1_1 - %f)*(b1_2 - %f)) >> h1_b1b2_2_2", b1_1_pedestal, b1_2_pedestal);

	h1_b1b2_2_2->Smooth(Smooth_index);	 // Сглаживание


	ntpl->Draw(var4_2, cut_bb_3, "", N_number);

	//	h1_b1b2_2_divide = (TH1F*)h1_b1b2_2_2->Clone();



	//	gStyle->SetOptFit(1111);
	//	h1_b1b2_2_2->Fit("landau", "", "", 0, 600);
	//	h1_b1b2_2_2->Draw();


	MyC_8->cd(5);

	TString var5_2 = Form("sqrt((b2_1 - %f)*(b2_2 - %f)) >> h2_b1b2_2_2", b2_1_pedestal, b2_2_pedestal);

	h2_b1b2_2_2->Smooth(Smooth_index);	 // Сглаживание


	ntpl->Draw(var5_2, cut_bb_3, "", N_number);

	//	h2_b1b2_2_divide = (TH1F*)h2_b1b2_2_2->Clone();

	//	gStyle->SetOptFit(1111);
	//	h2_b1b2_2_2->Fit("landau", "", "", 0, 600);
	//	h2_b1b2_2_2->Draw();


	MyC_8->cd(6);

	TString var6_2 = Form("sqrt((b3_1 - %f)*(b3_2 - %f)) >> h3_b1b2_2_2", b3_1_pedestal, b3_2_pedestal);

	h3_b1b2_2_2->Smooth(Smooth_index);	 // Сглаживание


	ntpl->Draw(var6_2, cut_bb_1 && cut_bb_2, "", N_number);

	//	h3_b1b2_2_divide = (TH1F*)h3_b1b2_2_2->Clone();

	//	gStyle->SetOptFit(1111);
	//	h3_b1b2_2_2->Fit("landau", "", "", 0, 600);
	//	h3_b1b2_2_2->Draw();



	MyC_8->cd(7);

	h1_b1b2_2_divide = (TH1F*)h1_b1b2_2_2->Clone();
	h1_b1b2_2_divide->GetXaxis()->SetTitle(" ");
	h1_b1b2_2_divide->GetYaxis()->SetTitle(" ");
	h1_b1b2_2_divide->SetTitle("h1/h2");
	h1_b1b2_2_divide->Divide(h1_b1b2_2_1);
	//	h1_b1b2_2_divide->Fit("","0");

	h1_b1b2_2_divide->Smooth(Smooth_index);	 // Сглаживание

	h1_b1b2_2_divide->Draw("COL");


	MyC_8->cd(8);

	h2_b1b2_2_divide = (TH1F*)h2_b1b2_2_2->Clone();
	h2_b1b2_2_divide->GetXaxis()->SetTitle(" ");
	h2_b1b2_2_divide->GetYaxis()->SetTitle(" ");
	h2_b1b2_2_divide->SetTitle("h1/h2");
	h2_b1b2_2_divide->Divide(h2_b1b2_2_1);

	h2_b1b2_2_divide->Smooth(Smooth_index);	 // Сглаживание


	h2_b1b2_2_divide->Draw("COL");

	MyC_8->cd(9);

	h3_b1b2_2_divide = (TH1F*)h3_b1b2_2_2->Clone();
	h3_b1b2_2_divide->GetXaxis()->SetTitle(" ");
	h3_b1b2_2_divide->GetYaxis()->SetTitle(" ");
	h3_b1b2_2_divide->SetTitle("h1/h2");
	h3_b1b2_2_divide->Divide(h3_b1b2_2_1);

	h3_b1b2_2_divide->Smooth(Smooth_index);	 // Сглаживание


	h3_b1b2_2_divide->Draw("COL");


	MyC_8->Update();



}

// для определения порого дискриминатора (в каналах амплитуды) // по каждому каналу амплитуды
void MyC_9() {



	TH1F* hb1_1_no_cut = new TH1F("hb1_1_no_cut", "b1_1 -pedestal without cuts", 1480, 20.0, 1500.0);
	TH1F* hb1_2_no_cut = new TH1F("hb1_2_no_cut", "b1_2 -pedestal without cuts", 1480, 20.0, 1500.0);
	TH1F* hb2_1_no_cut = new TH1F("hb2_1_no_cut", "b2_1 -pedestal without cuts", 1480, 20.0, 1500.0);
	TH1F* hb2_2_no_cut = new TH1F("hb2_2_no_cut", "b2_2 -pedestal without cuts", 1480, 20.0, 1500.0);
	TH1F* hb3_1_no_cut = new TH1F("hb3_1_no_cut", "b3_1 -pedestal without cuts", 1480, 20.0, 1500.0);
	TH1F* hb3_2_no_cut = new TH1F("hb3_2_no_cut", "b3_2 -pedestal without cuts", 1480, 20.0, 1500.0);


	TH1F* hb1_1_cut = new TH1F("hb1_1_cut", "b1_1 with cuts", 1480, 20.0, 1500.0);
	TH1F* hb1_2_cut = new TH1F("hb1_2_cut", "b1_2 with cuts", 1480, 20.0, 1500.0);
	TH1F* hb2_1_cut = new TH1F("hb2_1_cut", "b2_1 with cuts", 1480, 20.0, 1500.0);
	TH1F* hb2_2_cut = new TH1F("hb2_2_cut", "b2_2 with cuts", 1480, 20.0, 1500.0);
	TH1F* hb3_1_cut = new TH1F("hb3_1_cut", "b3_1 with cuts", 1480, 20.0, 1500.0);
	TH1F* hb3_2_cut = new TH1F("hb3_2_cut", "b3_2 with cuts", 1480, 20.0, 1500.0);


	TH1F* hb1_1_divide = new TH1F("hb1_1_divide", "Divide (h with cuts)/h, h=b1_1", 195, 5.0, 200.0);
	TH1F* hb1_2_divide = new TH1F("hb1_2_divide", "Divide (h with cuts)/h, h=b1_2", 195, 5.0, 200.0);
	TH1F* hb2_1_divide = new TH1F("hb2_1_divide", "Divide (h with cuts)/h, h=b2_1", 195, 5.0, 200.0);
	TH1F* hb2_2_divide = new TH1F("hb2_2_divide", "Divide (h with cuts)/h, h=b2_2", 195, 5.0, 200.0);
	TH1F* hb3_1_divide = new TH1F("hb3_1_divide", "Divide (h with cuts)/h, h=b3_1", 195, 5.0, 200.0);
	TH1F* hb3_2_divide = new TH1F("hb3_2_divide", "Divide (h with cuts)/h, h=b3_2", 195, 5.0, 200.0);


	const TCut cut_t1_1 = "t1_1>10&&t1_1<2000";
	const TCut cut_t1_2 = "t1_2>10&&t1_2<2000";
	const TCut cut_t2_1 = "t2_1>10&&t2_1<2000";
	const TCut cut_t2_2 = "t2_2>10&&t2_2<2000";
	const TCut cut_t3_1 = "t3_1>10&&t3_1<2000";
	const TCut cut_t3_2 = "t3_2>10&&t3_2<2000";

	TCanvas* MyC_9 = new TCanvas("MyC_9", "MyC_9", 0, 0, 2000, 1000);	  


	MyC_9->Divide(6, 3);


	MyC_9->cd(1);

	TString varb1_1_1 = Form("b1_1 - %f >> hb1_1_no_cut", b1_1_pedestal);

	hb1_1_no_cut->Smooth(Smooth_index);	 // Сглаживание

	
	ntpl->Draw(varb1_1_1, "", "", N_number);

	MyC_9->cd(2);

	TString varb1_2_1 = Form("b1_2 - %f >> hb1_2_no_cut", b1_2_pedestal);

	hb1_2_no_cut->Smooth(Smooth_index);	 // Сглаживание


	ntpl->Draw(varb1_2_1, "", "", N_number);


	MyC_9->cd(3);

	TString varb2_1_1 = Form("b2_1 - %f >> hb2_1_no_cut", b2_1_pedestal);

	hb2_1_no_cut->Smooth(Smooth_index);	 // Сглаживание


	ntpl->Draw(varb2_1_1, "", "", N_number);


	MyC_9->cd(4);

	TString varb2_2_1 = Form("b2_2 - %f >> hb2_2_no_cut", b2_2_pedestal);

	hb2_2_no_cut->Smooth(Smooth_index);	 // Сглаживание


	ntpl->Draw(varb2_2_1, "", "", N_number);


	MyC_9->cd(5);

	TString varb3_1_1 = Form("b3_1 - %f >> hb3_1_no_cut", b3_1_pedestal);

	hb3_1_no_cut->Smooth(Smooth_index);	 // Сглаживание


	ntpl->Draw(varb3_1_1, "", "", N_number);


	MyC_9->cd(6);

	TString varb3_2_1 = Form("b3_2 - %f >> hb3_2_no_cut", b3_2_pedestal);

	hb3_2_no_cut->Smooth(Smooth_index);	 // Сглаживание


	ntpl->Draw(varb3_2_1, "", "", N_number);


	MyC_9->cd(7);

	TString varb1_1_2 = Form("b1_1 - %f >> hb1_1_cut", b1_1_pedestal);

	hb1_1_cut->Smooth(Smooth_index);	 // Сглаживание


	ntpl->Draw(varb1_1_2, cut_t1_1, "", N_number);


	MyC_9->cd(8);

	TString varb1_2_2 = Form("b1_2 - %f >> hb1_2_cut", b1_2_pedestal);

	hb1_2_cut->Smooth(Smooth_index);	 // Сглаживание


	ntpl->Draw(varb1_2_2, cut_t1_2, "", N_number);


	MyC_9->cd(9);

	TString varb2_1_2 = Form("b2_1 - %f >> hb2_1_cut", b2_1_pedestal);

	hb2_1_cut->Smooth(Smooth_index);	 // Сглаживание


	ntpl->Draw(varb2_1_2, cut_t2_1, "", N_number);


	MyC_9->cd(10);

	TString varb2_2_2 = Form("b2_2 - %f >> hb2_2_cut", b2_2_pedestal);

	hb2_2_cut->Smooth(Smooth_index);	 // Сглаживание


	ntpl->Draw(varb2_2_2, cut_t2_2, "", N_number);


	MyC_9->cd(11);

	TString varb3_1_2 = Form("b3_1 - %f >> hb3_1_cut", b3_1_pedestal);

	hb3_1_cut->Smooth(Smooth_index);	 // Сглаживание


	ntpl->Draw(varb3_1_2, cut_t3_1, "", N_number);
	MyC_9->cd(12);

	TString varb3_2_2 = Form("b3_2 - %f >> hb3_2_cut", b3_2_pedestal);

	hb3_2_cut->Smooth(Smooth_index);	 // Сглаживание


	ntpl->Draw(varb3_2_2, cut_t3_2, "", N_number);

	MyC_9->cd(13);

	hb1_1_divide = (TH1F*)hb1_1_cut->Clone();
	hb1_1_divide->GetXaxis()->SetTitle(" ");
	hb1_1_divide->GetYaxis()->SetTitle(" ");
	hb1_1_divide->SetTitle("h1/h2");
	hb1_1_divide->Divide(hb1_1_no_cut);
	//	hb1_1_divide->Fit("","0");

	hb1_1_divide->Smooth(Smooth_index);	 // Сглаживание


	hb1_1_divide->Draw("COL");

	int i = 1;
	for (i = 1; i < 100; i++) {

		if (hb1_1_divide->GetBinContent(i) > 0.5) break;
	}

	double y2 = hb1_1_divide->GetBinContent(i);
	double y1 = hb1_1_divide->GetBinContent(i-1);
	double x2 = hb1_1_divide->GetBinCenter(i);
	double x1 = hb1_1_divide->GetBinCenter(i - 1);
	double k = (x2 - x1) / (y2 - y1);
	double p = x1 + (0.5 - y1) * k;

//	TLatex* tl = new TLatex();
//	tl->DrawLatex(10, 0.55, Form("threshold=%4.1f", p));

	cout << "1_1: threshold-pedestal=" << p << endl;
	cout << "1_1: threshold=" << p + b1_1_pedestal << endl;
	cout << "1_1: pedestal=" << b1_1_pedestal << endl;
	cout << endl;



	MyC_9->cd(14);

	hb1_2_divide = (TH1F*)hb1_2_cut->Clone();
	hb1_2_divide->GetXaxis()->SetTitle(" ");
	hb1_2_divide->GetYaxis()->SetTitle(" ");
	hb1_2_divide->SetTitle("h1/h2");
	hb1_2_divide->Divide(hb1_2_no_cut);
	//	hb1_2_divide->Fit("","0");

	hb1_2_divide->Smooth(Smooth_index);	 // Сглаживание


	hb1_2_divide->Draw("COL");

	for (i = 1; i < 100; i++) {

		if (hb1_2_divide->GetBinContent(i) > 0.5) break;
	}

	y2 = hb1_2_divide->GetBinContent(i);
	y1 = hb1_2_divide->GetBinContent(i - 1);
	x2 = hb1_2_divide->GetBinCenter(i);
	x1 = hb1_2_divide->GetBinCenter(i - 1);
	k = (x2 - x1) / (y2 - y1);
	p = x1 + (0.5 - y1) * k;

	//	TLatex* tl = new TLatex();
	//	tl->DrawLatex(10, 0.55, Form("threshold=%4.1f", p));

	cout << "1_2: threshold-pedestal=" << p << endl;
	cout << "1_2: threshold=" << p + b1_2_pedestal << endl;
	cout << "1_2: pedestal=" << b1_2_pedestal << endl;
	cout << endl;

	MyC_9->cd(15);

	hb2_1_divide = (TH1F*)hb2_1_cut->Clone();
	hb2_1_divide->GetXaxis()->SetTitle(" ");
	hb2_1_divide->GetYaxis()->SetTitle(" ");
	hb2_1_divide->SetTitle("h1/h2");
	hb2_1_divide->Divide(hb2_1_no_cut);
	//	hb2_1_divide->Fit("","0");

	hb2_1_divide->Smooth(Smooth_index);	 // Сглаживание


	hb2_1_divide->Draw("COL");

	for (i = 1; i < 100; i++) {

		if (hb2_1_divide->GetBinContent(i) > 0.5) break;
	}

	y2 = hb2_1_divide->GetBinContent(i);
	y1 = hb2_1_divide->GetBinContent(i - 1);
	x2 = hb2_1_divide->GetBinCenter(i);
	x1 = hb2_1_divide->GetBinCenter(i - 1);
	k = (x2 - x1) / (y2 - y1);
	p = x1 + (0.5 - y1) * k;

	//	TLatex* tl = new TLatex();
	//	tl->DrawLatex(10, 0.55, Form("threshold=%4.1f", p));

	cout << "2_1: threshold-pedestal=" << p << endl;
	cout << "2_1: threshold=" << p + b2_1_pedestal << endl;
	cout << "2_1: pedestal=" << b2_1_pedestal << endl;
	cout << endl;



	MyC_9->cd(16);
	hb2_2_divide = (TH1F*)hb2_2_cut->Clone();
	hb2_2_divide->GetXaxis()->SetTitle(" ");
	hb2_2_divide->GetYaxis()->SetTitle(" ");
	hb2_2_divide->SetTitle("h1/h2");
	hb2_2_divide->Divide(hb2_2_no_cut);
	//	hb2_2_divide->Fit("","0");

	hb2_2_divide->Smooth(Smooth_index);	 // Сглаживание


	hb2_2_divide->Draw("COL");


	for (i = 1; i < 100; i++) {

		if (hb2_2_divide->GetBinContent(i) > 0.5) break;
	}

	y2 = hb2_2_divide->GetBinContent(i);
	y1 = hb2_2_divide->GetBinContent(i - 1);
	x2 = hb2_2_divide->GetBinCenter(i);
	x1 = hb2_2_divide->GetBinCenter(i - 1);
	k = (x2 - x1) / (y2 - y1);
	p = x1 + (0.5 - y1) * k;

	//	TLatex* tl = new TLatex();
	//	tl->DrawLatex(10, 0.55, Form("threshold=%4.1f", p));

	cout << "2_2: threshold-pedestal=" << p << endl;
	cout << "2_2: threshold=" << p + b2_2_pedestal << endl;
	cout << "2_2: pedestal=" << b2_2_pedestal << endl;
	cout << endl;


	MyC_9->cd(17);

	hb3_1_divide = (TH1F*)hb3_1_cut->Clone();
	hb3_1_divide->GetXaxis()->SetTitle(" ");
	hb3_1_divide->GetYaxis()->SetTitle(" ");
	hb3_1_divide->SetTitle("h1/h2");
	hb3_1_divide->Divide(hb3_1_no_cut);
	//	hb3_1_divide->Fit("","0");

	hb3_1_divide->Smooth(Smooth_index);	 // Сглаживание


	hb3_1_divide->Draw("COL");
  
	for (i = 1; i < 100; i++) {

		if (hb3_1_divide->GetBinContent(i) > 0.5) break;
	}

	y2 = hb3_1_divide->GetBinContent(i);
	y1 = hb3_1_divide->GetBinContent(i - 1);
	x2 = hb3_1_divide->GetBinCenter(i);
	x1 = hb3_1_divide->GetBinCenter(i - 1);
	k = (x2 - x1) / (y2 - y1);
	p = x1 + (0.5 - y1) * k;

	//	TLatex* tl = new TLatex();
	//	tl->DrawLatex(10, 0.55, Form("threshold=%4.1f", p));

	cout << "3_1: threshold-pedestal=" << p << endl;
	cout << "3_1: threshold=" << p + b3_1_pedestal << endl;
	cout << "3_1: pedestal=" << b3_1_pedestal << endl;
	cout << endl;

	MyC_9->cd(18);

	hb3_2_divide = (TH1F*)hb3_2_cut->Clone();
	hb3_2_divide->GetXaxis()->SetTitle(" ");
	hb3_2_divide->GetYaxis()->SetTitle(" ");
	hb3_2_divide->SetTitle("h1/h2");
	hb3_2_divide->Divide(hb3_2_no_cut);
	//	hb3_2_divide->Fit("","0");

	hb3_2_divide->Smooth(Smooth_index);	 // Сглаживание


	hb3_2_divide->Draw("COL");

	for (i = 1; i < 100; i++) {

		if (hb3_2_divide->GetBinContent(i) > 0.5) break;
	}

	y2 = hb3_2_divide->GetBinContent(i);
	y1 = hb3_2_divide->GetBinContent(i - 1);
	x2 = hb3_2_divide->GetBinCenter(i);
	x1 = hb3_2_divide->GetBinCenter(i - 1);
	k = (x2 - x1) / (y2 - y1);
	p = x1 + (0.5 - y1) * k;

	//	TLatex* tl = new TLatex();
	//	tl->DrawLatex(10, 0.55, Form("threshold=%4.1f", p));

	cout << "3_2: threshold-pedestal=" << p << endl;
	cout << "3_2: threshold=" << p + b3_2_pedestal << endl;
	cout << "3_2: pedestal=" << b3_2_pedestal << endl;
	cout << endl;


	MyC_9->Update();

}

// для определения числа ф.э.
void MyC_10() {

	TCanvas* MyC_10 = new TCanvas("MyC_10", "MyC_10", 0, 0, 2000, 1000);	  



	


	TH1F* h1_b1b2_phe_1 = new TH1F("h1_b1b2_phe_1", "A1_1/A1_2 -pedestal", 200, -2.0, 20.0);
	TH1F* h2_b1b2_phe_1 = new TH1F("h2_b1b2_phe_1", "A2_1/A2_2 -pedestal", 200, -2.0, 20.0);
	TH1F* h3_b1b2_phe_1 = new TH1F("h3_b1b2_phe_1", "A3_1/A3_2 -pedestal", 200, -2.0, 20.0);

	TH1F* h1_b1b2_phe_2 = new TH1F("h1_b1b2_phe_2", "(A1_1-A1_2)/(A1_1+A1_2) -pedestal", 200, -2.0, 2.0);
	TH1F* h2_b1b2_phe_2 = new TH1F("h2_b1b2_phe_2", "(A2_1-A2_2)/(A2_1+A2_2) -pedestal", 200, -2.0, 2.0);
	TH1F* h3_b1b2_phe_2 = new TH1F("h3_b1b2_phe_2", "(A3_1-A3_2)/(A3_1+A3_2) -pedestal", 200, -2.0, 2.0);


//	TF1* g1 = (TF1*)h1_b1b2_phe_1->GetListOfFunctions()->FindObject("landau");
//	TF1* g2 = (TF1*)h2_b1b2_phe_1->GetListOfFunctions()->FindObject("landau");
//	TF1* g3 = (TF1*)h3_b1b2_phe_1->GetListOfFunctions()->FindObject("landau");
	

		// set the parameters to the mean and RMS of the histogram
	func_1->SetParameters(10.0, 2.0, -1.0, 2.0, -1.0);
	
	// give the parameters meaningful names
	func_1->SetParNames("A", "S1", "X1", "S2", "X2");
	


	Double_t sigma = 0.0;
	
	MyC_10->Divide(3, 2);

	MyC_10->cd(1);


	TString var1 = Form("(b1_1 - %f) / (b1_2 - %f) >> h1_b1b2_phe_1", b1_1_pedestal, b1_2_pedestal);

	h1_b1b2_phe_1->Smooth(Smooth_index);	 // Сглаживание
		
	ntpl_select->Draw(var1);

	h1_b1b2_phe_1->GetXaxis()->SetTitle("A1/A2");
	h1_b1b2_phe_1->GetYaxis()->SetTitle("Events number");
	h1_b1b2_phe_1->GetXaxis()->CenterTitle();
	h1_b1b2_phe_1->GetYaxis()->CenterTitle();

	gStyle->SetOptFit(1111);
	h1_b1b2_phe_1->Fit("landau", "", "", 0.0, 10.0);
	
	TF1* g1 = (TF1*)h1_b1b2_phe_1->GetListOfFunctions()->FindObject("landau");
		
	sigma = g1->GetParameter(2);

	printf("N_ph.e.=%f \n", 2.0 / pow(sigma, 2.0));

	h1_b1b2_phe_1->Draw();
	

	MyC_10->cd(2);

	TString var2 = Form("(b2_1 - %f) / (b2_2 - %f) >> h2_b1b2_phe_1", b2_1_pedestal, b2_2_pedestal);

	h2_b1b2_phe_1->Smooth(Smooth_index);	 // Сглаживание

	ntpl_select->Draw(var2);

	h2_b1b2_phe_1->GetXaxis()->SetTitle("A1/A2");
	h2_b1b2_phe_1->GetYaxis()->SetTitle("Events number");
	h2_b1b2_phe_1->GetXaxis()->CenterTitle();
	h2_b1b2_phe_1->GetYaxis()->CenterTitle();

	gStyle->SetOptFit(1111);
	h2_b1b2_phe_1->Fit("landau", "", "", 0.0, 10.0);
		
	TF1* g2 = (TF1*)h2_b1b2_phe_1->GetListOfFunctions()->FindObject("landau");
	
	sigma = g2->GetParameter(2);

	printf("N_ph.e.=%f \n", 2.0 / pow(sigma, 2.0));

	h2_b1b2_phe_1->Draw();


	MyC_10->cd(3);

	TString var3 = Form("(b3_1 - %f) / (b3_2 - %f) >> h3_b1b2_phe_1", b3_1_pedestal, b3_2_pedestal);

	h3_b1b2_phe_1->Smooth(Smooth_index);	 // Сглаживание

	ntpl_select->Draw(var3);

	h3_b1b2_phe_1->GetXaxis()->SetTitle("A1/A2");
	h3_b1b2_phe_1->GetYaxis()->SetTitle("Events number");
	h3_b1b2_phe_1->GetXaxis()->CenterTitle();
	h3_b1b2_phe_1->GetYaxis()->CenterTitle();

	gStyle->SetOptFit(1111);
	h3_b1b2_phe_1->Fit("landau", "", "", 0.0, 10.0);
		
	TF1* g3 = (TF1*)h3_b1b2_phe_1->GetListOfFunctions()->FindObject("landau");
	
	sigma = g3->GetParameter(2);

	printf("N_ph.e.=%f \n", 2.0 / pow(sigma, 2.0));
	
	h3_b1b2_phe_1->Draw();

	MyC_10->cd(4);

	TString var4 = Form("((b1_1 - %f) - (b1_2 - %f)) / ((b1_1 - %f) + (b1_2 - %f)) >> h1_b1b2_phe_2", b1_1_pedestal, b1_2_pedestal, b1_1_pedestal, b1_2_pedestal);

	h1_b1b2_phe_2->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(var4);

	h1_b1b2_phe_2->GetXaxis()->SetTitle("(A1-A2)/(A1+A2)");
	h1_b1b2_phe_2->GetYaxis()->SetTitle("Events number");
	h1_b1b2_phe_2->GetXaxis()->CenterTitle();
	h1_b1b2_phe_2->GetYaxis()->CenterTitle();

	gStyle->SetOptFit(1111);
	h1_b1b2_phe_2->Fit("my_fit_1", "", "", -2, 2);

	S1 = func_1->GetParameter(1);
	X1 = func_1->GetParameter(2);
	S2 = func_1->GetParameter(3);
	X2 = func_1->GetParameter(4);

	printf("N_ph.e.=%f \n", 2.0 / pow(delta_L(S1, X1, S2, X2), 2.0));


	h1_b1b2_phe_2->Draw();

	MyC_10->cd(5);

	TString var5 = Form("((b2_1 - %f) - (b2_2 - %f)) / ((b2_1 - %f) + (b2_2 - %f)) >> h2_b1b2_phe_2", b2_1_pedestal, b2_2_pedestal, b2_1_pedestal, b2_2_pedestal);

	h2_b1b2_phe_2->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(var5);

	h2_b1b2_phe_2->GetXaxis()->SetTitle("(A1-A2)/(A1+A2)");
	h2_b1b2_phe_2->GetYaxis()->SetTitle("Events number");
	h2_b1b2_phe_2->GetXaxis()->CenterTitle();
	h2_b1b2_phe_2->GetYaxis()->CenterTitle();

	gStyle->SetOptFit(1111);
	h2_b1b2_phe_2->Fit("my_fit_1", "", "", -2, 2);

	S1 = func_1->GetParameter(1);
	X1 = func_1->GetParameter(2);
	S2 = func_1->GetParameter(3);
	X2 = func_1->GetParameter(4);

	printf("N_ph.e.=%f \n", 2.0 / pow(delta_L(S1, X1, S2, X2), 2.0));


	h2_b1b2_phe_2->Draw();

	MyC_10->cd(6);

	TString var6 = Form("((b3_1 - %f) - (b3_2 - %f)) / ((b3_1 - %f) + (b3_2 - %f)) >> h3_b1b2_phe_2", b3_1_pedestal, b3_2_pedestal, b3_1_pedestal, b3_2_pedestal);

	h3_b1b2_phe_2->Smooth(Smooth_index);	 // Сглаживание

	
	ntpl_select->Draw(var6);

	h3_b1b2_phe_2->GetXaxis()->SetTitle("(A1-A2)/(A1+A2)");
	h3_b1b2_phe_2->GetYaxis()->SetTitle("Events number");
	h3_b1b2_phe_2->GetXaxis()->CenterTitle();
	h3_b1b2_phe_2->GetYaxis()->CenterTitle();

	gStyle->SetOptFit(1111);
	h3_b1b2_phe_2->Fit("my_fit_1", "", "", -2, 2);

	S1 = func_1->GetParameter(1);
	X1 = func_1->GetParameter(2);
	S2 = func_1->GetParameter(3);
	X2 = func_1->GetParameter(4);

	printf("N_ph.e.=%f \n", 2.0 / pow(delta_L(S1, X1, S2, X2), 2.0));
	
	h3_b1b2_phe_2->Draw();

	MyC_10->Update();


}

// временное разрешение, решение системы уравнений - поправка на амплитуду
void MyC_11() {


	TH1F* h_time_tr_1 = new TH1F("h_time_tr_1", "((t1_1+t1_2)*0.5-(t2_1+t2_2)*0.5 WITH AMPLITUDE CORRECTION", 200, -300, 300);
	TH1F* h_time_tr_2 = new TH1F("h_time_tr_2", "(t1_1-t1_2)*0.5-(t2_1-t2_2)*0.5 WITH AMPLITUDE CORRECTION", 200, -300, 300);

	TH1F* h_time_31_1 = new TH1F("h_time_31_1", "(t3_1+t3_2)*0.5-(t1_1+t1_2)*0.5 WITH AMPLITUDE CORRECTION", 200, -300, 300);
	TH1F* h_time_31_2 = new TH1F("h_time_31_2", "(t3_1-t3_2)*0.5-(t1_1-t1_2)*0.5 WITH AMPLITUDE CORRECTION", 200, -300, 300);

	TH1F* h_time_32_1 = new TH1F("h_time_32_1", "(t3_1+t3_2)*0.5-(t2_1+t2_2)*0.5 WITH AMPLITUDE CORRECTION", 200, -300, 300);
	TH1F* h_time_32_2 = new TH1F("h_time_32_2", "(t3_1-t3_2)*0.5-(t2_1-t2_2)*0.5 WITH AMPLITUDE CORRECTION", 200, -300, 300);


	Double_t sigma_tr_1 = 0.0, sigma_31_1 = 0.0, sigma_32_1 = 0.0;
	Double_t sigma_tr_2 = 0.0, sigma_31_2 = 0.0, sigma_32_2 = 0.0;



	TCanvas* MyC_11 = new TCanvas("MyC_11", "MyC_11", 0, 0, 2000, 1000);	  // для определения числа ф.э.

	MyC_11->Divide(3, 2);


	MyC_11->cd(1);

	TString vartime_tr_1 = Form("((t1_1-delta_t(b1_1, %f, %f, %f))+(t1_2-delta_t(b1_2, %f, %f, %f)))*0.5-((t2_1-delta_t(b2_1, %f, %f, %f))+(t2_2-delta_t(b2_2, %f, %f, %f)))*0.5>>h_time_tr_1", par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right, par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);


	ntpl_select->Draw(vartime_tr_1);

	h_time_tr_1->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h_time_tr_1->Fit("gaus", "", "", -200, 200);


	TF1* g1_1 = (TF1*)h_time_tr_1->GetListOfFunctions()->FindObject("gaus");

	sigma_tr_1 = g1_1->GetParameter(2);


	h_time_tr_1->GetXaxis()->SetTitle("Time, channel number");
	h_time_tr_1->GetYaxis()->SetTitle("Events number");
	h_time_tr_1->GetXaxis()->CenterTitle();
	h_time_tr_1->GetYaxis()->CenterTitle();

	h_time_tr_1->Draw();


	MyC_11->cd(2);

	TString vartime_31_1 = Form("((t3_1-delta_t(b3_1, %f, %f, %f))+(t3_2-delta_t(b3_2, %f, %f, %f)))*0.5-((t1_1-delta_t(b1_1, %f, %f, %f))+(t1_2-delta_t(b1_2, %f, %f, %f)))*0.5>>h_time_31_1", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right);


	ntpl_select->Draw(vartime_31_1);

	h_time_31_1->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h_time_31_1->Fit("gaus", "", "", -200, 200);


	TF1* g2_1 = (TF1*)h_time_31_1->GetListOfFunctions()->FindObject("gaus");

	sigma_31_1 = g2_1->GetParameter(2);


	h_time_31_1->GetXaxis()->SetTitle("Time, channel number");
	h_time_31_1->GetYaxis()->SetTitle("Events number");
	h_time_31_1->GetXaxis()->CenterTitle();
	h_time_31_1->GetYaxis()->CenterTitle();

	h_time_31_1->Draw();


	MyC_11->cd(3);

	TString vartime_32_1 = Form("((t3_1-delta_t(b3_1, %f, %f, %f))+(t3_2-delta_t(b3_2, %f, %f, %f)))*0.5-((t2_1-delta_t(b2_1, %f, %f, %f))+(t2_2-delta_t(b2_2, %f, %f, %f)))*0.5>>h_time_32_1", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);


	ntpl_select->Draw(vartime_32_1);

	h_time_32_1->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h_time_32_1->Fit("gaus", "", "", -200, 200);


	TF1* g3_1 = (TF1*)h_time_32_1->GetListOfFunctions()->FindObject("gaus");

	sigma_32_1 = g3_1->GetParameter(2);


	h_time_32_1->GetXaxis()->SetTitle("Time, channel number");
	h_time_32_1->GetYaxis()->SetTitle("Events number");
	h_time_32_1->GetXaxis()->CenterTitle();
	h_time_32_1->GetYaxis()->CenterTitle();

	h_time_32_1->Draw();


	printf("Sigma t #3=%.3f\n", sigma_t(sigma_31_1, sigma_32_1, sigma_tr_1));


	MyC_11->cd(4);

	TString vartime_tr_2 = Form("((t1_1-delta_t(b1_1, %f, %f, %f))-(t1_2-delta_t(b1_2, %f, %f, %f)))*0.5-((t2_1-delta_t(b2_1, %f, %f, %f))-(t2_2-delta_t(b2_2, %f, %f, %f)))*0.5>>h_time_tr_2", par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right, par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);


	ntpl_select->Draw(vartime_tr_2);

	h_time_tr_2->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h_time_tr_2->Fit("gaus", "", "", -200, 200);


	TF1* g1_2 = (TF1*)h_time_tr_2->GetListOfFunctions()->FindObject("gaus");

	sigma_tr_2 = g1_2->GetParameter(2);


	h_time_tr_2->GetXaxis()->SetTitle("Time, channel number");
	h_time_tr_2->GetYaxis()->SetTitle("Events number");
	h_time_tr_2->GetXaxis()->CenterTitle();
	h_time_tr_2->GetYaxis()->CenterTitle();

	h_time_tr_2->Draw();


	MyC_11->cd(5);

	TString vartime_31_2 = Form("((t3_1-delta_t(b3_1, %f, %f, %f))-(t3_2-delta_t(b3_2, %f, %f, %f)))*0.5-((t1_1-delta_t(b1_1, %f, %f, %f))-(t1_2-delta_t(b1_2, %f, %f, %f)))*0.5>>h_time_31_2", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right);


	ntpl_select->Draw(vartime_31_2);

	h_time_31_2->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h_time_31_2->Fit("gaus", "", "", -200, 200);


	TF1* g2_2 = (TF1*)h_time_31_2->GetListOfFunctions()->FindObject("gaus");

	sigma_31_2 = g2_2->GetParameter(2);


	h_time_31_2->GetXaxis()->SetTitle("Time, channel number");
	h_time_31_2->GetYaxis()->SetTitle("Events number");
	h_time_31_2->GetXaxis()->CenterTitle();
	h_time_31_2->GetYaxis()->CenterTitle();

	h_time_31_2->Draw();


	MyC_11->cd(6);

	TString vartime_32_2 = Form("((t3_1-delta_t(b3_1, %f, %f, %f))-(t3_2-delta_t(b3_2, %f, %f, %f)))*0.5-((t2_1-delta_t(b2_1, %f, %f, %f))-(t2_2-delta_t(b2_2, %f, %f, %f)))*0.5>>h_time_32_2", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);


	ntpl_select->Draw(vartime_32_2);

	h_time_32_2->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h_time_32_2->Fit("gaus", "", "", -200, 200);


	TF1* g3_2 = (TF1*)h_time_32_2->GetListOfFunctions()->FindObject("gaus");

	sigma_32_2 = g3_2->GetParameter(2);


	h_time_32_2->GetXaxis()->SetTitle("Time, channel number");
	h_time_32_2->GetYaxis()->SetTitle("Events number");
	h_time_32_2->GetXaxis()->CenterTitle();
	h_time_32_2->GetYaxis()->CenterTitle();

	h_time_32_2->Draw();


	printf("Sigma t #3=%.3f\n", sigma_t(sigma_31_2, sigma_32_2, sigma_tr_2));





	MyC_11->Update();

}

// временное разрешение, решение системы уравнений без поправки на амплитуду
void MyC_12() {


	TH1F* h_time_tr_1 = new TH1F("h_time_tr_1", "(t1_1+t1_2)*0.5-(t2_1+t2_2)*0.5 WITHOUT AMPLITUDE CORRECTION", 200, -300, 300);
	TH1F* h_time_tr_2 = new TH1F("h_time_tr_2", "(t1_1-t1_2)*0.5-(t2_1-t2_2)*0.5 WITHOUT AMPLITUDE CORRECTION", 200, -300, 300);

	TH1F* h_time_31_1 = new TH1F("h_time_31_1", "(t3_1+t3_2)*0.5-(t1_1+t1_2)*0.5 WITHOUT AMPLITUDE CORRECTION", 200, -300, 300);
	TH1F* h_time_31_2 = new TH1F("h_time_31_2", "(t3_1-t3_2)*0.5-(t1_1-t1_2)*0.5 WITHOUT AMPLITUDE CORRECTION", 200, -300, 300);

	TH1F* h_time_32_1 = new TH1F("h_time_32_1", "(t3_1+t3_2)*0.5-(t2_1+t2_2)*0.5 WITHOUT AMPLITUDE CORRECTION", 200, -300, 300);
	TH1F* h_time_32_2 = new TH1F("h_time_32_2", "(t3_1-t3_2)*0.5-(t2_1-t2_2)*0.5 WITHOUT AMPLITUDE CORRECTION", 200, -300, 300);


	Double_t sigma_tr_1 = 0.0, sigma_31_1 = 0.0, sigma_32_1 = 0.0;
	Double_t sigma_tr_2 = 0.0, sigma_31_2 = 0.0, sigma_32_2 = 0.0;



	TCanvas* MyC_12 = new TCanvas("MyC_12", "MyC_12", 0, 0, 2000, 1000);	  // для определения числа ф.э.

	MyC_12->Divide(3, 2);


	MyC_12->cd(1);

	TString vartime_tr_1 = Form("(t1_1+t1_2)*0.5-(t2_1+t2_2)*0.5>>h_time_tr_1");


	ntpl_select->Draw(vartime_tr_1);

	h_time_tr_1->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h_time_tr_1->Fit("gaus", "", "", -200, 200);


	TF1* g1_1 = (TF1*)h_time_tr_1->GetListOfFunctions()->FindObject("gaus");

	sigma_tr_1 = g1_1->GetParameter(2);


	h_time_tr_1->GetXaxis()->SetTitle("Time, channel number");
	h_time_tr_1->GetYaxis()->SetTitle("Events number");
	h_time_tr_1->GetXaxis()->CenterTitle();
	h_time_tr_1->GetYaxis()->CenterTitle();

	h_time_tr_1->Draw();


	MyC_12->cd(2);

	TString vartime_31_1 = Form("(t3_1+t3_2)*0.5-(t1_1+t1_2)*0.5>>h_time_31_1");


	ntpl_select->Draw(vartime_31_1);

	h_time_31_1->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h_time_31_1->Fit("gaus", "", "", -200, 200);


	TF1* g2_1 = (TF1*)h_time_31_1->GetListOfFunctions()->FindObject("gaus");

	sigma_31_1 = g2_1->GetParameter(2);


	h_time_31_1->GetXaxis()->SetTitle("Time, channel number");
	h_time_31_1->GetYaxis()->SetTitle("Events number");
	h_time_31_1->GetXaxis()->CenterTitle();
	h_time_31_1->GetYaxis()->CenterTitle();

	h_time_31_1->Draw();


	MyC_12->cd(3);

	TString vartime_32_1 = Form("(t3_1+t3_2)*0.5-(t2_1+t2_2)*0.5>>h_time_32_1");


	ntpl_select->Draw(vartime_32_1);

	h_time_32_1->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h_time_32_1->Fit("gaus", "", "", -200, 200);


	TF1* g3_1 = (TF1*)h_time_32_1->GetListOfFunctions()->FindObject("gaus");

	sigma_32_1 = g3_1->GetParameter(2);


	h_time_32_1->GetXaxis()->SetTitle("Time, channel number");
	h_time_32_1->GetYaxis()->SetTitle("Events number");
	h_time_32_1->GetXaxis()->CenterTitle();
	h_time_32_1->GetYaxis()->CenterTitle();

	h_time_32_1->Draw();


	printf("Sigma t #3=%.3f\n", sigma_t(sigma_31_1, sigma_32_1, sigma_tr_1));


	MyC_12->cd(4);

	TString vartime_tr_2 = Form("(t1_1-t1_2)*0.5-(t2_1-t2_2)*0.5>>h_time_tr_2");


	ntpl_select->Draw(vartime_tr_2);

	h_time_tr_2->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h_time_tr_2->Fit("gaus", "", "", -200, 200);


	TF1* g1_2 = (TF1*)h_time_tr_2->GetListOfFunctions()->FindObject("gaus");

	sigma_tr_2 = g1_2->GetParameter(2);


	h_time_tr_2->GetXaxis()->SetTitle("Time, channel number");
	h_time_tr_2->GetYaxis()->SetTitle("Events number");
	h_time_tr_2->GetXaxis()->CenterTitle();
	h_time_tr_2->GetYaxis()->CenterTitle();

	h_time_tr_2->Draw();


	MyC_12->cd(5);


	TString vartime_31_2 = Form("(t3_1-t3_2)*0.5-(t1_1-t1_2)*0.5>>h_time_31_2");


	ntpl_select->Draw(vartime_31_2);

	h_time_31_2->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h_time_31_2->Fit("gaus", "", "", -200, 200);


	TF1* g2_2 = (TF1*)h_time_31_2->GetListOfFunctions()->FindObject("gaus");

	sigma_31_2 = g2_2->GetParameter(2);


	h_time_31_2->GetXaxis()->SetTitle("Time, channel number");
	h_time_31_2->GetYaxis()->SetTitle("Events number");
	h_time_31_2->GetXaxis()->CenterTitle();
	h_time_31_2->GetYaxis()->CenterTitle();

	h_time_31_2->Draw();


	MyC_12->cd(6);


	TString vartime_32_2 = Form("(t3_1-t3_2)*0.5-(t2_1-t2_2)*0.5>>h_time_32_2");


	ntpl_select->Draw(vartime_32_2);

	h_time_32_2->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h_time_32_2->Fit("gaus", "", "", -200, 200);


	TF1* g3_2 = (TF1*)h_time_32_2->GetListOfFunctions()->FindObject("gaus");

	sigma_32_2 = g3_2->GetParameter(2);


	h_time_32_2->GetXaxis()->SetTitle("Time, channel number");
	h_time_32_2->GetYaxis()->SetTitle("Events number");
	h_time_32_2->GetXaxis()->CenterTitle();
	h_time_32_2->GetYaxis()->CenterTitle();

	h_time_32_2->Draw();


	printf("Sigma t #3=%.3f\n", sigma_t(sigma_31_2, sigma_32_2, sigma_tr_2));





	MyC_12->Update();

}

// различные комбинации времени - поправка на амплитуду
void MyC_13() {

	
	///
	
	TH1F* h1 = new TH1F("h1", "0.5*((t3_1+t3_2)*0.5-0.5*((t1_1-t1_2)*0.5+(t2_1-t2_2)*0.5)) WITH AMPLITUDE CORRECTION", 200, -300, 300);

	TH1F* h2 = new TH1F("h2", "0.5*((t3_1+t3_2)*0.5-0.5*((t1_1-t1_2)*0.5-(t2_1-t2_2)*0.5)) WITH AMPLITUDE CORRECTION", 200, -300, 300);

	TH1F* h3 = new TH1F("h3", "0.5*((t3_1+t3_2)*0.5-0.5*((t1_1+t1_2)*0.5-(t2_1+t2_2)*0.5)) WITH AMPLITUDE CORRECTION", 200, -300, 300);

	TH1F* h4 = new TH1F("h4", "0.5*((t3_1+t3_2)*0.5-0.5*((t1_1+t1_2)*0.5+(t2_1+t2_2)*0.5)) WITH AMPLITUDE CORRECTION", 200, -300, 300);
	///
	
	TH1F* h5 = new TH1F("h5", "0.5*((t3_1-t3_2)*0.5-0.5*((t1_1-t1_2)*0.5+(t2_1-t2_2)*0.5)) WITH AMPLITUDE CORRECTION", 200, -300, 300);

	TH1F* h6 = new TH1F("h6", "((t3_1-t3_2)*0.5-0.5*((t1_1-t1_2)*0.5+(t2_1-t2_2)*0.5))/sqrt(3) WITH AMPLITUDE CORRECTION", 200, -300, 300);
	
	TH1F* h7 = new TH1F("h7", "0.5*((t3_1-t3_2)*0.5-0.5*((t1_1-t1_2)*0.5-(t2_1-t2_2)*0.5)) WITH AMPLITUDE CORRECTION", 200, -300, 300);

	TH1F* h8 = new TH1F("h8", "0.5*((t3_1-t3_2)*0.5-0.5*((t1_1+t1_2)*0.5-(t2_1+t2_2)*0.5)) WITH AMPLITUDE CORRECTION", 200, -300, 300);

	TH1F* h9 = new TH1F("h9", "0.5*((t3_1-t3_2)*0.5-0.5*((t1_1+t1_2)*0.5+(t2_1+t2_2)*0.5)) WITH AMPLITUDE CORRECTION", 200, -300, 300);

	TCanvas* MyC_13 = new TCanvas("MyC_13", "MyC_13", 0, 0, 2000, 1000);

	MyC_13->Divide(3, 3);


	MyC_13->cd(1);

	TString var1 = Form("0.5*((((t3_1 - delta_t(b3_1, % f, % f, % f)) + (t3_2 - delta_t(b3_2, % f, % f, % f))) * 0.5 - 0.5 * (((t1_1 - delta_t(b1_1, % f, % f, % f)) - (t1_2 - delta_t(b1_2, % f, % f, % f))) * 0.5 + ((t2_1 - delta_t(b2_1, % f, % f, % f)) - (t2_2 - delta_t(b2_2, % f, % f, % f))) * 0.5)))>>h1", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right, par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);

	ntpl_select->Draw(var1);

	h1->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h1->Fit("gaus", "", "", -100, 100);



	h1->GetXaxis()->SetTitle("Time, channel number");
	h1->GetYaxis()->SetTitle("Events number");
	h1->GetXaxis()->CenterTitle();
	h1->GetYaxis()->CenterTitle();

	h1->Draw();

	MyC_13->cd(2);

	TString var2 = Form("0.5*((((t3_1 - delta_t(b3_1, % f, % f, % f)) + (t3_2 - delta_t(b3_2, % f, % f, % f))) * 0.5 - 0.5 * (((t1_1 - delta_t(b1_1, % f, % f, % f)) - (t1_2 - delta_t(b1_2, % f, % f, % f))) * 0.5 - ((t2_1 - delta_t(b2_1, % f, % f, % f)) - (t2_2 - delta_t(b2_2, % f, % f, % f))) * 0.5)))>>h2", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right, par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);

	ntpl_select->Draw(var2);

	h2->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h2->Fit("gaus", "", "", -100, 100);



	h2->GetXaxis()->SetTitle("Time, channel number");
	h2->GetYaxis()->SetTitle("Events number");
	h2->GetXaxis()->CenterTitle();
	h2->GetYaxis()->CenterTitle();

	h2->Draw();

	MyC_13->cd(3);

	TString var3 = Form("0.5*((((t3_1 - delta_t(b3_1, % f, % f, % f)) + (t3_2 - delta_t(b3_2, % f, % f, % f))) * 0.5 - 0.5 * (((t1_1 - delta_t(b1_1, % f, % f, % f)) + (t1_2 - delta_t(b1_2, % f, % f, % f))) * 0.5 - ((t2_1 - delta_t(b2_1, % f, % f, % f)) + (t2_2 - delta_t(b2_2, % f, % f, % f))) * 0.5)))>>h3", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right, par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);

	ntpl_select->Draw(var3);

	h3->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h3->Fit("gaus", "", "", -100, 100);



	h3->GetXaxis()->SetTitle("Time, channel number");
	h3->GetYaxis()->SetTitle("Events number");
	h3->GetXaxis()->CenterTitle();
	h3->GetYaxis()->CenterTitle();

	h3->Draw();

	MyC_13->cd(4);

	TString var4 = Form("0.5*((((t3_1 - delta_t(b3_1, % f, % f, % f)) + (t3_2 - delta_t(b3_2, % f, % f, % f))) * 0.5 - 0.5 * (((t1_1 - delta_t(b1_1, % f, % f, % f)) + (t1_2 - delta_t(b1_2, % f, % f, % f))) * 0.5 + ((t2_1 - delta_t(b2_1, % f, % f, % f)) + (t2_2 - delta_t(b2_2, % f, % f, % f))) * 0.5)))>>h4", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right, par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);

	ntpl_select->Draw(var4);

	h4->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h4->Fit("gaus", "", "", -100, 100);



	h4->GetXaxis()->SetTitle("Time, channel number");
	h4->GetYaxis()->SetTitle("Events number");
	h4->GetXaxis()->CenterTitle();
	h4->GetYaxis()->CenterTitle();

	h4->Draw();


	/////



	MyC_13->cd(5);

	TString var5 = Form("0.5*((((t3_1 - delta_t(b3_1, % f, % f, % f)) - (t3_2 - delta_t(b3_2, % f, % f, % f))) * 0.5 - 0.5 * (((t1_1 - delta_t(b1_1, % f, % f, % f)) - (t1_2 - delta_t(b1_2, % f, % f, % f))) * 0.5 + ((t2_1 - delta_t(b2_1, % f, % f, % f)) - (t2_2 - delta_t(b2_2, % f, % f, % f))) * 0.5)))>>h5", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right, par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);

	ntpl_select->Draw(var5);

	h5->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h5->Fit("gaus", "", "", -100, 100);



	h5->GetXaxis()->SetTitle("Time, channel number");
	h5->GetYaxis()->SetTitle("Events number");
	h5->GetXaxis()->CenterTitle();
	h5->GetYaxis()->CenterTitle();

	h5->Draw();

	MyC_13->cd(6);

	TString var6 = Form("(((t3_1 - delta_t(b3_1, % f, % f, % f)) - (t3_2 - delta_t(b3_2, % f, % f, % f))) * 0.5 - 0.5 * (((t1_1 - delta_t(b1_1, % f, % f, % f)) - (t1_2 - delta_t(b1_2, % f, % f, % f))) * 0.5 + ((t2_1 - delta_t(b2_1, % f, % f, % f)) - (t2_2 - delta_t(b2_2, % f, % f, % f))) * 0.5)) / % f>>h6", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right, par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right, sqrt(3));

	ntpl_select->Draw(var6);

	h6->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h6->Fit("gaus", "", "", -100, 100);



	h6->GetXaxis()->SetTitle("Time, channel number");
	h6->GetYaxis()->SetTitle("Events number");
	h6->GetXaxis()->CenterTitle();
	h6->GetYaxis()->CenterTitle();

	h6->Draw();


	MyC_13->cd(7);

	TString var7 = Form("0.5*((((t3_1 - delta_t(b3_1, % f, % f, % f)) - (t3_2 - delta_t(b3_2, % f, % f, % f))) * 0.5 - 0.5 * (((t1_1 - delta_t(b1_1, % f, % f, % f)) - (t1_2 - delta_t(b1_2, % f, % f, % f))) * 0.5 - ((t2_1 - delta_t(b2_1, % f, % f, % f)) - (t2_2 - delta_t(b2_2, % f, % f, % f))) * 0.5)))>>h7", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right, par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);

	ntpl_select->Draw(var7);

	h7->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h7->Fit("gaus", "", "", -100, 100);



	h7->GetXaxis()->SetTitle("Time, channel number");
	h7->GetYaxis()->SetTitle("Events number");
	h7->GetXaxis()->CenterTitle();
	h7->GetYaxis()->CenterTitle();

	h7->Draw();

	MyC_13->cd(8);

	TString var8 = Form("0.5*((((t3_1 - delta_t(b3_1, % f, % f, % f)) - (t3_2 - delta_t(b3_2, % f, % f, % f))) * 0.5 - 0.5 * (((t1_1 - delta_t(b1_1, % f, % f, % f)) + (t1_2 - delta_t(b1_2, % f, % f, % f))) * 0.5 - ((t2_1 - delta_t(b2_1, % f, % f, % f)) + (t2_2 - delta_t(b2_2, % f, % f, % f))) * 0.5)))>>h8", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right, par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);

	ntpl_select->Draw(var8);

	h8->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h8->Fit("gaus", "", "", -100, 100);



	h8->GetXaxis()->SetTitle("Time, channel number");
	h8->GetYaxis()->SetTitle("Events number");
	h8->GetXaxis()->CenterTitle();
	h8->GetYaxis()->CenterTitle();

	h8->Draw();

	MyC_13->cd(9);

	TString var9 = Form("0.5*((((t3_1 - delta_t(b3_1, % f, % f, % f)) - (t3_2 - delta_t(b3_2, % f, % f, % f))) * 0.5 - 0.5 * (((t1_1 - delta_t(b1_1, % f, % f, % f)) + (t1_2 - delta_t(b1_2, % f, % f, % f))) * 0.5 + ((t2_1 - delta_t(b2_1, % f, % f, % f)) + (t2_2 - delta_t(b2_2, % f, % f, % f))) * 0.5)))>>h9", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right, par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);

	ntpl_select->Draw(var9);

	h9->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h9->Fit("gaus", "", "", -100, 100);



	h9->GetXaxis()->SetTitle("Time, channel number");
	h9->GetYaxis()->SetTitle("Events number");
	h9->GetXaxis()->CenterTitle();
	h9->GetYaxis()->CenterTitle();

	h9->Draw();




}

// различные комбинации времени без поправки на амплитуду
void MyC_14() {


	///

	TH1F* h1 = new TH1F("h1", "0.5*((t3_1+t3_2)*0.5-0.5*((t1_1-t1_2)*0.5+(t2_1-t2_2)*0.5)) WITHOUT AMPLITUDE CORRECTION", 200, -3000, 3000);

	TH1F* h2 = new TH1F("h2", "0.5*((t3_1+t3_2)*0.5-0.5*((t1_1-t1_2)*0.5-(t2_1-t2_2)*0.5)) WITHOUT AMPLITUDE CORRECTION", 200, -3000, 3000);

	TH1F* h3 = new TH1F("h3", "0.5*((t3_1+t3_2)*0.5-0.5*((t1_1+t1_2)*0.5-(t2_1+t2_2)*0.5)) WITHOUT AMPLITUDE CORRECTION", 200, -3000, 3000);

	TH1F* h4 = new TH1F("h4", "0.5*((t3_1+t3_2)*0.5-0.5*((t1_1+t1_2)*0.5+(t2_1+t2_2)*0.5)) WITHOUT AMPLITUDE CORRECTION", 200, -3000, 3000);
	///

	TH1F* h5 = new TH1F("h5", "0.5*((t3_1-t3_2)*0.5-0.5*((t1_1-t1_2)*0.5+(t2_1-t2_2)*0.5)) WITHOUT AMPLITUDE CORRECTION", 200, -3000, 3000);

	TH1F* h6 = new TH1F("h6", "((t3_1-t3_2)*0.5-0.5*((t1_1-t1_2)*0.5+(t2_1-t2_2)*0.5))/sqrt(3) WITHOUT AMPLITUDE CORRECTION", 200, -3000, 3000);

	TH1F* h7 = new TH1F("h7", "0.5*((t3_1-t3_2)*0.5-0.5*((t1_1-t1_2)*0.5-(t2_1-t2_2)*0.5)) WITHOUT AMPLITUDE CORRECTION", 200, -3000, 3000);

	TH1F* h8 = new TH1F("h8", "0.5*((t3_1-t3_2)*0.5-0.5*((t1_1+t1_2)*0.5-(t2_1+t2_2)*0.5)) WITHOUT AMPLITUDE CORRECTION", 200, -3000, 3000);

	TH1F* h9 = new TH1F("h9", "0.5*((t3_1-t3_2)*0.5-0.5*((t1_1+t1_2)*0.5+(t2_1+t2_2)*0.5)) WITHOUT AMPLITUDE CORRECTION", 200, -3000, 3000);

	TCanvas* MyC_14 = new TCanvas("MyC_14", "MyC_14", 0, 0, 2000, 1000);

	MyC_14->Divide(3, 3);


	MyC_14->cd(1);

	TString var1 = Form("0.5*(((t3_1 + t3_2) * 0.5 - 0.5 * ((t1_1 - t1_2) * 0.5 + (t2_1 - t2_2) * 0.5)))>>h1");

	ntpl_select->Draw(var1);

	h1->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h1->Fit("gaus", "", "", -100, 100);



	h1->GetXaxis()->SetTitle("Time, channel number");
	h1->GetYaxis()->SetTitle("Events number");
	h1->GetXaxis()->CenterTitle();
	h1->GetYaxis()->CenterTitle();

	h1->Draw();

	MyC_14->cd(2);


	TString var2 = Form("0.5*(((t3_1 + t3_2) * 0.5 - 0.5 * ((t1_1 - t1_2) * 0.5 - (t2_1 - t2_2) * 0.5)))>>h2");


	ntpl_select->Draw(var2);

	h2->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h2->Fit("gaus", "", "", -100, 100);



	h2->GetXaxis()->SetTitle("Time, channel number");
	h2->GetYaxis()->SetTitle("Events number");
	h2->GetXaxis()->CenterTitle();
	h2->GetYaxis()->CenterTitle();

	h2->Draw();

	MyC_14->cd(3);

	TString var3 = Form("0.5*(((t3_1 + t3_2) * 0.5 - 0.5 * ((t1_1 + t1_2) * 0.5 - (t2_1 + t2_2) * 0.5)))>>h3");

	ntpl_select->Draw(var3);

	h3->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h3->Fit("gaus", "", "", -100, 100);



	h3->GetXaxis()->SetTitle("Time, channel number");
	h3->GetYaxis()->SetTitle("Events number");
	h3->GetXaxis()->CenterTitle();
	h3->GetYaxis()->CenterTitle();

	h3->Draw();

	MyC_14->cd(4);

	TString var4 = Form("0.5*(((t3_1 + t3_2) * 0.5 - 0.5 * ((t1_1 + t1_2) * 0.5 + (t2_1 + t2_2) * 0.5)))>>h4");

	ntpl_select->Draw(var4);

	h4->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h4->Fit("gaus", "", "", -100, 100);



	h4->GetXaxis()->SetTitle("Time, channel number");
	h4->GetYaxis()->SetTitle("Events number");
	h4->GetXaxis()->CenterTitle();
	h4->GetYaxis()->CenterTitle();

	h4->Draw();


	/////



	MyC_14->cd(5);

	TString var5 = Form("0.5*(((t3_1 - t3_2) * 0.5 - 0.5 * ((t1_1 - t1_2) * 0.5 + (t2_1 - t2_2) * 0.5)))>>h5");


	ntpl_select->Draw(var5);

	h5->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h5->Fit("gaus", "", "", -100, 100);



	h5->GetXaxis()->SetTitle("Time, channel number");
	h5->GetYaxis()->SetTitle("Events number");
	h5->GetXaxis()->CenterTitle();
	h5->GetYaxis()->CenterTitle();

	h5->Draw();

	MyC_14->cd(6);


	TString var6 = Form("(((t3_1 - t3_2) * 0.5 - 0.5 * ((t1_1 - t1_2) * 0.5 + (t2_1 - t2_2) * 0.5))) / % f>>h6", sqrt(3));

	
	ntpl_select->Draw(var6);

	h6->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h6->Fit("gaus", "", "", -100, 100);



	h6->GetXaxis()->SetTitle("Time, channel number");
	h6->GetYaxis()->SetTitle("Events number");
	h6->GetXaxis()->CenterTitle();
	h6->GetYaxis()->CenterTitle();

	h6->Draw();


	MyC_14->cd(7);

	TString var7 = Form("0.5*(((t3_1 - t3_2) * 0.5 - 0.5 * ((t1_1 - t1_2) * 0.5 - (t2_1 - t2_2) * 0.5)))>>h7");


	ntpl_select->Draw(var7);

	h7->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h7->Fit("gaus", "", "", -100, 100);



	h7->GetXaxis()->SetTitle("Time, channel number");
	h7->GetYaxis()->SetTitle("Events number");
	h7->GetXaxis()->CenterTitle();
	h7->GetYaxis()->CenterTitle();

	h7->Draw();

	MyC_14->cd(8);

	TString var8 = Form("0.5*(((t3_1 - t3_2) * 0.5 - 0.5 * ((t1_1 + t1_2) * 0.5 - (t2_1 + t2_2) * 0.5)))>>h8");

	ntpl_select->Draw(var8);

	h8->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h8->Fit("gaus", "", "", -100, 100);



	h8->GetXaxis()->SetTitle("Time, channel number");
	h8->GetYaxis()->SetTitle("Events number");
	h8->GetXaxis()->CenterTitle();
	h8->GetYaxis()->CenterTitle();

	h8->Draw();

	MyC_14->cd(9);

	TString var9 = Form("0.5*(((t3_1 - t3_2) * 0.5 - 0.5 * ((t1_1 + t1_2) * 0.5 + (t2_1 + t2_2) * 0.5)))>>h9");

	ntpl_select->Draw(var9);

	h9->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h9->Fit("gaus", "", "", -100, 100);



	h9->GetXaxis()->SetTitle("Time, channel number");
	h9->GetYaxis()->SetTitle("Events number");
	h9->GetXaxis()->CenterTitle();
	h9->GetYaxis()->CenterTitle();

	h9->Draw();




}


// временное разрешение	+ координата
void MyC_1_(Int_t channel) {


	TH1F* h1_time = new TH1F("h1_time", "((t3_1-t3_2)-0.5*((t1_1-t1_2)-(t2_1-t2_2)))/sqrt(3) WITH AMPLITUDE CORRECTION", 200, -300, 300);

	TH2* h_time_amplitude = new TH2D("h_time_amplitude", "0.5*(t3_1-t3_2):0.5*ln(A3_1/A3_2)+Cut WITH AMPLITUDE CORRECTION", 100, -5.0, 5.0, 1000, -500.0, 500.0);

	TH1F* h1_t1t2_1 = new TH1F("h1_t1t2_1", "0.5*(t1_1+t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	TH1F* h2_t1t2_1 = new TH1F("h2_t1t2_1", "0.5*(t2_1+t2_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	TH1F* h3_t1t2_1 = new TH1F("h3_t1t2_1", "0.5*(t3_1+t3_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);

	TH1F* h1_t1t2_2 = new TH1F("h1_t1t2_2", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	TH1F* h2_t1t2_2 = new TH1F("h2_t1t2_2", "0.5*(t2_1-t2_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	TH1F* h3_t1t2_2 = new TH1F("h3_t1t2_2", "0.5*(t3_1-t3_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);



	TCanvas* MyC_1 = new TCanvas("MyC_1", "MyC_1", 0, 0, 2000, 1000);	  // для определения числа ф.э.

	MyC_1->Divide(3, 3);


	MyC_1->cd(1);

	TString vartime = Form("((t3_1-delta_t(b3_1, %f, %f, %f))-(t3_2-delta_t(b3_2, %f, %f, %f))-0.5*(((t1_1-delta_t(b1_1, %f, %f, %f))-(t1_2-delta_t(b1_2, %f, %f, %f)))-((t2_1-delta_t(b2_1, %f, %f, %f))-(t2_2-delta_t(b2_2, %f, %f, %f)))))/(%f)>>h1_time", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right, par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right, sqrt(3));


	ntpl_select->Draw(vartime);

	h1_time->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h1_time->Fit("gaus", "", "", -200, 200);



	h1_time->GetXaxis()->SetTitle("Time, channel number");
	h1_time->GetYaxis()->SetTitle("Events number");
	h1_time->GetXaxis()->CenterTitle();
	h1_time->GetYaxis()->CenterTitle();

	h1_time->Draw();


	MyC_1->cd(2);



	TString vartb = Form("0.5*((t3_1-delta_t(b3_1, %f, %f, %f))-(t3_2-delta_t(b3_2, %f, %f, %f))):0.5*log((b3_1-%f)/(b3_2-%f))>>h_time_amplitude", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, b3_1_pedestal, b3_2_pedestal);

	ntpl_select->Draw(vartb);


	h_time_amplitude->GetXaxis()->SetTitle("Ln(A3_1/A3_2)/2");
	h_time_amplitude->GetYaxis()->SetTitle("Time (T3_1-T3_2)/2, channel number");
	h_time_amplitude->GetXaxis()->CenterTitle();
	h_time_amplitude->GetYaxis()->CenterTitle();

	TF1* f1 = new TF1("f1", "[0]+[1]*x", -5, 5);
	f1->SetParameters(0., 1.);
	f1->SetLineColor(kRed);
	h_time_amplitude->Fit(f1);
	h_time_amplitude->Draw();
	f1->Draw("same");


	MyC_1->cd(3);




	MyC_1->cd(4);



	TString vart1_1 = Form("0.5*((t1_1-delta_t(b1_1, %f, %f, %f))+(t1_2-delta_t(b1_2, %f, %f, %f)))>>h1_t1t2_1", par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right);


	ntpl_select->Draw(vart1_1);

	h1_t1t2_1->Smooth(Smooth_index);	 // Сглаживание

	h1_t1t2_1->Fit("gaus", "", "", -200, 200);

	h1_t1t2_1->GetXaxis()->SetTitle("Time, channel number");
	h1_t1t2_1->GetYaxis()->SetTitle("Events number");
	h1_t1t2_1->GetXaxis()->CenterTitle();
	h1_t1t2_1->GetYaxis()->CenterTitle();


	h1_t1t2_1->Draw();

	MyC_1->cd(5);



	TString vart2_1 = Form("0.5*((t2_1-delta_t(b2_1, %f, %f, %f))+(t2_2-delta_t(b2_2, %f, %f, %f)))>>h2_t1t2_1", par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);


	ntpl_select->Draw(vart2_1);

	h2_t1t2_1->Smooth(Smooth_index);	 // Сглаживание

	h2_t1t2_1->Fit("gaus", "", "", -200, 200);

	h2_t1t2_1->GetXaxis()->SetTitle("Time, channel number");
	h2_t1t2_1->GetYaxis()->SetTitle("Events number");
	h2_t1t2_1->GetXaxis()->CenterTitle();
	h2_t1t2_1->GetYaxis()->CenterTitle();


	h2_t1t2_1->Draw();

	MyC_1->cd(6);



	TString vart3_1 = Form("0.5*((t3_1-delta_t(b3_1, %f, %f, %f))+(t3_2-delta_t(b3_2, %f, %f, %f)))>>h3_t1t2_1", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right);


	ntpl_select->Draw(vart3_1);

	h3_t1t2_1->Smooth(Smooth_index);	 // Сглаживание

	h3_t1t2_1->Fit("gaus", "", "", -200, 200);

	h3_t1t2_1->GetXaxis()->SetTitle("Time, channel number");
	h3_t1t2_1->GetYaxis()->SetTitle("Events number");
	h3_t1t2_1->GetXaxis()->CenterTitle();
	h3_t1t2_1->GetYaxis()->CenterTitle();


	h3_t1t2_1->Draw();

	MyC_1->cd(7);


	TString vart1_2 = Form("0.5*((t1_1-delta_t(b1_1, %f, %f, %f))-(t1_2-delta_t(b1_2, %f, %f, %f)))>>h1_t1t2_2", par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right);

	h1_t1t2_2->Smooth(Smooth_index);	 // Сглаживание

	ntpl_select->Draw(vart1_2);



	h1_t1t2_2->GetXaxis()->SetTitle("Time, channel number");
	h1_t1t2_2->GetYaxis()->SetTitle("Events number");
	h1_t1t2_2->GetXaxis()->CenterTitle();
	h1_t1t2_2->GetYaxis()->CenterTitle();

	gStyle->SetOptFit(1111);
	h1_t1t2_2->Fit("gaus", "", "", -150, 150);
	h1_t1t2_2->Fit("my_fit_2", "+", "", -150, 150);


	h1_t1t2_2->Draw();

	S1 = func_2->GetParameter(1);
	X1 = func_2->GetParameter(2);
	S2 = func_2->GetParameter(3);
	X2 = func_2->GetParameter(4);


	printf("v_scint= %.3f ns/m\n", k_TDC * delta_L(S1, X1, S2, X2) / (0.01 * L_scint));

	MyC_1->cd(8);



	TString vart2_2 = Form("0.5*((t2_1-delta_t(b2_1, %f, %f, %f))-(t2_2-delta_t(b2_2, %f, %f, %f)))>>h2_t1t2_2", par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);

	h2_t1t2_2->Smooth(Smooth_index);	 // Сглаживание

	ntpl_select->Draw(vart2_2);



	h2_t1t2_2->GetXaxis()->SetTitle("Time, channel number");
	h2_t1t2_2->GetYaxis()->SetTitle("Events number");
	h2_t1t2_2->GetXaxis()->CenterTitle();
	h2_t1t2_2->GetYaxis()->CenterTitle();

	gStyle->SetOptFit(1111);
	h2_t1t2_2->Fit("gaus", "", "", -150, 150);
	h2_t1t2_2->Fit("my_fit_2", "+", "", -150, 150);
	h2_t1t2_2->Draw();

	S1 = func_2->GetParameter(1);
	X1 = func_2->GetParameter(2);
	S2 = func_2->GetParameter(3);
	X2 = func_2->GetParameter(4);


	printf("v_scint= %.3f ns/m\n", k_TDC * delta_L(S1, X1, S2, X2) / (0.01 * L_scint));

	MyC_1->cd(9);



	TString vart3_2 = Form("0.5*((t3_1-delta_t(b3_1, %f, %f, %f))-(t3_2-delta_t(b3_2,  %f, %f, %f)))>>h3_t1t2_2", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right);

	h3_t1t2_2->Smooth(Smooth_index);	 // Сглаживание

	ntpl_select->Draw(vart3_2);

	h3_t1t2_2->GetXaxis()->SetTitle("Time, channel number");
	h3_t1t2_2->GetYaxis()->SetTitle("Events number");
	h3_t1t2_2->GetXaxis()->CenterTitle();
	h3_t1t2_2->GetYaxis()->CenterTitle();

	gStyle->SetOptFit(1111);
	h3_t1t2_2->Fit("gaus", "", "", -150, 150);
	h3_t1t2_2->Fit("my_fit_2", "+", "", -150, 150);
	h3_t1t2_2->Draw();

	S1 = func_2->GetParameter(1);
	X1 = func_2->GetParameter(2);
	S2 = func_2->GetParameter(3);
	X2 = func_2->GetParameter(4);

	printf("v_scint= %.3f ns/m\n", k_TDC * delta_L(S1, X1, S2, X2) / (0.01 * L_scint));

	MyC_1->Update();

}





// геометрическая эффективность
void eff() {

	TTree* ntpl_off = ntpl->CloneTree(N_number);


	///
	///
	///
	/// Расчет эффективностей
	///
	///
	///
	///
	///


	printf("eff (1&2&3)/(1&2)= %.3f %%\n", (((Double_t)ntpl_off->GetEntries(cut_bb_1 && cut_bb_2 && cut_bb_3)) / ((Double_t)ntpl_off->GetEntries(cut_bb_1 && cut_bb_2))) * 100.0);

	printf("eff (1&2)/(1)= %.3f %%\n", (((Double_t)ntpl_off->GetEntries(cut_bb_1 && cut_bb_2)) / ((Double_t)ntpl_off->GetEntries(cut_bb_1))) * 100.0);

	printf("eff (1&2)/(2)= %.3f %%\n", (((Double_t)ntpl_off->GetEntries(cut_bb_1 && cut_bb_2)) / ((Double_t)ntpl_off->GetEntries(cut_bb_2))) * 100.0);

	printf("eff (1&3)/(1)= %.3f %%\n", (((Double_t)ntpl_off->GetEntries(cut_bb_1 && cut_bb_3)) / ((Double_t)ntpl_off->GetEntries(cut_bb_1))) * 100.0);

	printf("eff (2&3)/(2)= %.3f %%\n", (((Double_t)ntpl_off->GetEntries(cut_bb_2 && cut_bb_3)) / ((Double_t)ntpl_off->GetEntries(cut_bb_2))) * 100.0);

	printf("eff (1&2&3)/(3)= %.3f %%\n", (((Double_t)ntpl_off->GetEntries(cut_bb_1 && cut_bb_2 && cut_bb_3)) / ((Double_t)ntpl_off->GetEntries(cut_bb_3))) * 100.0);

	printf("1-eff (1&2)/(3)= %.3f %%\n", (1.0 - ((Double_t)ntpl_off->GetEntries(cut_bb_1 && cut_bb_2)) / ((Double_t)ntpl_off->GetEntries(cut_bb_3))) * 100.0);

}


int main() {

	
	
	
	//gROOT->SetBatch();

	gStyle->SetOptStat(111111); // draw statistics on plots,
							// (0) for no output
	gStyle->SetOptFit(1111);    // draw fit results on plot,


	// The number of workers
	const UInt_t nWorkers = 4U;

	Int_t nthreads = 8;
	ROOT::EnableImplicitMT(nthreads);


	// Create the task group and give work to it
	ROOT::Experimental::TTaskGroup tg;

	ROOT::EnableThreadSafety();

	//ROOT::EnableImplicitMT();

	
	//for (Long64_t i = 0; i < N_number; ++i) {
	//	ntpl->GetEntry(i); // параллельное чтение
	//}

	//f->cd();

   
	// set the parameters to the mean and RMS of the histogram
	func_1->SetParameters(10.0, 5.0, -5.0, 5.0, -5.0);
	func_2->SetParameters(10.0, 1.0, -100.0, 1.0, -100.0);
	// give the parameters meaningful names
	func_1->SetParNames("A", "S1", "X1", "S2", "X2");
	func_2->SetParNames("A", "S1", "X1", "S2", "X2");

	
	func_3->SetParameters(10000.0, 1000.0, 10.0);
	func_3->SetParNames("k", "t0", "b0");


	



	///
	/// 
	/// 
	/// Отрисовка канвы и гистограмм
	/// 
	/// 
	/// 
	///



	//f.Close();

	 	//ntpl_select->Delete("all");

	






	tg.Run(workItem0);
	tg.Run([&]() {
				
		printf("Running workItem1...\n"); 
		
		MyC_1();
		MyC_2();
		MyC_3();
		MyC_4();
		MyC_5();
		MyC_6();
		MyC_7();
		MyC_8();
		MyC_9();
		MyC_10();
		MyC_11();
		MyC_12();
		MyC_13();
		MyC_14();

		});


	printf("Running something in the \"main\" thread\n");

	// Wait until all items are complete
	tg.Wait();

	printf("All work completed.\n");

   
	





	return 0;
}
