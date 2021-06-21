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


#define Z_CHANNEL 0   // канал Z=1...16 стрипа для ограничения координаты вдоль сцинтиллятора
					   // если Z_CHANNEL ==0, то координатную зависимость не учитывать, беруться все координаты вдоль счетчика

#define X_CHANNEL -1	  // если X_CHANNEL ==0, то координатную зависимость не учитывать, беруться все координаты поперек счетчика
						  // если X_CHANNEL ==-1, то вообще выбросить это условие
#define on_off_one 0   // убрать +1 к регистру

#define on_off_aZ_conditions 0   // убрать условия на амплитуды стрипов	если ==0

#define MIN_CHANNEL 1	 // диапозон каналов для прогона в цикле

#define MAX_CHANNEL 16

const Int_t N_number = 1560000000; // 4884983 число событий из файла

const Int_t Smooth_index = 2;	// количество прогонов сглаживания гистограммы

const Double_t L_scint = 100.0 + 2.0 * (6.0 + 12.0); // размер сцинтиллятора + световоды с двух сторон, см, счетчики 120x120 мм

//const Double_t L_scint = 106.0 + 2.0 * (13.0 + 12.0); // размер сцинтиллятора + световоды с двух сторон, см, счетчики 200x200 мм

const Double_t L_scint_without = 106.0; // размер сцинтиллятора без световодов с двух сторон, см, счетчики 200x200 мм


//const Double_t L_scint = 106.0; // размер сцинтиллятора 200x200 мм

// const Double_t L_scint = 8.0; // ширина стрипа


const Double_t k_TDC = 0.1; // 1 канал ЗЦП=0,1 нс

const Double_t b1_1_pedestal = 44.0;	   // пьедисталы для вычета из амплитуд
const Double_t b1_2_pedestal = 34.0;
const Double_t b2_1_pedestal = 54.0;
const Double_t b2_2_pedestal = 68.0;
const Double_t b3_1_pedestal = 64.0;
const Double_t b3_2_pedestal = 45.0;

// диапозон по времени для фита

const Double_t t1_1_left = 1300.0 - 0 * 400.0;
const Double_t t1_1_right = 1600.0 - 0 * 400.0;

const Double_t t1_2_left = 1300.0 - 0 * 400.0;
const Double_t t1_2_right = 1600.0 - 0 * 400.0;

const Double_t t2_1_left = 1450.0 - 0 * 400.0;
const Double_t t2_1_right = 1650.0 - 0 * 400.0;

const Double_t t2_2_left = 1350.0 - 0 * 400.0;
const Double_t t2_2_right = 1500.0 - 0 * 400.0;

const Double_t t3_1_left = 1300.0 - 0 * 400.0;
const Double_t t3_1_right = 1500.0 - 0 * 400.0;

const Double_t t3_2_left = 1450.0 - 0 * 400.0;
const Double_t t3_2_right = 1700.0 - 0 * 400.0;



//// время-амплитудная поправка

///
const Double_t par0_tb1_1_left = 15000.0;
const Double_t par1_tb1_1_left = 1400.0 + 90.0 * 1 - 0 * 400.0 + 0 * 50.0 - 250.0 - 50.0;
const Double_t par2_tb1_1_left = b1_1_pedestal;

const Double_t par0_tb1_1_right = 15000.0;
const Double_t par1_tb1_1_right = 1400.0 + 90.0 * 1 - 0 * 400.0 + 0 * 50.0 - 250.0 - 50.0;
const Double_t par2_tb1_1_right = b1_1_pedestal;


const Double_t par0_tb1_2_left = 15000.0;
const Double_t par1_tb1_2_left = 1400.0 + 50.0 * 1 - 0 * 400.0 - 0 * 30.0 - 250.0 - 100.0 + 30.0 + 50.0 - 20.0;
const Double_t par2_tb1_2_left = b1_2_pedestal;

const Double_t par0_tb1_2_right = 15000.0;
const Double_t par1_tb1_2_right = 1400.0 + 50.0 * 1 - 0 * 400.0 - 0 * 30.0 - 250.0 - 100.0 + 30.0 + 50.0 - 20.0;
const Double_t par2_tb1_2_right = b1_2_pedestal;

///
const Double_t par0_tb2_1_left = 15000.0;
const Double_t par1_tb2_1_left = 1400.0 + 90.0 * 1 - 0 * 400.0 + 0 * 50.0 - 250.0 - 30.0;
const Double_t par2_tb2_1_left = b2_1_pedestal;

const Double_t par0_tb2_1_right = 15000.0;
const Double_t par1_tb2_1_right = 1400.0 + 90.0 * 1 - 0 * 400.0 + 0 * 50.0 - 250.0 - 30.0;
const Double_t par2_tb2_1_right = b2_1_pedestal;


const Double_t par0_tb2_2_left = 15000.0;
const Double_t par1_tb2_2_left = 1400.0 + 150.0 * 1 - 0 * 400.0 + 0 * 100.0 - 250.0 + 150.0 - 30.0;
const Double_t par2_tb2_2_left = b2_2_pedestal;

const Double_t par0_tb2_2_right = 15000.0;
const Double_t par1_tb2_2_right = 1400.0 + 150.0 * 1 - 0 * 400.0 + 0 * 100.0 - 250.0 + 150.0 - 30.0;
const Double_t par2_tb2_2_right = b2_2_pedestal;

///
const Double_t par0_tb3_1_left = 15000.0;
const Double_t par1_tb3_1_left = 1400.0 + 50.0 * 1 - 0 * 400.0 - 0 * 30.0 - 250.0 + 30.0 - 300.0 + 30.0;
const Double_t par2_tb3_1_left = b3_1_pedestal;

const Double_t par0_tb3_1_right = 15000.0;
const Double_t par1_tb3_1_right = 1400.0 + 50.0 * 1 - 0 * 400.0 - 0 * 30.0 - 250.0 + 30.0 - 300.0 + 30.0;
const Double_t par2_tb3_1_right = b3_1_pedestal;


const Double_t par0_tb3_2_left = 15000.0;
const Double_t par1_tb3_2_left = 1400.0 + 20.0 * 0 - 400.0 - 0 * 50.0 + 250.0 - 300.0 + 30.0;
const Double_t par2_tb3_2_left = b3_2_pedestal;

const Double_t par0_tb3_2_right = 15000.0;
const Double_t par1_tb3_2_right = 1400.0 + 20.0 * 0 - 400.0 - 0 * 50.0 + 250.0 - 300.0 + 30.0;
const Double_t par2_tb3_2_right = b3_2_pedestal;


/////
/////
///// Каты
/////
/////

// измерение 1: 410  330  500	   (счетчики 200x200)

 // измерение 2: 520  370  440	   (счетчики 200x200)

 // измерение 1: 400  360  290	   (счетчики 120x120)

 // измерение 2: 180  240  380	   (счетчики 120x120)

// измерение 3: 170  360  130	   (счетчики 120x120 и 200x200)

TCut cut_bb_1 = Form("sqrt((aF[0]-%f)*(aF[1]-%f))>170.0&&sqrt((aF[0]-%f)*(aF[1]-%f))<1900.0", 1.0 * b1_1_pedestal, 1.0 * b1_2_pedestal, 1.0 * b1_1_pedestal, 1.0 * b1_2_pedestal);
TCut cut_bb_2 = Form("sqrt((aF[2]-%f)*(aF[3]-%f))>360.0&&sqrt((aF[2]-%f)*(aF[3]-%f))<1900.0", 1.0 * b2_1_pedestal, 1.0 * b2_2_pedestal, 1.0 * b2_1_pedestal, 1.0 * b2_2_pedestal);
TCut cut_bb_3 = Form("sqrt((aF[4]-%f)*(aF[5]-%f))>130.0&&sqrt((aF[4]-%f)*(aF[5]-%f))<1900.0", 1.0 * b3_1_pedestal, 1.0 * b3_2_pedestal, 1.0 * b3_1_pedestal, 1.0 * b3_2_pedestal);

//TCut cut_bb_1 = "";
//TCut cut_bb_2 = "";
//TCut cut_bb_3 = "";



//const TCut cut_time_1 = "t1_1>900.0&&t1_2>900.0&&t1_1<1400.0&&t1_2<1400.0";
//const TCut cut_time_2 = "t2_1>1000.0&&t2_2>1000.0&&t2_1<1500.0&&t2_2<1500.0";
//const TCut cut_time_3 = "t3_1>1200.0&&t3_2>1550.0&&t3_1<2000.0&&t3_2<2000.0";

//const TCut cut_time_3 = "t3_2<1480";

//TCut cut_time_1 = "t1_1>1040.0&&t1_2>1070.0&&t1_1<1070.0&&t1_2<1100.0";
//TCut cut_time_2 = "t2_1>1070.0&&t2_2>1070.0&&t2_1<1120.0&&t2_2<1130.0";
//TCut cut_time_3 = "t3_1>1010.0&&t3_2>1020.0&&t3_1<1050.0&&t3_2<1060.0";


//Char_t cz[17];

//TCut *cz = new TCut[17];   // массив катов

//TCut cz[17];


/*


TCut cz2 = "(regZ==3)&&(16384-aZ[2]<13000)&&(16384-aZ[2]>2700)";
TCut cz3 = "(regZ==5)&&(16384-aZ[3]<13000)&&(16384-aZ[3]>2500)";
TCut cz4 = "(regZ==9)&&(16384-aZ[4]<13000)&&(16384-aZ[4]>2000)";
TCut cz5 = "(regZ==17)&&(16384-aZ[5]<13000)&&(16384-aZ[5]>2500)";
TCut cz6 = "(regZ==33)&&(16384-aZ[6]<13000)&&(16384-aZ[6]>2500)";
TCut cz7 = "(regZ==65)&&(16384-aZ[7]<16000)&&(16384-aZ[7]>4000)";
TCut cz8 = "(regZ==129)&&(16384-aZ[8]<16000)&&(16384-aZ[8]>4000)";
TCut cz9 = "(regZ==257)&&(16384-aZ[9]<13000)&&(16384-aZ[9]>600)";
TCut cz10 = "(regZ==513)&&(16384-aZ[10]<10000)&&(16384-aZ[10]>800)";
TCut cz11 = "(regZ==1025)&&(16384-aZ[11]<16000)&&(16384-aZ[11]>3500)";
TCut cz12 = "(regZ==2049)&&(16384-aZ[12]<16000)&&(16384-aZ[12]>4000)";
TCut cz13 = "(regZ==4097)&&(16384-aZ[13]<8000)&&(16384-aZ[13]>1000)";
TCut cz14 = "(regZ==8193)&&(16384-aZ[14]<16000)&&(16384-aZ[14]>4000)";
TCut cz15 = "(regZ==16385)&&(16384-aZ[15]<16000)&&(16384-aZ[15]>3500)";
TCut cz16 = "(regZ==32769)&&(16384-aZ[16]<16000)&&(16384-aZ[16]>3500)";

*/

/*



Char_t* cz2 = (Char_t*)"(regZ==3)&&(16384-aZ[2]<13000)&&(16384-aZ[2]>2700)";
Char_t* cz3 = (Char_t*)"(regZ==5)&&(16384-aZ[3]<13000)&&(16384-aZ[3]>2500)";
Char_t* cz4 = (Char_t*)"(regZ==9)&&(16384-aZ[4]<13000)&&(16384-aZ[4]>2000)";
Char_t* cz5 = (Char_t*)"(regZ==17)&&(16384-aZ[5]<13000)&&(16384-aZ[5]>2500)";
Char_t* cz6 = (Char_t*)"(regZ==33)&&(16384-aZ[6]<13000)&&(16384-aZ[6]>2500)";
Char_t* cz7 = (Char_t*)"(regZ==65)&&(16384-aZ[7]<16000)&&(16384-aZ[7]>4000)";
Char_t* cz8 = (Char_t*)"(regZ==129)&&(16384-aZ[8]<16000)&&(16384-aZ[8]>4000)";
Char_t* cz9 = (Char_t*)"(regZ==257)&&(16384-aZ[9]<13000)&&(16384-aZ[9]>600)";
Char_t* cz10 = (Char_t*)"(regZ==513)&&(16384-aZ[10]<10000)&&(16384-aZ[10]>800)";
Char_t* cz11 = (Char_t*)"(regZ==1025)&&(16384-aZ[11]<16000)&&(16384-aZ[11]>3500)";
Char_t* cz12 = (Char_t*)"(regZ==2049)&&(16384-aZ[12]<16000)&&(16384-aZ[12]>4000)";
Char_t* cz13 = (Char_t*)"(regZ==4097)&&(16384-aZ[13]<8000)&&(16384-aZ[13]>1000)";
Char_t* cz14 = (Char_t*)"(regZ==8193)&&(16384-aZ[14]<16000)&&(16384-aZ[14]>4000)";
Char_t* cz15 = (Char_t*)"(regZ==16385)&&(16384-aZ[15]<16000)&&(16384-aZ[15]>3500)";
Char_t* cz16 = (Char_t*)"(regZ==32769)&&(16384-aZ[16]<16000)&&(16384-aZ[16]>3500)";

*/


/*



cz[2] = "(regZ==3)&&(16384-aZ[2]<13000)&&(16384-aZ[2]>2700)";
cz[3] = "(regZ==5)&&(16384-aZ[3]<13000)&&(16384-aZ[3]>2500)";
cz[4] = "(regZ==9)&&(16384-aZ[4]<13000)&&(16384-aZ[4]>2000)";
cz[5] = "(regZ==17)&&(16384-aZ[5]<13000)&&(16384-aZ[5]>2500)";
cz[6] = "(regZ==33)&&(16384-aZ[6]<13000)&&(16384-aZ[6]>2500)";
cz[7] = "(regZ==65)&&(16384-aZ[7]<16000)&&(16384-aZ[7]>4000)";
cz[8] = "(regZ==129)&&(16384-aZ[8]<16000)&&(16384-aZ[8]>4000)";
cz[9] = "(regZ==257)&&(16384-aZ[9]<13000)&&(16384-aZ[9]>600)";
cz[10] = "(regZ==513)&&(16384-aZ[10]<10000)&&(16384-aZ[10]>800)";
cz[11] = "(regZ==1025)&&(16384-aZ[11]<16000)&&(16384-aZ[11]>3500)";
cz[12] = "(regZ==2049)&&(16384-aZ[12]<16000)&&(16384-aZ[12]>4000)";
cz[13] = "(regZ==4097)&&(16384-aZ[13]<8000)&&(16384-aZ[13]>1000)";
cz[14] = "(regZ==8193)&&(16384-aZ[14]<16000)&&(16384-aZ[14]>4000)";
cz[15] = "(regZ==16385)&&(16384-aZ[15]<16000)&&(16384-aZ[15]>3500)";
cz[16] = "(regZ==32769)&&(16384-aZ[16]<16000)&&(16384-aZ[16]>3500)";


*/
const TCut cut_time_3 = "";


const TCut cut_aplitude_1 = "aF[0]>10.0&&aF[1]>10.0&&aF[0]<1900.0&&aF[1]<1900.0";
const TCut cut_aplitude_2 = "aF[2]>10.0&&aF[3]>10.0&&aF[2]<1900.0&&aF[3]<1900.0";
const TCut cut_aplitude_3 = "aF[4]>10.0&&aF[5]>10.0&&aF[4]<1900.0&&aF[5]<1900.0";
const TCut cut_bb_pedestal = Form("aF[0]>%f && aF[1]>%f && aF[2]>%f && aF[3]>%f && aF[4]>%f && aF[5]>%f", b1_1_pedestal, b1_2_pedestal, b2_1_pedestal, b2_2_pedestal, b3_1_pedestal, b3_2_pedestal);

Int_t x_1_1 = 8;	   // ограничение по X калориметра для счетчика №1 (200x200)
Int_t x_1_2 = 9;
Int_t x_1_3 = 10;
Int_t x_1_4 = 11;
Int_t x_1_5 = 12;

Int_t x_2_1 = 5;	   // ограничение по X калориметра для счетчика №2 (200x200)
Int_t x_2_2 = 6;
Int_t x_2_3 = 7;
Int_t x_2_4 = 8;
Int_t x_2_5 = 9;

Int_t x_3_1 = 2;	   // ограничение по X калориметра для счетчика №3 (200x200)
Int_t x_3_2 = 3;
Int_t x_3_3 = 4;
Int_t x_3_4 = 5;					
Int_t x_3_5 = 6;

TCut cx_1 = Form("((regX==(pow(2, %i - 1)+1*%i))||(regX==(pow(2, %i - 1)+1*%i))||(regX==(pow(2, %i - 1)+1*%i))||(regX==(pow(2, %i - 1)+1*%i))||(regX==(pow(2, %i - 1)+1*%i))) && (regZ)", x_1_1, on_off_one, x_1_2, on_off_one, x_1_3, on_off_one, x_1_4, on_off_one, x_1_5, on_off_one);

TCut cx_2 = Form("((regX==(pow(2, %i - 1)+1*%i))||(regX==(pow(2, %i - 1)+1*%i))||(regX==(pow(2, %i - 1)+1*%i))||(regX==(pow(2, %i - 1)+1*%i))||(regX==(pow(2, %i - 1)+1*%i))) && (regZ)", x_2_1, on_off_one, x_2_2, on_off_one, x_2_3, on_off_one, x_2_4, on_off_one, x_2_5, on_off_one);

TCut cx_3 = Form("((regX==(pow(2, %i - 1)+1*%i))||(regX==(pow(2, %i - 1)+1*%i))||(regX==(pow(2, %i - 1)+1*%i))||(regX==(pow(2, %i - 1)+1*%i))||(regX==(pow(2, %i - 1)+1*%i))) && (regZ)", x_3_1, on_off_one, x_3_2, on_off_one, x_3_3, on_off_one, x_3_4, on_off_one, x_3_5, on_off_one);



//const TCut reg_HHC = Form("(regZ&pow(2, %i)) && (regX)", Z_CHANNEL);

TCut reg_HHC;

TCut cut_full = "";

TCut cut_full_1 = "";

TCut cut_full_2 = "";

TCut cut_full_3 = "";

TCut* cz = new TCut[17];

void name_cut() {




	/*


	cz[2] = "(regZ==3)&&(16384-aZ[2]<13000)&&(16384-aZ[2]>2700)";
	cz[3] = "(regZ==5)&&(16384-aZ[3]<13000)&&(16384-aZ[3]>2500)";
	cz[4] = "(regZ==9)&&(16384-aZ[4]<13000)&&(16384-aZ[4]>2000)";
	cz[5] = "(regZ==17)&&(16384-aZ[5]<13000)&&(16384-aZ[5]>2500)";
	cz[6] = "(regZ==33)&&(16384-aZ[6]<13000)&&(16384-aZ[6]>2500)";
	cz[7] = "(regZ==65)&&(16384-aZ[7]<16000)&&(16384-aZ[7]>4000)";
	cz[8] = "(regZ==129)&&(16384-aZ[8]<16000)&&(16384-aZ[8]>4000)";
	cz[9] = "(regZ==257)&&(16384-aZ[9]<13000)&&(16384-aZ[9]>600)";
	cz[10] = "(regZ==513)&&(16384-aZ[10]<10000)&&(16384-aZ[10]>800)";
	cz[11] = "(regZ==1025)&&(16384-aZ[11]<16000)&&(16384-aZ[11]>3500)";
	cz[12] = "(regZ==2049)&&(16384-aZ[12]<16000)&&(16384-aZ[12]>4000)";
	cz[13] = "(regZ==4097)&&(16384-aZ[13]<8000)&&(16384-aZ[13]>1000)";
	cz[14] = "(regZ==8193)&&(16384-aZ[14]<16000)&&(16384-aZ[14]>4000)";
	cz[15] = "(regZ==16385)&&(16384-aZ[15]<16000)&&(16384-aZ[15]>3500)";
	cz[16] = "(regZ==32769)&&(16384-aZ[16]<16000)&&(16384-aZ[16]>3500)";
  */


	if (on_off_aZ_conditions == 0) {

		cz[1] = "";
		cz[2] = "";
		cz[3] = "";
		cz[4] = "";
		cz[5] = "";
		cz[6] = "";
		cz[7] = "";
		cz[8] = "";
		cz[9] = "";
		cz[10] = "";
		cz[11] = "";
		cz[12] = "";
		cz[13] = "";
		cz[14] = "";
		cz[15] = "";
		cz[16] = "";
	}
	else {
		cz[1] = "";
		cz[2] = "(16384-aZ[2]<13000)&&(16384-aZ[2]>2700)";
		cz[3] = "(16384-aZ[3]<13000)&&(16384-aZ[3]>2500)";
		cz[4] = "(16384-aZ[4]<13000)&&(16384-aZ[4]>2000)";
		cz[5] = "(16384-aZ[5]<13000)&&(16384-aZ[5]>2500)";
		cz[6] = "(16384-aZ[6]<13000)&&(16384-aZ[6]>2500)";
		cz[7] = "(16384-aZ[7]<16000)&&(16384-aZ[7]>4000)";
		cz[8] = "(16384-aZ[8]<16000)&&(16384-aZ[8]>4000)";
		cz[9] = "(16384-aZ[9]<13000)&&(16384-aZ[9]>600)";
		cz[10] = "(16384-aZ[10]<10000)&&(16384-aZ[10]>800)";
		cz[11] = "(16384-aZ[11]<16000)&&(16384-aZ[11]>3500)";
		cz[12] = "(16384-aZ[12]<16000)&&(16384-aZ[12]>4000)";
		cz[13] = "(16384-aZ[13]<8000)&&(16384-aZ[13]>850)";
		cz[14] = "(16384-aZ[14]<16000)&&(16384-aZ[14]>4000)";
		cz[15] = "(16384-aZ[15]<16000)&&(16384-aZ[15]>3500)";
		cz[16] = "(16384-aZ[16]<16000)&&(16384-aZ[16]>3500)";

	}

	


	Char_t* name;

	name = Form("cz%i", Z_CHANNEL);

	//cout << name << endl;

	//TCut reg_HHC = Form("(regZ==(pow(2, %i - 1)+1)) && (regX) && %s", Z_CHANNEL, name);

	TCut reg_HHC_temp = "";

	if (X_CHANNEL == -1) {
		reg_HHC_temp = Form("regZ==(pow(2, %i - 1)+1*%i)", Z_CHANNEL, on_off_one);
	}
	else {

		reg_HHC_temp = Form("(regZ==(pow(2, %i - 1)+1*%i)) && (regX)", Z_CHANNEL, on_off_one);
		//reg_HHC_temp = Form("(regZ==(pow(2, %i - 1)+1*%i))", Z_CHANNEL, on_off_one);

		
	}

	reg_HHC = reg_HHC_temp && cz[Z_CHANNEL];


	if (Z_CHANNEL == 0) {
		reg_HHC = "";
	}

	if (X_CHANNEL == 0 || X_CHANNEL == -1) {
		cx_1 = "";
		cx_2 = "";
		cx_3 = "";
	}

	cut_full_1 = cut_bb_1 && cut_aplitude_1 && cx_1 && reg_HHC;

	cut_full_2 = cut_bb_2 && cut_aplitude_2 && cx_2 && reg_HHC;

	cut_full_3 = cut_bb_3 && cut_aplitude_3 && cx_3 && reg_HHC;
}






//const TCut reg_HHC = "";

// const TCut cut_full = cut_aplitude_1 && cut_aplitude_2 && cut_aplitude_3 && cut_time_1 && cut_time_2 && cut_time_3 && cut_bb_1 && cut_bb_2 && cut_bb_3 && cut_bb_pedestal;

//	TCut cut_full = cut_aplitude_1 && cut_aplitude_2 && cut_aplitude_3 && cut_bb_1 && cut_bb_2 && cut_bb_3 && cut_bb_pedestal;

// const TCut cut_full = cut_bb_1 && cut_bb_2 && cut_bb_3 && cut_time_3 && reg_HHC;


// const TCut cut_full = "regX && regZ";








//	TCut cut_full = cut_bb_1 && cut_bb_2;

//	TCut cut_full = cut_bb_3;

//const TCut cut_b1b2 = "b2_1>80.0&&b2_2>110.0";

//const TCut cut_full = cut_bb_1&&cut_bb_3&&cut_b1b2;

void zero_fit_par() {


}

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


	Double_t fitval = par[0] / (x[0] - par[1]) + par[2];

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

// TFile* f = TFile::Open("../macros/stend_07.04.21.root"); // со смазкой, счетчик 6 в центре	    3 - сверху, 1 - снизу   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс), вариант триггера №2

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

// TFile* f = TFile::Open("../macros/tqmu_19.05.21.root"); // со смазкой, счетчик №1 счетчик 1, №2 счетчик 3, №3 счетчик 6   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс)

//TFile* f = TFile::Open("../macros/tqmu_21.05.21.root"); // со смазкой, счетчик №1 счетчик 5, №2 счетчик 4, №3 счетчик 2   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс)

// ТЕСТИРОВАНИЕ 8 СЧЕТЧИКОВ 120x120 мм + калориметр для измерения по зонам

//TFile* f = TFile::Open("../macros/tqmu_25.05.21.root"); // со смазкой, счетчик №1 счетчик 8, №2 счетчик 1, №3 счетчик 3   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс)

//TFile* f = TFile::Open("../macros/tqmu_26.05.21.root"); // со смазкой, счетчик №1 счетчик 8, №2 счетчик 1, №3 счетчик 3   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс)	передвинули к правому краю, НЕПРАВИЛЬНО

//TFile* f = TFile::Open("../macros/tqmu_27.05.21.root"); // со смазкой, счетчик №1 счетчик 8, №2 счетчик 1, №3 счетчик 3   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс) передвинули к левому краю, ПРАВИЛЬНО

// TFile* f = TFile::Open("../macros/tqmu_28.05.21.root"); // со смазкой, счетчик №1 счетчик 4, №2 счетчик 2, №3 счетчик 7   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс) передвинули к левому краю, ПРАВИЛЬНО  (возможго 12, 13, 14 каналы по Z как-то перепутаны ??)

TFile* f = TFile::Open("../macros/tqmu_31.05.21.root"); // со смазкой, счетчик №1 счетчик 6 (120x120), №2 счетчик 5 (120x120), №3 счетчик 3 (200x200)   ПОМЕНЯЛИ БЛОК ПИТАНИЯ ВВИ, поменяли дискриминатор на другой (разрешение электроники 100 пс) передвинули к левому краю, ПРАВИЛЬНО  (возможго 12, 13, 14 каналы по Z как-то перепутаны ??)



TTree* ntpl = (TTree*)f->Get("T");

TTree* ntpl_select = (TTree*)f->Get("T");

//TTree* ntpl_select = ntpl->CopyTree(cut_full, "", N_number);	  // создает новое дерево с условиями отбора


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



void time(Int_t number_count) {

	name_cut();

	
	TH1F* h_temp = new TH1F("h_temp", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);

	
	//TCut* cz = new TCut[17];

	/*
	
	

	TH1F** h_temp = new TH1F*[17];


	h_temp[5] = new TH1F("h_temp5", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	h_temp[6] = new TH1F("h_temp6", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	h_temp[7] = new TH1F("h_temp7", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	h_temp[8] = new TH1F("h_temp8", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	h_temp[9] = new TH1F("h_temp9", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	h_temp[10] = new TH1F("h_temp10", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	h_temp[11] = new TH1F("h_temp11", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	h_temp[12] = new TH1F("h_temp12", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	h_temp[13] = new TH1F("h_temp13", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	h_temp[14] = new TH1F("h_temp14", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	h_temp[15] = new TH1F("h_temp15", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	h_temp[16] = new TH1F("h_temp16", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);

   */

	/*



	if (number_count == 1) {
		TH1F* h_temp = new TH1F("h_temp", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);

	}
	else
		if (number_count == 2) {
			TH1F* h_temp = new TH1F("h_temp", "0.5*(t2_1-t2_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);

		}
		else	if (number_count == 3) {
			TH1F* h_temp = new TH1F("h_temp", "0.5*(t3_1-t3_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);

		}

   */



	Double_t sigma_temp_ = 0.0;


	/*
	
	

	TF1** fit_temp = new TF1 * [17];
	
	fit_temp[5] = new TF1("fit_temp5", "gaus", -100, 100);
	fit_temp[6] = new TF1("fit_temp6", "gaus", -100, 100);
	fit_temp[7] = new TF1("fit_temp7", "gaus", -100, 100);
	fit_temp[8] = new TF1("fit_temp8", "gaus", -100, 100);
	fit_temp[9] = new TF1("fit_temp9", "gaus", -100, 100);
	fit_temp[10] = new TF1("fit_temp10", "gaus", -100, 100);
	fit_temp[11] = new TF1("fit_temp11", "gaus", -100, 100);
	fit_temp[12] = new TF1("fit_temp12", "gaus", -100, 100);
	fit_temp[13] = new TF1("fit_temp13", "gaus", -100, 100);
	fit_temp[14] = new TF1("fit_temp14", "gaus", -100, 100);
	fit_temp[15] = new TF1("fit_temp15", "gaus", -100, 100);
	fit_temp[16] = new TF1("fit_temp16", "gaus", -100, 100);

	*/
	
	TF1* fit_temp = new TF1("fit_temp", "gaus", -100, 100);


	TCut reg_HHC_temp = "";

	TCut reg_HHC_temp_temp = "";

	TCut cut_full_temp = "";

	Double_t sigma_arr[16];

//	h_temp->GetXaxis()->SetTitle("Time, channel number");
//	h_temp->GetYaxis()->SetTitle("Events number");
//	h_temp->GetXaxis()->CenterTitle();
//	h_temp->GetYaxis()->CenterTitle();

	TCanvas* MyC_time = new TCanvas("MyC_time", "MyC_time", 0, 0, 2000, 1000);	  // для определения числа ф.э.
	//TString vart_temp = "";
	
	
	
	
	TString vart_temp = Form("0.5*((tF[0]-delta_t(aF[0], %f, %f, %f))-(tF[1]-delta_t(aF[1], %f, %f, %f)))>>h_temp", par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right);


	if (number_count == 1) {
		vart_temp = Form("0.5*((tF[0]-delta_t(aF[0], %f, %f, %f))-(tF[1]-delta_t(aF[1], %f, %f, %f)))>>h_temp", par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right);

	}
	else if (number_count == 2) {
		vart_temp = Form("0.5*((tF[2]-delta_t(aF[2], %f, %f, %f))-(tF[3]-delta_t(aF[3], %f, %f, %f)))>>h_temp", par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);

	}
	else if (number_count == 3) {
		vart_temp = Form("0.5*((tF[4]-delta_t(aF[4], %f, %f, %f))-(tF[5]-delta_t(aF[5], %f, %f, %f)))>>h_temp", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right);

	}
   
	//cout << vart_temp << endl;

	MyC_time->Divide(1, 1);

	MyC_time->cd(1);

	for (Int_t i = MIN_CHANNEL; i <= MAX_CHANNEL; i = i + 1) {

		cout << endl;
		cout << i << endl;
		cout << endl;

		reg_HHC_temp = "";

		reg_HHC_temp_temp = "";

		cut_full_temp = "";


		//h_temp->Reset();

		//h_temp->Clear();

	//	h_temp = new TH1F("h_temp", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);



		

		if (X_CHANNEL == -1) {
			reg_HHC_temp_temp = Form("regZ==(pow(2, %i - 1)+1*%i)", i, on_off_one);
		}
		else {

			reg_HHC_temp_temp = Form("(regZ==(pow(2, %i - 1)+1*%i)) && (regX)", i, on_off_one);
			//reg_HHC_temp_temp = Form("(regZ==(pow(2, %i - 1)+1*%i))", i, on_off_one);


		}


		reg_HHC_temp = reg_HHC_temp_temp && cz[i];



		if (number_count == 1) {
			cut_full_temp = cut_bb_1 && cut_aplitude_1 && cx_1 && reg_HHC_temp;
			//vart_temp = Form("0.5*((tF[0]-delta_t(aF[0], %f, %f, %f))-(tF[1]-delta_t(aF[1], %f, %f, %f)))>>h_temp[%i]", par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right, i);

		}
		else if (number_count == 2) {
			cut_full_temp = cut_bb_2 && cut_aplitude_2 && cx_2 && reg_HHC_temp;
			//vart_temp = Form("0.5*((tF[2]-delta_t(aF[2], %f, %f, %f))-(tF[3]-delta_t(aF[3], %f, %f, %f)))>>h_temp[%i]", par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right, i);

		}
		else if (number_count == 3) {
			cut_full_temp = cut_bb_3 && cut_aplitude_3 && cx_3 && reg_HHC_temp;
			//vart_temp = Form("0.5*((tF[4]-delta_t(aF[4], %f, %f, %f))-(tF[5]-delta_t(aF[5], %f, %f, %f)))>>h_temp[%i]", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right, i);

		}





		//	cut_full_temp = cut_bb_1 && cut_aplitude_1 && reg_HHC_temp;

		ntpl_select->Draw(vart_temp, cut_full_temp);

		h_temp->Smooth(Smooth_index);	 // Сглаживание

		gStyle->SetOptFit(1111);
		h_temp->Fit(fit_temp, "", "", -100, 100);

		sigma_temp_ = fit_temp->GetParameter(2);
		sigma_arr[i] = sigma_temp_;

		h_temp->Draw();

		//h_temp->Reset();

		//h_temp->Clear();

		//TF1* fit_temp = new TF1("fit_temp", "gaus", -100, 100);


	}

	for (Int_t i = MIN_CHANNEL; i <= MAX_CHANNEL; i = i + 1) {

		cout << sigma_arr[i] * k_TDC << endl;	  // в нс	


	}

	MyC_time->Update();

}



void time_(Int_t number_count, Int_t i) {

	name_cut();

	TH1F* h_temp = new TH1F("h_temp", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);


	/*



	  if (number_count == 1) {
		  TH1F* h_temp = new TH1F("h_temp", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	  }
	  else
		  if (number_count == 2) {
			  TH1F* h_temp = new TH1F("h_temp", "0.5*(t2_1-t2_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
		  }
		  else	if (number_count == 3) {
			  TH1F* h_temp = new TH1F("h_temp", "0.5*(t3_1-t3_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);
		  }


	 */

	Double_t sigma_temp_ = 0.0;

	TF1* fit_temp = new TF1("fit_temp", "gaus", -100, 100);


	TCut reg_HHC_temp = "";

	TCut reg_HHC_temp_temp = "";

	TCut cut_full_temp = "";

	Double_t sigma_arr[17];

	h_temp->GetXaxis()->SetTitle("Time, channel number");
	h_temp->GetYaxis()->SetTitle("Events number");
	h_temp->GetXaxis()->CenterTitle();
	h_temp->GetYaxis()->CenterTitle();

	TCanvas* MyC_time = new TCanvas("MyC_time", "MyC_time", 0, 0, 2000, 1000);	  // для определения числа ф.э.

	TString vart_temp = Form("0.5*((tF[0]-delta_t(aF[0], %f, %f, %f))-(tF[1]-delta_t(aF[1], %f, %f, %f)))>>h_temp", par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right);


	if (number_count == 1) {
		vart_temp = Form("0.5*((tF[0]-delta_t(aF[0], %f, %f, %f))-(tF[1]-delta_t(aF[1], %f, %f, %f)))>>h_temp", par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right);
		//	cout << 11111 << endl;
	}
	else if (number_count == 2) {
		vart_temp = Form("0.5*((tF[2]-delta_t(aF[2], %f, %f, %f))-(tF[3]-delta_t(aF[3], %f, %f, %f)))>>h_temp", par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);

	}
	else if (number_count == 3) {
		vart_temp = Form("0.5*((tF[4]-delta_t(aF[4], %f, %f, %f))-(tF[5]-delta_t(aF[5], %f, %f, %f)))>>h_temp", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right);

	}


	//cout << vart_temp << endl;

	MyC_time->Divide(1, 1);


	MyC_time->cd(1);


	//	h_temp->Reset();

	//	h_temp->Clear();

	//	h_temp = new TH1F("h_temp", "0.5*(t1_1-t1_2) WITH AMPLITUDE CORRECTION", 800, -400.0, 400.0);


	

	if (X_CHANNEL == -1) {
		reg_HHC_temp_temp = Form("regZ==(pow(2, %i - 1)+1*%i)", i, on_off_one);
	}
	else {

		reg_HHC_temp_temp = Form("(regZ==(pow(2, %i - 1)+1*%i)) && (regX)", i, on_off_one);
	}


	reg_HHC_temp = reg_HHC_temp_temp && cz[i];

	//		cut_full_temp = cut_bb_1 && cut_aplitude_1 && cx_1 && reg_HHC_temp;


	if (number_count == 1) {
		cut_full_temp = cut_bb_1 && cut_aplitude_1 && cx_1 && reg_HHC_temp;
		//		cout << 1111 << endl;
	}
	else if (number_count == 2) {
		cut_full_temp = cut_bb_2 && cut_aplitude_2 && cx_2 && reg_HHC_temp;

	}
	else if (number_count == 3) {
		cut_full_temp = cut_bb_3 && cut_aplitude_3 && cx_3 && reg_HHC_temp;

	}





	ntpl_select->Draw(vart_temp, cut_full_temp);

	h_temp->Smooth(Smooth_index);	 // Сглаживание

	gStyle->SetOptFit(1111);
	h_temp->Fit(fit_temp, "", "", -100, 100);

	sigma_temp_ = fit_temp->GetParameter(2);
	sigma_arr[i] = sigma_temp_;

	h_temp->Draw();


	cout << sigma_arr[i] * k_TDC << endl;	  // в нс	

	MyC_time->Update();

}





// амплитуда-амплитуда
void amplit(Int_t number_count) {

	name_cut();

	TH1F* h1_b1b2_1_temp1 = new TH1F("h1_b1b2_1_temp1", "0.5*ln(A1_1/A1_2) -pedestal", 200, -2.0, 2.0);


	//	TH1F* h1_b1b2_2 = new TH1F("h1_b1b2_2", "sqrt(A1_1*A1_2) -pedestal", 1500, 0, 1500.0);



	//	TH2* h4_1 = new TH2D("h4_1", "b1_1:b1_2 -pedestal", 2100, 0, 2100, 2100, 0, 2100);

	//	TH2* h4_1_same = new TH2D("h4_1_same", "b1_1:b1_2+Cut -pedestal", 2100, 0, 2100, 2100, 0, 2100);


		// set the parameters to the mean and RMS of the histogram
	func_1->SetParameters(10.0, 2.0, -2.0, 2.0, -0.5);
	func_2->SetParameters(10.0, 1.0, -100.0, 1.0, -100.0);
	// give the parameters meaningful names
	func_1->SetParNames("A", "S1", "X1", "S2", "X2");
	func_2->SetParNames("A", "S1", "X1", "S2", "X2");


	func_3->SetParameters(10000.0, 1000.0, 10.0);
	func_3->SetParNames("k", "t0", "b0");


	TCanvas* MyC_amplit_ = new TCanvas("MyC_amplit_", "MyC_amplit_", 0, 0, 1000, 1000);



	TCut reg_HHC_temp = "";

	TCut reg_HHC_temp_temp = "";

	TCut cut_full_temp = "";

	Double_t lambda_arr[17];

	Double_t k_lambda_arr[17];

	Int_t n_number_arr[17];	 // массив для определения засветки

	Int_t n_number_temp = 0;

	//TString var1_temp = Form("(aF[0] - %f) : (aF[1] - %f) >> h4_1", b1_1_pedestal, b1_2_pedestal);
	//TString var1_same_temp = Form("(aF[0] - %f) : (aF[1] - %f) >> h4_1_same", b1_1_pedestal, b1_2_pedestal);

	//TString var2_temp = Form("sqrt((aF[0] - %f)*(aF[1] - %f)) >> h1_b1b2_2", b1_1_pedestal, b1_2_pedestal);

	TString var3_temp = Form("0.5*log((aF[0] - %f)/(aF[1] - %f)) >> h1_b1b2_1_temp1", b1_1_pedestal, b1_2_pedestal);




	if (number_count == 1) {
		//	var1_temp = Form("(aF[0] - %f) : (aF[1] - %f) >> h4_1", b1_1_pedestal, b1_2_pedestal);
		//	var1_same_temp = Form("(aF[0] - %f) : (aF[1] - %f) >> h4_1_same", b1_1_pedestal, b1_2_pedestal);

		//	var2_temp = Form("sqrt((aF[0] - %f)*(aF[1] - %f)) >> h1_b1b2_2", b1_1_pedestal, b1_2_pedestal);

		var3_temp = Form("0.5*log((aF[0] - %f)/(aF[1] - %f)) >> h1_b1b2_1_temp1", b1_1_pedestal, b1_2_pedestal);
	}
	else if (number_count == 2) {
		//	var1_temp = Form("(aF[2] - %f) : (aF[3] - %f) >> h4_1", b2_1_pedestal, b2_2_pedestal);
		//	var1_same_temp = Form("(aF[2] - %f) : (aF[3] - %f) >> h4_1_same", b2_1_pedestal, b2_2_pedestal);

		//	var2_temp = Form("sqrt((aF[2] - %f)*(aF[3] - %f)) >> h1_b1b2_2", b2_1_pedestal, b2_2_pedestal);

		var3_temp = Form("0.5*log((aF[2] - %f)/(aF[3] - %f)) >> h1_b1b2_1_temp1", b2_1_pedestal, b2_2_pedestal);
	}
	else if (number_count == 3) {
		//	var1_temp = Form("(aF[4] - %f) : (aF[5] - %f) >> h4_1", b3_1_pedestal, b3_2_pedestal);
		//	var1_same_temp = Form("(aF[4] - %f) : (aF[5] - %f) >> h4_1_same", b3_1_pedestal, b3_2_pedestal);

		//	var2_temp = Form("sqrt((aF[4] - %f)*(aF[5] - %f)) >> h1_b1b2_2", b3_1_pedestal, b3_2_pedestal);

		var3_temp = Form("0.5*log((aF[4] - %f)/(aF[5] - %f)) >> h1_b1b2_1_temp1", b3_1_pedestal, b3_2_pedestal);

	}



	MyC_amplit_->Divide(1, 3);

	/*



	MyC_amplit_->cd(1);




	ntpl->Draw(var1_temp, "", "", N_number);
	ntpl_select->Draw(var1_same_temp, cut_full_temp);

	h4_1->GetXaxis()->SetTitle("Amplitude, channel number");
	h4_1->GetYaxis()->SetTitle("Amplitude, channel number");
	h4_1->GetXaxis()->CenterTitle();
	h4_1->GetYaxis()->CenterTitle();


	h4_1->Draw();
	h4_1_same->SetMarkerColor(kRed);
	h4_1_same->Draw("same");



	MyC_amplit_->cd(2);


	h1_b1b2_2->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(var2_temp, cut_full_temp);

	gStyle->SetOptFit(1111);
	h1_b1b2_2->Fit("landau", "", "", 100.0, 800.0);

	h1_b1b2_2->SetFillColor(kGreen);
	h1_b1b2_2->Draw();
   */

	MyC_amplit_->cd(3);

	for (Int_t i = MIN_CHANNEL; i <= MAX_CHANNEL; i = i + 1) {

		cout << endl;
		cout << i << endl;
		cout << endl;

		reg_HHC_temp_temp = "";
		reg_HHC_temp = "";
		cut_full_temp = "";

		

		if (X_CHANNEL == -1) {
			reg_HHC_temp_temp = Form("regZ==(pow(2, %i - 1)+1*%i)", i, on_off_one);
		}
		else {

			reg_HHC_temp_temp = Form("(regZ==(pow(2, %i - 1)+1*%i)) && (regX)", i, on_off_one);
			//reg_HHC_temp_temp = Form("(regZ==(pow(2, %i - 1)+1*%i))", i, on_off_one);


		}

		reg_HHC_temp = reg_HHC_temp_temp && cz[i];

		//cut_full_temp = cut_bb_1 && cut_aplitude_1 && cx_1 && reg_HHC_temp;


		if (number_count == 1) {
			cut_full_temp = cut_bb_1 && cut_aplitude_1 && cx_1 && reg_HHC_temp;
		}
		else if (number_count == 2) {
			cut_full_temp = cut_bb_2 && cut_aplitude_2 && cx_2 && reg_HHC_temp;
		}
		else if (number_count == 3) {
			cut_full_temp = cut_bb_3 && cut_aplitude_3 && cx_3 && reg_HHC_temp;

		}

		h1_b1b2_1_temp1->Smooth(Smooth_index);	 // Сглаживание


		n_number_temp = ntpl_select->Draw(var3_temp, cut_full_temp);

		gStyle->SetOptFit(1111);
		
		func_1->SetParameters(10.0, 2.0, -2.0, 2.0, -0.5);

		
		h1_b1b2_1_temp1->Fit("my_fit_1", "", "", -2, 2);
		h1_b1b2_1_temp1->Fit("my_fit_1", "", "", -2, 2);
		h1_b1b2_1_temp1->Draw();

		
	


		S1 = func_1->GetParameter(1);
		X1 = func_1->GetParameter(2);
		S2 = func_1->GetParameter(3);
		X2 = func_1->GetParameter(4);

		Double_t lambda = L_scint / delta_L(S1, X1, S2, X2);
		Double_t k_lambda = exp(L_scint_without / lambda);

		printf("lambda_scint= %.3f cm\n", lambda);
		printf("oslablenye= %.3f raz\n", k_lambda);
		printf("number of events= %i \n", (int) ntpl_select->GetEntries(cut_full_temp));

		lambda_arr[i] = lambda;

		k_lambda_arr[i] = k_lambda;

		n_number_arr[i] = n_number_temp;
		
		

	//	h1_b1b2_1->Reset();

	//	h1_b1b2_1->Clear();

	}

	for (Int_t i = MIN_CHANNEL; i <= MAX_CHANNEL; i = i + 1) {

		cout << lambda_arr[i] << endl;	  // в см	


	}

	cout << endl;

	for (Int_t i = MIN_CHANNEL; i <= MAX_CHANNEL; i = i + 1) {

		cout << k_lambda_arr[i] << endl;	  // в см	


	}

	cout << endl;

	for (Int_t i = MIN_CHANNEL; i <= MAX_CHANNEL; i = i + 1) {

		cout << n_number_arr[i] << endl;	  // в см	


	}


	MyC_amplit_->Update();

}


// амплитуда-амплитуда
void amplit_(Int_t number_count, Int_t i) {

	name_cut();

	TH1F* h1_b1b2_1_temp2 = new TH1F("h1_b1b2_1_temp2", "0.5*ln(A1_1/A1_2) -pedestal", 200, -2.0, 2.0);


	TH1F* h1_b1b2_2 = new TH1F("h1_b1b2_2", "sqrt(A1_1*A1_2) -pedestal", 1500, 0, 1500.0);



	TH2* h4_1 = new TH2D("h4_1", "b1_1:b1_2 -pedestal", 2100, 0, 2100, 2100, 0, 2100);

	TH2* h4_1_same = new TH2D("h4_1_same", "b1_1:b1_2+Cut -pedestal", 2100, 0, 2100, 2100, 0, 2100);


	// set the parameters to the mean and RMS of the histogram
	func_1->SetParameters(10.0, 2.0, -2.0, 2.0, -0.5);
	func_2->SetParameters(10.0, 1.0, -100.0, 1.0, -100.0);
	// give the parameters meaningful names
	func_1->SetParNames("A", "S1", "X1", "S2", "X2");
	func_2->SetParNames("A", "S1", "X1", "S2", "X2");


	func_3->SetParameters(10000.0, 1000.0, 10.0);
	func_3->SetParNames("k", "t0", "b0");


	TCanvas* MyC_amplit_ = new TCanvas("MyC_amplit_", "MyC_amplit_", 0, 0, 2000, 1000);



	TCut reg_HHC_temp = "";

	TCut reg_HHC_temp_temp = "";

	TCut cut_full_temp = "";

	Double_t lambda_arr[17];

	Double_t k_lambda_arr[17];

	Int_t n_number_temp = 0;



	TString var1_temp = Form("(aF[0] - %f) : (aF[1] - %f) >> h4_1", b1_1_pedestal, b1_2_pedestal);
	TString var1_same_temp = Form("(aF[0] - %f) : (aF[1] - %f) >> h4_1_same", b1_1_pedestal, b1_2_pedestal);

	TString var2_temp = Form("sqrt((aF[0] - %f)*(aF[1] - %f)) >> h1_b1b2_2", b1_1_pedestal, b1_2_pedestal);

	TString var3_temp = Form("0.5*log((aF[0] - %f)/(aF[1] - %f)) >> h1_b1b2_1", b1_1_pedestal, b1_2_pedestal);




	if (number_count == 1) {
		var1_temp = Form("(aF[0] - %f) : (aF[1] - %f) >> h4_1", b1_1_pedestal, b1_2_pedestal);
		var1_same_temp = Form("(aF[0] - %f) : (aF[1] - %f) >> h4_1_same", b1_1_pedestal, b1_2_pedestal);

		var2_temp = Form("sqrt((aF[0] - %f)*(aF[1] - %f)) >> h1_b1b2_2", b1_1_pedestal, b1_2_pedestal);

		var3_temp = Form("0.5*log((aF[0] - %f)/(aF[1] - %f)) >> h1_b1b2_1_temp2", b1_1_pedestal, b1_2_pedestal);
	}
	else if (number_count == 2) {
		var1_temp = Form("(aF[2] - %f) : (aF[3] - %f) >> h4_1", b2_1_pedestal, b2_2_pedestal);
		var1_same_temp = Form("(aF[2] - %f) : (aF[3] - %f) >> h4_1_same", b2_1_pedestal, b2_2_pedestal);

		var2_temp = Form("sqrt((aF[2] - %f)*(aF[3] - %f)) >> h1_b1b2_2", b2_1_pedestal, b2_2_pedestal);

		var3_temp = Form("0.5*log((aF[2] - %f)/(aF[3] - %f)) >> h1_b1b2_1_temp2", b2_1_pedestal, b2_2_pedestal);
	}
	else if (number_count == 3) {
		var1_temp = Form("(aF[4] - %f) : (aF[5] - %f) >> h4_1", b3_1_pedestal, b3_2_pedestal);
		var1_same_temp = Form("(aF[4] - %f) : (aF[5] - %f) >> h4_1_same", b3_1_pedestal, b3_2_pedestal);

		var2_temp = Form("sqrt((aF[4] - %f)*(aF[5] - %f)) >> h1_b1b2_2", b3_1_pedestal, b3_2_pedestal);

		var3_temp = Form("0.5*log((aF[4] - %f)/(aF[5] - %f)) >> h1_b1b2_1_temp2", b3_1_pedestal, b3_2_pedestal);

	}




	if (X_CHANNEL == -1) {
		reg_HHC_temp_temp = Form("regZ==(pow(2, %i - 1)+1*%i)", i, on_off_one);
	}
	else {

	  	reg_HHC_temp_temp = Form("(regZ==(pow(2, %i - 1)+1*%i)) && (regX)", i, on_off_one);


	}


	reg_HHC_temp = reg_HHC_temp_temp && cz[i];

	//cut_full_temp = cut_bb_1 && cut_aplitude_1 && cx_1 && reg_HHC_temp;


	if (number_count == 1) {
		cut_full_temp = cut_bb_1 && cut_aplitude_1 && cx_1 && reg_HHC_temp;
	}
	else if (number_count == 2) {
		cut_full_temp = cut_bb_2 && cut_aplitude_2 && cx_2 && reg_HHC_temp;
	}
	else if (number_count == 3) {
		cut_full_temp = cut_bb_3 && cut_aplitude_3 && cx_3 && reg_HHC_temp;

	}

	MyC_amplit_->Divide(1, 3);

	MyC_amplit_->cd(1);




	ntpl->Draw(var1_temp, "", "", N_number);
	 
	
	ntpl_select->Draw(var1_same_temp, cut_full_temp);

	h4_1->GetXaxis()->SetTitle("Amplitude, channel number");
	h4_1->GetYaxis()->SetTitle("Amplitude, channel number");
	h4_1->GetXaxis()->CenterTitle();
	h4_1->GetYaxis()->CenterTitle();


	h4_1->Draw();
	h4_1_same->SetMarkerColor(kRed);
	h4_1_same->Draw("same");



	MyC_amplit_->cd(2);


	h1_b1b2_2->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(var2_temp, cut_full_temp);

	gStyle->SetOptFit(1111);
	h1_b1b2_2->Fit("landau", "", "", 100.0, 800.0);

	h1_b1b2_2->SetFillColor(kGreen);
	h1_b1b2_2->Draw();


	MyC_amplit_->cd(3);


	h1_b1b2_1_temp2->Smooth(Smooth_index);	 // Сглаживание


	n_number_temp = ntpl_select->Draw(var3_temp, cut_full_temp);

	gStyle->SetOptFit(1111);

	func_1->SetParameters(10.0, 2.0, -2.0, 2.0, -0.5);


	h1_b1b2_1_temp2->Fit("my_fit_1", "", "", -2, 2);
	h1_b1b2_1_temp2->Fit("my_fit_1", "", "", -2, 2);
	h1_b1b2_1_temp2->Draw();

	S1 = func_1->GetParameter(1);
	X1 = func_1->GetParameter(2);
	S2 = func_1->GetParameter(3);
	X2 = func_1->GetParameter(4);

	Double_t lambda = L_scint / delta_L(S1, X1, S2, X2);
	Double_t k_lambda = exp(L_scint_without / lambda);

	printf("lambda_scint= %.3f cm\n", lambda);
	printf("oslablenye= %.3f raz\n", k_lambda);
	printf("events nymber %i raz\n", n_number_temp);

	lambda_arr[i] = lambda;

	k_lambda_arr[i] = k_lambda;

	MyC_amplit_->Update();

}



// временное разрешение с коррекцией
void MyC_1() {

	name_cut();

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







	MyC_1->cd(2);



	MyC_1->cd(3);




	MyC_1->cd(4);



	TString vart1_1 = Form("0.5*((tF[0]-delta_t(aF[0], %f, %f, %f))+(tF[1]-delta_t(aF[1], %f, %f, %f)))>>h1_t1t2_1", par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right);


	ntpl_select->Draw(vart1_1, cut_full_1);

	h1_t1t2_1->Smooth(Smooth_index);	 // Сглаживание

	h1_t1t2_1->Fit("gaus", "", "", -200, 200);

	h1_t1t2_1->GetXaxis()->SetTitle("Time, channel number");
	h1_t1t2_1->GetYaxis()->SetTitle("Events number");
	h1_t1t2_1->GetXaxis()->CenterTitle();
	h1_t1t2_1->GetYaxis()->CenterTitle();


	h1_t1t2_1->Draw();

	MyC_1->cd(5);



	TString vart2_1 = Form("0.5*((tF[2]-delta_t(aF[2], %f, %f, %f))+(tF[3]-delta_t(aF[3], %f, %f, %f)))>>h2_t1t2_1", par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);


	ntpl_select->Draw(vart2_1, cut_full_2);

	h2_t1t2_1->Smooth(Smooth_index);	 // Сглаживание

	h2_t1t2_1->Fit("gaus", "", "", -200, 200);

	h2_t1t2_1->GetXaxis()->SetTitle("Time, channel number");
	h2_t1t2_1->GetYaxis()->SetTitle("Events number");
	h2_t1t2_1->GetXaxis()->CenterTitle();
	h2_t1t2_1->GetYaxis()->CenterTitle();


	h2_t1t2_1->Draw();

	MyC_1->cd(6);



	TString vart3_1 = Form("0.5*((tF[4]-delta_t(aF[4], %f, %f, %f))+(tF[5]-delta_t(aF[5], %f, %f, %f)))>>h3_t1t2_1", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right);


	ntpl_select->Draw(vart3_1, cut_full_3);

	h3_t1t2_1->Smooth(Smooth_index);	 // Сглаживание

	h3_t1t2_1->Fit("gaus", "", "", -200, 200);

	h3_t1t2_1->GetXaxis()->SetTitle("Time, channel number");
	h3_t1t2_1->GetYaxis()->SetTitle("Events number");
	h3_t1t2_1->GetXaxis()->CenterTitle();
	h3_t1t2_1->GetYaxis()->CenterTitle();


	h3_t1t2_1->Draw();

	MyC_1->cd(7);


	TString vart1_2 = Form("0.5*((tF[0]-delta_t(aF[0], %f, %f, %f))-(tF[1]-delta_t(aF[1], %f, %f, %f)))>>h1_t1t2_2", par0_tb1_1_right, par1_tb1_1_right, par2_tb1_1_right, par0_tb1_2_right, par1_tb1_2_right, par2_tb1_2_right);

	h1_t1t2_2->Smooth(Smooth_index);	 // Сглаживание

	ntpl_select->Draw(vart1_2, cut_full_1);



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



	TString vart2_2 = Form("0.5*((tF[2]-delta_t(aF[2], %f, %f, %f))-(tF[3]-delta_t(aF[3], %f, %f, %f)))>>h2_t1t2_2", par0_tb2_1_right, par1_tb2_1_right, par2_tb2_1_right, par0_tb2_2_right, par1_tb2_2_right, par2_tb2_2_right);

	h2_t1t2_2->Smooth(Smooth_index);	 // Сглаживание

	ntpl_select->Draw(vart2_2, cut_full_2);



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



	TString vart3_2 = Form("0.5*((tF[4]-delta_t(aF[4], %f, %f, %f))-(tF[5]-delta_t(aF[5],  %f, %f, %f)))>>h3_t1t2_2", par0_tb3_1_right, par1_tb3_1_right, par2_tb3_1_right, par0_tb3_2_right, par1_tb3_2_right, par2_tb3_2_right);

	h3_t1t2_2->Smooth(Smooth_index);	 // Сглаживание

	ntpl_select->Draw(vart3_2, cut_full_3);

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

	name_cut();

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


	TString var1 = Form("(aF[0] - %f) : (aF[1] - %f) >> h4_1", b1_1_pedestal, b1_2_pedestal);
	TString var1_same = Form("(aF[0] - %f) : (aF[1] - %f) >> h4_1_same", b1_1_pedestal, b1_2_pedestal);

	ntpl->Draw(var1, "", "", N_number);
	ntpl_select->Draw(var1_same, cut_full_1);

	h4_1->GetXaxis()->SetTitle("Amplitude, channel number");
	h4_1->GetYaxis()->SetTitle("Amplitude, channel number");
	h4_1->GetXaxis()->CenterTitle();
	h4_1->GetYaxis()->CenterTitle();


	h4_1->Draw();
	h4_1_same->SetMarkerColor(kRed);
	h4_1_same->Draw("same");

	MyC_2->cd(2);

	TString var2 = Form("(aF[2] - %f) : (aF[3] - %f) >> h5_1", b2_1_pedestal, b2_2_pedestal);
	TString var2_same = Form("(aF[2] - %f) : (aF[3] - %f) >> h5_1_same", b2_1_pedestal, b2_2_pedestal);



	ntpl->Draw(var2, "", "", N_number);
	ntpl_select->Draw(var2_same, cut_full_2);

	h5_1->GetXaxis()->SetTitle("Amplitude, channel number");
	h5_1->GetYaxis()->SetTitle("Amplitude, channel number");
	h5_1->GetXaxis()->CenterTitle();
	h5_1->GetYaxis()->CenterTitle();

	h5_1->Draw();
	h5_1_same->SetMarkerColor(kRed);
	h5_1_same->Draw("same");


	MyC_2->cd(3);

	TString var3 = Form("(aF[4] - %f) : (aF[5] - %f) >> h6_1", b3_1_pedestal, b3_2_pedestal);
	TString var3_same = Form("(aF[4] - %f) : (aF[5] - %f) >> h6_1_same", b3_1_pedestal, b3_2_pedestal);


	ntpl->Draw(var3, "", "", N_number);
	ntpl_select->Draw(var3_same, cut_full_3);


	h6_1->GetXaxis()->SetTitle("Amplitude, channel number");
	h6_1->GetYaxis()->SetTitle("Amplitude, channel number");
	h6_1->GetXaxis()->CenterTitle();
	h6_1->GetYaxis()->CenterTitle();

	h6_1->Draw();
	h6_1_same->SetMarkerColor(kRed);
	h6_1_same->SetMarkerStyle(6);
	h6_1_same->Draw("same");

	MyC_2->cd(4);

	TString var4 = Form("sqrt((aF[0] - %f)*(aF[1] - %f)) >> h1_b1b2_2", b1_1_pedestal, b1_2_pedestal);

	h1_b1b2_2->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(var4, cut_full_1);

	gStyle->SetOptFit(1111);
	h1_b1b2_2->Fit("landau", "", "", 100.0, 800.0);

	h1_b1b2_2->SetFillColor(kGreen);
	h1_b1b2_2->Draw();

	MyC_2->cd(5);

	TString var5 = Form("sqrt((aF[2] - %f)*(aF[3] - %f)) >> h2_b1b2_2", b2_1_pedestal, b2_2_pedestal);

	h2_b1b2_2->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(var5, cut_full_2);
	gStyle->SetOptFit(1111);
	h2_b1b2_2->Fit("landau", "", "", 100.0, 800.0);

	h2_b1b2_2->SetFillColor(kGreen);
	h2_b1b2_2->Draw();

	MyC_2->cd(6);

	TString var6 = Form("sqrt((aF[4] - %f)*(aF[5] - %f)) >> h3_b1b2_2", b3_1_pedestal, b3_2_pedestal);

	h3_b1b2_2->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(var6, cut_full_3);
	gStyle->SetOptFit(1111);
	h3_b1b2_2->Fit("landau", "", "", 100.0, 800.0);

	h3_b1b2_2->SetFillColor(kGreen);
	h3_b1b2_2->Draw();

	MyC_2->cd(7);

	TString var7 = Form("0.5*log((aF[0] - %f)/(aF[1] - %f)) >> h1_b1b2_1", b1_1_pedestal, b1_2_pedestal);

	h1_b1b2_1->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(var7, cut_full_1);

	gStyle->SetOptFit(1111);
	h1_b1b2_1->Fit("my_fit_1", "", "", -2, 2);
	h1_b1b2_1->Draw();

	S1 = func_1->GetParameter(1);
	X1 = func_1->GetParameter(2);
	S2 = func_1->GetParameter(3);
	X2 = func_1->GetParameter(4);

	printf("lambda_scint= %.3f cm\n", L_scint / delta_L(S1, X1, S2, X2));


	MyC_2->cd(8);

	TString var8 = Form("0.5*log((aF[2] - %f)/(aF[3] - %f)) >> h2_b1b2_1", b2_1_pedestal, b2_2_pedestal);

	h2_b1b2_1->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(var8, cut_full_2);

	gStyle->SetOptFit(1111);
	h2_b1b2_1->Fit("my_fit_1", "", "", -2, 2);
	h2_b1b2_1->Draw();

	S1 = func_1->GetParameter(1);
	X1 = func_1->GetParameter(2);
	S2 = func_1->GetParameter(3);
	X2 = func_1->GetParameter(4);

	printf("lambda_scint= %.3f cm\n", L_scint / delta_L(S1, X1, S2, X2));

	MyC_2->cd(9);

	TString var9 = Form("0.5*log((aF[4] - %f)/(aF[5] - %f)) >> h3_b1b2_1", b3_1_pedestal, b3_2_pedestal);

	h3_b1b2_1->Smooth(Smooth_index);	 // Сглаживание

	ntpl_select->Draw(var9, cut_full_3);

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

	name_cut();

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

	name_cut();

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

	name_cut();

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

	ntpl->Draw("aF[0]:tF[0]>>h1_1", "", "", N_number);
	ntpl_select->Draw("aF[0]:tF[0]>>h1_1_same", cut_full_1);

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

	ntpl->Draw("aF[1]:tF[1]>>h1_2", "", "", N_number);
	ntpl_select->Draw("aF[1]:tF[1]>>h1_2_same", cut_full_1);

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

	ntpl->Draw("aF[2]:tF[2]>>h2_1", "", "", N_number);
	ntpl_select->Draw("aF[2]:tF[2]>>h2_1_same", cut_full_2);

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

	ntpl->Draw("aF[3]:tF[3]>>h2_2", "", "", N_number);
	ntpl_select->Draw("aF[3]:tF[3]>>h2_2_same", cut_full_2);

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

	ntpl->Draw("aF[4]:tF[4]>>h3_1", "", "", N_number);
	ntpl_select->Draw("aF[4]:tF[4]>>h3_1_same", cut_full_3);

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


	ntpl->Draw("aF[5]:tF[5]>>h3_2", "", "", N_number);
	ntpl_select->Draw("aF[5]:tF[5]>>h3_2_same", cut_full_3);

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

	name_cut();

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


	name_cut();

	TH1F* h1_t1t2_1_no_corr = new TH1F("h1_t1t2_1_no_corr", "0.5*(t1_1+t1_2) WITHOUT AMPLITUDE CORRECTION", 1000, 500.0, 1500.0);
	TH1F* h2_t1t2_1_no_corr = new TH1F("h2_t1t2_1_no_corr", "0.5*(t2_1+t2_2) WITHOUT AMPLITUDE CORRECTION", 1000, 500.0, 1500.0);
	TH1F* h3_t1t2_1_no_corr = new TH1F("h3_t1t2_1_no_corr", "0.5*(t3_1+t3_2) WITHOUT AMPLITUDE CORRECTION", 1000, 500.0, 1500.0);

	TH1F* h1_t1t2_2_no_corr = new TH1F("h1_t1t2_2_no_corr", "0.5*(t1_1-t1_2) WITHOUT AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	TH1F* h2_t1t2_2_no_corr = new TH1F("h2_t1t2_2_no_corr", "0.5*(t2_1-t2_2) WITHOUT AMPLITUDE CORRECTION", 800, -400.0, 400.0);
	TH1F* h3_t1t2_2_no_corr = new TH1F("h3_t1t2_2_no_corr", "0.5*(t3_1-t3_2) WITHOUT AMPLITUDE CORRECTION", 800, -400.0, 400.0);


	TCanvas* MyC_7 = new TCanvas("MyC_7", "MyC_7", 0, 0, 2000, 1000);



	MyC_7->Divide(3, 3);


	MyC_7->cd(1);

	MyC_7->cd(2);

	MyC_7->cd(3);

	MyC_7->cd(4);

	TString vart1_1_no_corr = Form("0.5*(tF[0]+tF[1])>>h1_t1t2_1_no_corr");

	h1_t1t2_1_no_corr->Smooth(Smooth_index);	 // Сглаживание

	ntpl_select->Draw(vart1_1_no_corr, cut_full_1);

	h1_t1t2_1_no_corr->Fit("gaus", "", "", 700, 1300);

	h1_t1t2_1_no_corr->GetXaxis()->SetTitle("Time, channel number");
	h1_t1t2_1_no_corr->GetYaxis()->SetTitle("Events number");
	h1_t1t2_1_no_corr->GetXaxis()->CenterTitle();
	h1_t1t2_1_no_corr->GetYaxis()->CenterTitle();

	h1_t1t2_1_no_corr->Draw();

	MyC_7->cd(5);



	TString vart2_1_no_corr = Form("0.5*(tF[2]+tF[3])>>h2_t1t2_1_no_corr");

	h2_t1t2_1_no_corr->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(vart2_1_no_corr, cut_full_2);

	h2_t1t2_1_no_corr->Fit("gaus", "", "", 700, 1300);

	h2_t1t2_1_no_corr->GetXaxis()->SetTitle("Time, channel number");
	h2_t1t2_1_no_corr->GetYaxis()->SetTitle("Events number");
	h2_t1t2_1_no_corr->GetXaxis()->CenterTitle();
	h2_t1t2_1_no_corr->GetYaxis()->CenterTitle();

	h2_t1t2_1_no_corr->Draw();

	MyC_7->cd(6);



	TString vart3_1_no_corr = Form("0.5*(tF[4]+tF[5])>>h3_t1t2_1_no_corr");

	h3_t1t2_1_no_corr->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(vart3_1_no_corr, cut_full_3);

	h3_t1t2_1_no_corr->Fit("gaus", "", "", 700, 1300);

	h3_t1t2_1_no_corr->GetXaxis()->SetTitle("Time, channel number");
	h3_t1t2_1_no_corr->GetYaxis()->SetTitle("Events number");
	h3_t1t2_1_no_corr->GetXaxis()->CenterTitle();
	h3_t1t2_1_no_corr->GetYaxis()->CenterTitle();

	h3_t1t2_1_no_corr->Draw();

	MyC_7->cd(7);


	TString vart1_2_no_corr = Form("0.5*(tF[0]-tF[1])>>h1_t1t2_2_no_corr");

	h1_t1t2_2_no_corr->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(vart1_2_no_corr, cut_full_1);


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



	TString vart2_2_no_corr = Form("0.5*(tF[2]-tF[3])>>h2_t1t2_2_no_corr");

	h2_t1t2_2_no_corr->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(vart2_2_no_corr, cut_full_2);



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



	TString vart3_2_no_corr = Form("0.5*(tF[4]-tF[5])>>h3_t1t2_2_no_corr");

	h3_t1t2_2_no_corr->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(vart3_2_no_corr, cut_full_3);

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






// для определения порого дискриминатора (в каналах амплитуды) // по каждому каналу амплитуды
void MyC_9() {

	name_cut();

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
	double y1 = hb1_1_divide->GetBinContent(i - 1);
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

	name_cut();

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


	TString var1 = Form("(aF[0] - %f) / (aF[1] - %f) >> h1_b1b2_phe_1", b1_1_pedestal, b1_2_pedestal);

	h1_b1b2_phe_1->Smooth(Smooth_index);	 // Сглаживание

	ntpl_select->Draw(var1, cut_full_1);

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

	TString var2 = Form("(aF[2] - %f) / (aF[3] - %f) >> h2_b1b2_phe_1", b2_1_pedestal, b2_2_pedestal);

	h2_b1b2_phe_1->Smooth(Smooth_index);	 // Сглаживание

	ntpl_select->Draw(var2, cut_full_2);

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

	TString var3 = Form("(aF[4] - %f) / (aF[5] - %f) >> h3_b1b2_phe_1", b3_1_pedestal, b3_2_pedestal);

	h3_b1b2_phe_1->Smooth(Smooth_index);	 // Сглаживание

	ntpl_select->Draw(var3, cut_full_3);

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

	TString var4 = Form("((aF[0] - %f) - (aF[1] - %f)) / ((aF[0] - %f) + (aF[1] - %f)) >> h1_b1b2_phe_2", b1_1_pedestal, b1_2_pedestal, b1_1_pedestal, b1_2_pedestal);

	h1_b1b2_phe_2->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(var4, cut_full_1);

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

	TString var5 = Form("((aF[2] - %f) - (aF[3] - %f)) / ((aF[2] - %f) + (aF[3] - %f)) >> h2_b1b2_phe_2", b2_1_pedestal, b2_2_pedestal, b2_1_pedestal, b2_2_pedestal);

	h2_b1b2_phe_2->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(var5, cut_full_2);

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

	TString var6 = Form("((aF[4] - %f) - (aF[5] - %f)) / ((aF[4] - %f) + (aF[5] - %f)) >> h3_b1b2_phe_2", b3_1_pedestal, b3_2_pedestal, b3_1_pedestal, b3_2_pedestal);

	h3_b1b2_phe_2->Smooth(Smooth_index);	 // Сглаживание


	ntpl_select->Draw(var6, cut_full_3);

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


void a() {
	
		TH1F* h[17];
		int i;
		h[0] = new TH1F("h0", "h0", 100, 0., 8000.);
		h[0]->SetLineWidth(2);
		h[0]->SetLineColor(kBlue + 1);
		for (i = 1; i < 17; i++) {
			h[i] = (TH1F*)h[0]->Clone(Form("h%d", i));
			h[i]->SetTitle(Form("h%02d", i));
		}

		TCanvas* cc = new TCanvas("cc", "cc", 1000, 1000);
		gStyle->SetOptStat(10);
		cc->Divide(4, 4, 0.01, 0.01);

		for (i = 1; i < 17; i++) {
			cc->cd(i);
			ntpl->Draw(Form("16384-aZ[%d]>>h%d", i, i), Form("regZ&%d", 1 << (i - 1)));
			cc->Update();
			cc->Flush();
		}

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
		MyC_9();
		MyC_10();


		});


	printf("Running something in the \"main\" thread\n");

	// Wait until all items are complete
	tg.Wait();

	printf("All work completed.\n");








	return 0;
}
