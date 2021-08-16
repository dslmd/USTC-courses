#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "math.h"
#include "string.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3D.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TTree.h"
#include "TSystem.h"
#include "TRandom.h"


int MRPCcut1(Float_t mrpc11_t[32],Float_t mrpc12_t[32],Float_t mrpc21_t[32],Float_t mrpc22_t[32]){
	if(((mrpc11_t[0]>0&&mrpc12_t[31]>0)||(mrpc11_t[1]>0&&mrpc12_t[30]>0)||(mrpc11_t[2]>0&&mrpc12_t[29]>0)||(mrpc11_t[3]>0&&mrpc12_t[28]>0)||(mrpc11_t[4]>0&&mrpc12_t[27]>0)||(mrpc11_t[5]>0&&mrpc12_t[26]>0)||
	   (mrpc11_t[6]>0&&mrpc12_t[25]>0)||(mrpc11_t[7]>0&&mrpc12_t[24]>0)||(mrpc11_t[8]>0&&mrpc12_t[23]>0)||(mrpc11_t[9]>0&&mrpc12_t[22]>0)||(mrpc11_t[10]>0&&mrpc12_t[21]>0)||(mrpc11_t[11]>0&&mrpc12_t[20]>0)||
    ( mrpc11_t[12]>0&&mrpc12_t[19]>0)||(mrpc11_t[13]>0&&mrpc12_t[18]>0)||(mrpc11_t[14]>0&&mrpc12_t[17]>0)||(mrpc11_t[15]>0&&mrpc12_t[16]>0)||(mrpc11_t[16]>0&&mrpc12_t[15]>0)||(mrpc11_t[17]>0&&mrpc12_t[14]>0)||(mrpc11_t[18]>0&&mrpc12_t[13]>0)
     ||(mrpc11_t[19]>0&&mrpc12_t[12]>0)||(mrpc11_t[20]>0&&mrpc12_t[11]>0)||(mrpc11_t[21]>0&&mrpc12_t[10]>0)||(mrpc11_t[22]>0&&mrpc12_t[9]>0)||(mrpc11_t[23]>0&&mrpc12_t[8]>0)||(mrpc11_t[24]>0&&mrpc12_t[7]>0)||(mrpc11_t[25]>0&&mrpc12_t[6]>0)||(mrpc11_t[26]>0&&mrpc12_t[5]>0)
     ||(mrpc11_t[27]>0&&mrpc12_t[4]>0)||(mrpc11_t[28]>0&&mrpc12_t[3]>0)||(mrpc11_t[29]>0&&mrpc12_t[2]>0)||(mrpc11_t[30]>0&&mrpc12_t[1]>0)||(mrpc11_t[31]>0&&mrpc12_t[0]>0))&&

  ((mrpc21_t[0]>0&&mrpc22_t[31]>0)||(mrpc21_t[1]>0&&mrpc22_t[30]>0)||(mrpc21_t[2]>0&&mrpc22_t[29]>0)||(mrpc21_t[3]>0&&mrpc22_t[28]>0)||(mrpc21_t[4]>0&&mrpc22_t[27]>0)||(mrpc21_t[5]>0&&mrpc22_t[26]>0)||
    (mrpc21_t[6]>0&&mrpc22_t[25]>0)||(mrpc21_t[7]>0&&mrpc22_t[24]>0)||(mrpc21_t[8]>0&&mrpc22_t[23]>0)||(mrpc21_t[9]>0&&mrpc22_t[22]>0)||(mrpc21_t[10]>0&&mrpc22_t[21]>0)||(mrpc21_t[11]>0&&mrpc22_t[20]>0)||
   ( mrpc21_t[12]>0&&mrpc22_t[19]>0)||(mrpc21_t[13]>0&&mrpc22_t[18]>0)||(mrpc21_t[14]>0&&mrpc22_t[17]>0)||(mrpc21_t[15]>0&&mrpc22_t[16]>0)||(mrpc21_t[16]>0&&mrpc22_t[15]>0)||(mrpc21_t[17]>0&&mrpc22_t[14]>0)||(mrpc21_t[18]>0&&mrpc22_t[13]>0)
     ||(mrpc21_t[19]>0&&mrpc22_t[12]>0)||(mrpc21_t[20]>0&&mrpc22_t[11]>0)||(mrpc21_t[21]>0&&mrpc22_t[10]>0)||(mrpc21_t[22]>0&&mrpc22_t[9]>0)||(mrpc21_t[23]>0&&mrpc22_t[8]>0)||(mrpc21_t[24]>0&&mrpc22_t[7]>0)||(mrpc21_t[25]>0&&mrpc22_t[6]>0)||(mrpc21_t[26]>0&&mrpc22_t[5]>0)
     ||(mrpc21_t[27]>0&&mrpc22_t[4]>0)||(mrpc21_t[28]>0&&mrpc22_t[3]>0)||(mrpc21_t[29]>0&&mrpc22_t[2]>0)||(mrpc21_t[30]>0&&mrpc22_t[1]>0)||(mrpc21_t[31]>0&&mrpc22_t[0]>0)))return 0;

	else return 1;
}

int MRPCcut2(Float_t mrpc11_TOT[32],Float_t mrpc12_TOT[32],Float_t mrpc21_TOT[32],Float_t mrpc22_TOT[32])
{
  if(((mrpc11_TOT[0]>0&&mrpc12_TOT[31]>0)||(mrpc11_TOT[1]>0&&mrpc12_TOT[30]>0)||(mrpc11_TOT[2]>0&&mrpc12_TOT[29]>0)||(mrpc11_TOT[3]>0&&mrpc12_TOT[28]>0)||(mrpc11_TOT[4]>0&&mrpc12_TOT[27]>0)||(mrpc11_TOT[5]>0&&mrpc12_TOT[26]>0)||
  (mrpc11_TOT[6]>0&&mrpc12_TOT[25]>0)||(mrpc11_TOT[7]>0&&mrpc12_TOT[24]>0)||(mrpc11_TOT[8]>0&&mrpc12_TOT[23]>0)||(mrpc11_TOT[9]>0&&mrpc12_TOT[22]>0)||(mrpc11_TOT[10]>0&&mrpc12_TOT[21]>0)||(mrpc11_TOT[11]>0&&mrpc12_TOT[20]>0)||
 ( mrpc11_TOT[12]>0&&mrpc12_TOT[19]>0)||(mrpc11_TOT[13]>0&&mrpc12_TOT[18]>0)||(mrpc11_TOT[14]>0&&mrpc12_TOT[17]>0)||(mrpc11_TOT[15]>0&&mrpc12_TOT[16]>0))&&

((mrpc21_TOT[0]>0&&mrpc22_TOT[31]>0)||(mrpc21_TOT[1]>0&&mrpc22_TOT[30]>0)||(mrpc21_TOT[2]>0&&mrpc22_TOT[29]>0)||(mrpc21_TOT[3]>0&&mrpc22_TOT[28]>0)||(mrpc21_TOT[4]>0&&mrpc22_TOT[27]>0)||(mrpc21_TOT[5]>0&&mrpc22_TOT[26]>0)||
 (mrpc21_TOT[6]>0&&mrpc22_TOT[25]>0)||(mrpc21_TOT[7]>0&&mrpc22_TOT[24]>0)||(mrpc21_TOT[8]>0&&mrpc22_TOT[23]>0)||(mrpc21_TOT[9]>0&&mrpc22_TOT[22]>0)||(mrpc21_TOT[10]>0&&mrpc22_TOT[21]>0)||(mrpc21_TOT[11]>0&&mrpc22_TOT[20]>0)||
( mrpc21_TOT[12]>0&&mrpc22_TOT[19]>0)||(mrpc21_TOT[13]>0&&mrpc22_TOT[18]>0)||(mrpc21_TOT[14]>0&&mrpc22_TOT[17]>0)||(mrpc21_TOT[15]>0&&mrpc22_TOT[16]>0)))return 0;

else return 1;
}
int T0cut(Float_t mrpc11_TOT[32],Float_t mrpc12_TOT[32],Float_t mrpc21_TOT[32],Float_t mrpc22_TOT[32],Int_t padnum1,Int_t padnum2){
	if(!(mrpc11_TOT[padnum1]>=1500&&mrpc11_TOT[padnum1]<=6000)) return 1;
	if(!(mrpc12_TOT[31-padnum1]>=1500&&mrpc12_TOT[31-padnum1]<=6000)) return 1;
	if(!(mrpc21_TOT[padnum2]>=1500&&mrpc21_TOT[padnum2]<=6000)) return 1;
	if(!(mrpc22_TOT[63-padnum2]>=1500&&mrpc22_TOT[63-padnum2]<=6000)) return 1;


	else return 0;
}

int T0Res(){
//define variables*************************************


Float_t mrpc11_T[32],mrpc11_L[32],mrpc11_TOT[32],mrpc12_T[32],mrpc12_L[32],mrpc12_TOT[32];

Float_t mrpc21_T[32],mrpc21_L[32],mrpc21_TOT[32],mrpc22_T[32],mrpc22_L[32],mrpc22_TOT[32];


Float_t MQ11,MQ12,MQ21,MQ22;
const Int_t padnum1=16;
const Int_t padnum2=16;//padnum1,padnum2<32


const char *treefile="CBM_test_6000V";
char buf[1024];
Double_t detector1_eff,detector2_eff,detector3_eff;
Double_t all_event=0.;
Double_t real_event1=0.,real_event2=0.,real_event3=0.;
Int_t fireN1,fireN2;

//read data from TTree*********************
sprintf(buf,"%s.root",treefile);
TFile*mfile=new TFile(buf,"UPDATA");
TTree*mtree=(TTree*)mfile->Get("mtree");
	mtree->SetBranchAddress("mrpc11_T",&mrpc11_T[0]);
  mtree->SetBranchAddress("mrpc11_L",&mrpc11_L[0]);
  mtree->SetBranchAddress("mrpc11_TOT",&mrpc11_TOT[0]);

	mtree->SetBranchAddress("mrpc12_T",&mrpc12_T[0]);
  mtree->SetBranchAddress("mrpc12_L",&mrpc12_L[0]);
  mtree->SetBranchAddress("mrpc12_TOT",&mrpc12_TOT[0]);

  mtree->SetBranchAddress("mrpc21_T",&mrpc21_T[0]);
  mtree->SetBranchAddress("mrpc21_L",&mrpc21_L[0]);
  mtree->SetBranchAddress("mrpc21_TOT",&mrpc21_TOT[0]);

  mtree->SetBranchAddress("mrpc22_T",&mrpc22_T[0]);
  mtree->SetBranchAddress("mrpc22_L",&mrpc22_L[0]);
  mtree->SetBranchAddress("mrpc22_TOT",&mrpc22_TOT[0]);
Int_t entries=mtree->GetEntries();
cout<<"The total events are "<<entries<<endl;
//******************create canvas and histogram*******************
TCanvas*c1= new TCanvas("c1","c1",800,800);
  c1->SetFillColor(10);
  c1->SetFrameFillColor(10);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderSize(0);
  c1->SetLeftMargin(0.145);
  c1->SetRightMargin(0.065);

  TCanvas*c2= new TCanvas("c2","c2",800,800);
  c2->Divide(2,3);
  c2->SetFillColor(10);
  c2->SetFrameFillColor(10);
  c2->SetFrameBorderMode(0);
  c2->SetFrameBorderSize(0);
  c2->SetLeftMargin(0.145);
  c2->SetRightMargin(0.065);

  TCanvas*c3= new TCanvas("c3","c3",800,1600);
  c3->Divide(1,2);
  c3->SetFillColor(10);
  c3->SetFrameFillColor(10);
  c3->SetFrameBorderMode(0);
  c3->SetFrameBorderSize(0);
  c3->SetLeftMargin(0.145);
  c3->SetRightMargin(0.065);

  TCanvas*c4= new TCanvas("c4","c4",800,800);
  c4->Divide(1,2);
  c4->SetFillColor(10);
  c4->SetFrameFillColor(10);
  c4->SetFrameBorderMode(0);
  c4->SetFrameBorderSize(0);
  c4->SetLeftMargin(0.145);
  c4->SetRightMargin(0.065);


  TCanvas*c5= new TCanvas("c5","c5",800,800);
  c5->Divide(1,2);
  c5->SetFillColor(10);
  c5->SetFrameFillColor(10);
  c5->SetFrameBorderMode(0);
  c5->SetFrameBorderSize(0);
  c5->SetLeftMargin(0.145);
  c5->SetRightMargin(0.065);

  TCanvas*c6= new TCanvas("c6","c6",800,800);
  c6->SetFillColor(10);
  c6->SetFrameFillColor(10);
  c6->SetFrameBorderMode(0);
  c6->SetFrameBorderSize(0);
  c6->SetLeftMargin(0.145);
  c6->SetRightMargin(0.065);

  TCanvas*c7= new TCanvas("c7","c7",800,800);
  c7->Divide(1,2);
  c7->SetFillColor(10);
  c7->SetFrameFillColor(10);
  c7->SetFrameBorderMode(0);
  c7->SetFrameBorderSize(0);
  c7->SetLeftMargin(0.145);
  c7->SetRightMargin(0.065);





  gStyle->SetOptStat(11);
  gStyle->SetOptFit(111);

   TH1F *T0= new TH1F("T0","T0",100,-5000,5000);
     T0->GetXaxis()->SetTitle("(t1+t2-t3-t4)/2");
     T0->GetXaxis()->SetTitleOffset(1.6);
     T0->GetYaxis()->SetTitle("count");
     T0->GetYaxis()->SetTitleOffset(1.6);

  TH1F *Q1 =  new TH1F("Q1","Q1",100,0,8000);
    Q1->GetXaxis()->SetTitle("qdc");
    Q1->GetYaxis()->SetTitle("count");
    Q1->GetYaxis()->SetTitleOffset(1.6);
  TH1F *Q2 =  new TH1F("Q2","Q2",100,0,8000);
    Q2->GetXaxis()->SetTitle("qdc");
    Q2->GetYaxis()->SetTitle("count");
    Q2->GetYaxis()->SetTitleOffset(1.6);
  TH1F *Q3 =  new TH1F("Q3","Q3",100,0,8000);
    Q3->GetXaxis()->SetTitle("qdc");
    Q3->GetYaxis()->SetTitle("count");
    Q3->GetYaxis()->SetTitleOffset(1.6);
  TH1F *Q4 =  new TH1F("Q4","Q4",100,0,8000);
    Q4->GetXaxis()->SetTitle("qdc");
    Q4->GetYaxis()->SetTitle("count");
    Q4->GetYaxis()->SetTitleOffset(1.6);
  TH1F *Q1Q2 =  new TH1F("Q1Q2","Q1Q2",200,0,16000);
    Q1Q2->GetXaxis()->SetTitle("qdc");
    Q1Q2->GetYaxis()->SetTitle("count");
    Q1Q2->GetYaxis()->SetTitleOffset(1.6);
  TH1F *Q3Q4 =  new TH1F("Q3Q4","Q3Q4",200,0,16000);
    Q3Q4->GetXaxis()->SetTitle("qdc");
    Q3Q4->GetYaxis()->SetTitle("count");
    Q3Q4->GetYaxis()->SetTitleOffset(1.6);

 TH2D*Q1_Q2=new TH2D("Q1_Q2","Q1_Q2",200,0,12000,200,0,12000);
	 Q1_Q2->GetXaxis()->SetTitle("Q1");
     Q1_Q2->GetYaxis()->SetTitle("Q2");
 TH2D*Q3_Q4=new TH2D("Q3_Q4","Q3_Q4",200,0,12000,200,0,12000);
	 Q3_Q4->GetXaxis()->SetTitle("Q3");
     Q3_Q4->GetYaxis()->SetTitle("Q4");
 TH2D *T1vsTOT=new TH2D("T1vsTOT","T1vsTOT",200,0,12000,200,-3000,3000);
    T1vsTOT->GetXaxis()->SetTitle("(q1+q2)/2");
	T1vsTOT->GetYaxis()->SetTitle("(t1+t2-t3-t4)/2");
 TH2D *T2vsTOT=new TH2D("T2vsTOT","T2vsTOT",200,0,12000,200,-3000,3000);
    T2vsTOT->GetXaxis()->SetTitle("(q3+q4)/2");
  T2vsTOT->GetYaxis()->SetTitle("(t3+t4-t1-t2)/2");
  

  TH2D *T1vsTOT_correction=new TH2D("T1vsTOT_correction","T1vsTOT_correction",200,0,12000,200,-3000,3000);
  T1vsTOT_correction->GetXaxis()->SetTitle("(q1+q2)/2");
T1vsTOT_correction->GetYaxis()->SetTitle("((t1+t2-t3-t4)/2");
TH2D *T2vsTOT_correction=new TH2D("T2vsTOT_correction","T2vsTOT_correction",200,0,12000,200,-3000,3000);
  T2vsTOT_correction->GetXaxis()->SetTitle("(q3+q4)/2");
T2vsTOT_correction->GetYaxis()->SetTitle("((t3+t4-t1-t2)/2");

 TH1F*pT1vsTOT;
 TH1F*pT2vsTOT;
 TH2D *test=new TH2D("test","test",200,0,12000,200,-3000,3000);

 TF1 *mrpc1_fun;  //MRPC T-A correction
	  sprintf(buf,"mrpc1fun_%d",padnum1);
	  mrpc1_fun = new TF1(buf, "[0]+[1]/sqrt(x)+[2]/x+[3]/x/sqrt(x)+[4]/x/x", 1000, 8000);
	  mrpc1_fun->SetLineWidth(1);
	  mrpc1_fun->SetLineColor(4);
 TF1 *mrpc2_fun;  //MRPC T-A correction
	  sprintf(buf,"mrpc2fun_%d",padnum2);
	  mrpc2_fun = new TF1(buf, "[0]+[1]/sqrt(x)+[2]/x+[3]/x/sqrt(x)+[4]/x/x", 1000, 8000);
	  mrpc2_fun->SetLineWidth(1);
	  mrpc2_fun->SetLineColor(4);

  TF1 *g = new TF1("g", "gaus", -3000, 3000);
   g->SetLineColor(2);
   g->SetLineWidth(1);


//*********************caculate**********************************


for(int j=0;j<entries;j++){
  mtree->GetEntry(j);
  if(MRPCcut1(mrpc11_L,mrpc12_L,mrpc21_L,mrpc22_L)==1) continue;
  if(MRPCcut1(mrpc11_TOT,mrpc12_TOT,mrpc21_TOT,mrpc22_TOT)==1) continue;

fireN1 = 0;
fireN2 = 0;
  for(int i=0;i<32;i++)
  {
    if((mrpc11_L[i]>0)&&(mrpc12_L[31-i]>0))
    fireN1=fireN1+1;
    if((mrpc21_L[i]>0)&&(mrpc22_L[31-i]>0))
    fireN2++;
  }
if((fireN1<3)&&(fireN2<3)&&(mrpc11_TOT[padnum1]>0)&&(mrpc12_TOT[31-padnum1]>0)&&(mrpc21_TOT[padnum2]>0)&&(mrpc22_TOT[31-padnum2]>0))
 {
  Q1->Fill(mrpc11_TOT[padnum1]);
	Q2->Fill(mrpc12_TOT[31-padnum1]);
	Q3->Fill(mrpc21_TOT[padnum2]);
	Q4->Fill(mrpc22_TOT[31-padnum2]);

	Q1Q2->Fill(mrpc11_TOT[padnum1]+mrpc12_TOT[31-padnum1]);
	Q3Q4->Fill(mrpc21_TOT[padnum2]+mrpc22_TOT[31-padnum2]);

	Q1_Q2->Fill(mrpc11_TOT[padnum1],mrpc12_TOT[31-padnum1]);
	Q3_Q4->Fill(mrpc21_TOT[padnum2],mrpc22_TOT[31-padnum2]);
	T1vsTOT->Fill((mrpc11_TOT[padnum1]+mrpc12_TOT[31-padnum1])/2,(mrpc11_L[padnum1]+mrpc12_L[31-padnum1]-mrpc21_L[padnum2]-mrpc22_L[31-padnum2])/2);
	T2vsTOT->Fill((mrpc21_TOT[padnum2]+mrpc22_TOT[31-padnum2])/2,(mrpc21_L[padnum2]+mrpc22_L[31-padnum2]-mrpc11_L[padnum1]-mrpc12_L[31-padnum1])/2);


	T0->Fill((mrpc11_L[padnum1]+mrpc12_L[31-padnum1]-mrpc21_L[padnum2]-mrpc22_L[31-padnum2])/2);

}
}
c1->cd();
T0->Fit("g","","",-2000,2000);
//T0->Fit("g");
T0->Draw();
//c1->SaveAs("TORes.png");

c2->cd(1);
Q1->Draw();
c2->cd(2);
Q2->Draw();
c2->cd(3);
Q3->Draw();
c2->cd(4);
Q4->Draw();
c2->cd(5);
Q1Q2->Draw();
c2->cd(6);
Q3Q4->Draw();
//c2->SaveAs("qdc.png");

    c3->cd(1);
	Q1_Q2->Draw("colz");
	c3->cd(2);
	Q3_Q4->Draw("colz");
//	c3->SaveAs("Q_Q.png");

	c4->cd(1);
	T1vsTOT->Draw("colz");
	c4->cd(2);
	T2vsTOT->Draw("colz");
//	c4->SaveAs("TvsTOT.png");


	//--------------------------T-A correctioo--------------------------------------------
  T1vsTOT->ProfileX();
    pT1vsTOT = (TH1F*)gDirectory->Get("T1vsTOT_pfx");
   c5->cd(1);
    pT1vsTOT->SetName("T1vsTOT_pfx");
    pT1vsTOT->SetTitle("T1vsTOT_pfx");
    pT1vsTOT->GetXaxis()->SetTitle("TOT");
    pT1vsTOT->GetYaxis()->SetTitleOffset(1.6);
    pT1vsTOT->GetYaxis()->SetTitle("tdc");
	pT1vsTOT->GetYaxis()->SetRangeUser(-3000,3000);
	pT1vsTOT->Fit(mrpc1_fun,"R+","",1400,8000);
 pT1vsTOT->Draw();

 T2vsTOT->ProfileX();
    pT2vsTOT = (TH1F*)gDirectory->Get("T2vsTOT_pfx");
    c5->cd(2);
    pT2vsTOT->SetName("T2vsTOT_pfx");
    pT2vsTOT->SetTitle("T2vsTOT_pfx");
    pT2vsTOT->GetXaxis()->SetTitle("TOT");
    pT2vsTOT->GetYaxis()->SetTitleOffset(1.6);
    pT2vsTOT->GetYaxis()->SetTitle("tdc");
	pT2vsTOT->GetYaxis()->SetRangeUser(-3000,3000);
	pT2vsTOT->Fit(mrpc2_fun,"R+","",1000,8000);

  pT2vsTOT->Draw();
 
  
  
//  c5->SaveAs("TvsTOT_profile.png");

test->ProfileX();
//-------------------------------Resolution after mrpc T-A corection-------------
 
  Float_t mrpc1_cor;
  Float_t mrpc2_cor;

  TH1F *T02= new TH1F("T02","T02",100,-5000,5000);
     T02->GetXaxis()->SetTitle("(t1+t2-t3-t4)/2");
     T02->GetXaxis()->SetTitleOffset(1.6);
     T02->GetYaxis()->SetTitle("count");
     T02->GetYaxis()->SetTitleOffset(1.6);

  for(int j=0;j<entries;j++){
   mtree->GetEntry(j);

   if(MRPCcut1(mrpc11_L,mrpc12_L,mrpc21_L,mrpc22_L)==1) continue;
   if(MRPCcut1(mrpc11_TOT,mrpc12_TOT,mrpc21_TOT,mrpc22_TOT)==1) continue;

 fireN1 = 0;
 fireN2 = 0;
   for(int i=0;i<32;i++)
   {

     if((mrpc11_L[i]>0)&&(mrpc12_L[31-i]>0))
     fireN1=fireN1+1;
     if((mrpc21_L[i]>0)&&(mrpc22_L[31-i]>0))
     fireN2++;
   }
   if((fireN1<3)&&(fireN2<3)&&(mrpc11_TOT[padnum1]>0)&&(mrpc12_TOT[31-padnum1]>0)&&(mrpc21_TOT[padnum2]>0)&&(mrpc22_TOT[31-padnum2]>0))
{
  mrpc1_cor=mrpc1_fun->Eval((mrpc11_TOT[padnum1]+mrpc12_TOT[31-padnum1])/2);
  mrpc2_cor=mrpc2_fun->Eval((mrpc21_TOT[padnum2]+mrpc22_TOT[31-padnum2])/2);
  T1vsTOT_correction->Fill((mrpc11_TOT[padnum1]+mrpc12_TOT[31-padnum1])/2,(mrpc11_L[padnum1]+mrpc12_L[31-padnum1]-mrpc21_L[padnum2]-mrpc22_L[31-padnum2])/2-mrpc1_cor);
	T2vsTOT_correction->Fill((mrpc21_TOT[padnum2]+mrpc22_TOT[31-padnum2])/2,(mrpc21_L[padnum2]+mrpc22_L[31-padnum2]-mrpc11_L[padnum1]-mrpc12_L[31-padnum1])/2-mrpc2_cor);
  if(mrpc11_TOT[padnum1]<4200)
  T02->Fill(((mrpc11_L[padnum1]+mrpc12_L[31-padnum1])/2-(mrpc21_L[padnum2]+mrpc22_L[31-padnum2])/2-mrpc1_cor+mrpc2_cor));
 // else
  //{
  //  T02->Fill(((mrpc11_L[padnum1]+mrpc12_L[31-padnum1])/2-(mrpc21_L[padnum2]+mrpc22_L[31-padnum2])/2+mrpc2_cor-mrpc1_cor));
  //}
  
}
  }


  c6->cd();
  T02->Fit("g","","",-2000,2000);
  T02->Draw();
 // c6->SaveAs("TORes1.png");

  c7->cd(1);
	T1vsTOT_correction->Draw("colz");
	c7->cd(2);
	T2vsTOT_correction->Draw("colz");
//	c7->SaveAs("TvsTOT_correction.png");

  return 0;
}








