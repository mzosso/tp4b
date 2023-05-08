#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <typeinfo>
#include <TRandom.h>
#include <TVectorT.h>
using namespace std;

void create_plots()
{
    ifstream output("output_area.txt");
    ifstream outputT("output_T.txt");
    ifstream outputV("output_V.txt");
    //vector<ifstream> outputs={output, outputT, outputV};
    vector<string> var={"n_d", "T", "V"};
    string line;
    string line1;
    string line2;
    
    //TVectorD sigma_nd(2);
    //TVectorD nd(2);
    TVectorD areas(4);
    TVectorD sigma_area(4);
    TVectorD sigma_T(6);
    TVectorD Ts(6);
    
    TVectorD sigma_V(5);
    double arr[] = {110, 200, 430, 500, 600, 300};
    double arr_V[] = { 300, 1000, 1500, 2000, 600};
    TVectorD Vs(5, arr_V);
    TVectorD Tss(6, arr);
    int count=0;
    while (getline(output, line)) 
    {
        //cout<<"line="<<line<<endl;
        //stringstream ss(line);
       // string valueString;
        //ss >> valueString;
        //cout<<"valueString ="<<valueString<<endl;
       
        if(line.find("sigma")!=string::npos)
        {
            int pos = line.find("=");
            // Copy substring after pos
            string sub = line.substr(pos + 1);
            sigma_area[count]=stod(sub);
            cout<<"sigmas area="<<sigma_area[count]<<endl;
            count+=1;
        }
        if(line.find("area= ")!=string::npos)
        {
            int pos = line.find("=");
            // Copy substring after pos
            string sub = line.substr(pos + 1);
            areas[count]=stod(sub);
           //cout<<"areas ="<<areas[count]<<endl;
        }
        
    }
    count=0;
    while (getline(outputT, line1)) 
    {
        //cout<<"line="<<line<<endl;
        //stringstream ss(line);
       // string valueString;
        //ss >> valueString;
        //cout<<"valueString ="<<valueString<<endl;
       
        if(line1.find("sigma")!=string::npos)
        {
            int pos = line1.find("=");
            // Copy substring after pos
            string sub = line1.substr(pos + 1);
            sigma_T[count]=stod(sub);
           cout<<"sigmas T="<<sigma_T[count]<<endl;
           count+=1;
        }
        /*if(line.find("T= ")!=string::npos)
        {
            cout<<"in if T"<<endl;
            int pos = line1.find("=");
            // Copy substring after pos
            string sub = line1.substr(pos + 1);
            Ts[count]=stod(sub);
            //cout<<"T ="<<Ts[count]<<endl;
        }*/
        
    }
    count=0;
    while (getline(outputV, line2)) 
    {
        //cout<<"line="<<line<<endl;
        //stringstream ss(line);
       // string valueString;
        //ss >> valueString;
        //cout<<"valueString ="<<valueString<<endl;
        
        if(line2.find("sigma")!=string::npos)
        {
            int pos = line2.find("=");
            // Copy substring after pos
            string sub = line2.substr(pos + 1);
            sigma_V[count]=stod(sub);
            cout<<"sigmas V="<<sigma_V[count]<<endl;
            cout<<"count   "<<count<<endl;
            count+=1;
        }
        /*if(line2.find("V=")!=string::npos)
        {
            int pos = line2.find("=");
            // Copy substring after pos
            string sub = line2.substr(pos + 1);
            Vs[count]=stod(sub);
            //cout<<"V ="<<Vs[count]<<endl;
        }*/
        
    }
    



    TCanvas *canv_area = new TCanvas("canv_area", "area", 200, 10, 1000, 650);
	canv_area->SetGrid();
	
	TH2F *hpxnd = new TH2F("hpxnd", "area", 10000, 9*1e3, 3*1e6, 50, 0, 2);
	
	hpxnd->Draw();
	
	hpxnd->GetYaxis()->SetTitle("time resolution [ns]");
	
	
	hpxnd->GetYaxis()->SetLabelSize(0.05);
	hpxnd->GetXaxis()->SetLabelSize(0.05);
	hpxnd->GetYaxis()->SetTitleSize(0.05);
	hpxnd->GetXaxis()->SetTitleSize(0.05);
	hpxnd->GetXaxis()->SetNoExponent(0);
	hpxnd->SetTitle("time resolution, sigma_t vs A");
	hpxnd->GetYaxis()->CenterTitle();
	hpxnd->GetXaxis()->SetTitle("A [m^2]");
	hpxnd->GetXaxis()->CenterTitle();
	
	TGraph *gr_nd = new TGraph (areas, sigma_area);
	gr_nd ->SetMarkerStyle(20);
	gr_nd ->SetMarkerColor(4);
	gr_nd ->SetMarkerSize(1.0);
	gr_nd ->SetLineWidth(2);
	gr_nd ->SetLineColor(4);
	
	gr_nd->Draw("P");



    TCanvas *canv_T = new TCanvas("canv_T", "temperature", 200, 10, 1000, 650);
	canv_T->SetGrid();
	
	TH2F *hpxt = new TH2F("hpxt", "temperature", 300, 90, 650, 100, 0, 2);
	
	hpxt->Draw();
	
	hpxt->GetYaxis()->SetTitle("time resolution [ns]");
	
	
	hpxt->GetYaxis()->SetLabelSize(0.05);
	hpxt->GetXaxis()->SetLabelSize(0.05);
	hpxt->GetYaxis()->SetTitleSize(0.05);
	hpxt->GetXaxis()->SetTitleSize(0.05);
	hpxt->GetXaxis()->SetNoExponent();
	hpxt->SetTitle("time resolution, sigma_t vs T");
	hpxt->GetYaxis()->CenterTitle();
	hpxt->GetXaxis()->SetTitle("T [K]");
	hpxt->GetXaxis()->CenterTitle();
	
	TGraph *gr_T = new TGraph ( Tss, sigma_T);
	gr_T ->SetMarkerStyle(20);
	gr_T ->SetMarkerColor(4);
	gr_T ->SetMarkerSize(1.0);
	gr_T ->SetLineWidth(2);
	gr_T ->SetLineColor(4);
	
	gr_T->Draw("P");


    TCanvas *canv_V = new TCanvas("canv_ V", "Voltage", 200, 10, 1000, 650);
	canv_V->SetGrid();
	
	TH2F *hpxv = new TH2F("hpxv", "Voltage", 20, 200, 2000, 100, 0, 2.5);
	
	hpxv->Draw();
	
	hpxv->GetYaxis()->SetTitle("time resolution [ns]");
	
	
	hpxv->GetYaxis()->SetLabelSize(0.05);
	hpxv->GetXaxis()->SetLabelSize(0.05);
	hpxv->GetYaxis()->SetTitleSize(0.05);
	hpxv->GetXaxis()->SetTitleSize(0.05);
	hpxv->GetXaxis()->SetNoExponent();
	hpxv->SetTitle("time resolution, sigma_t vs V");
	hpxv->GetYaxis()->CenterTitle();
	hpxv->GetXaxis()->SetTitle("U [V]");
	hpxv->GetXaxis()->CenterTitle();
	
	TGraph *gr_V = new TGraph ( Vs, sigma_V);
	gr_V ->SetMarkerStyle(20);
	gr_V ->SetMarkerColor(4);
	gr_V ->SetMarkerSize(1.0);
	gr_V ->SetLineWidth(2);
	gr_V ->SetLineColor(4);
	
	gr_V->Draw("P");
}