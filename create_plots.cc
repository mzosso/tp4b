#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <typeinfo>
#include <TRandom.h>
using namespace std;

void create_plots()
{
    string line;
    ifstream output("output_n_d.txt");
    TVectorD sigma_nd(2);
    TVectorD nd(2);
    while (getline(output, line)) 
    {
        //cout<<"line="<<line<<endl;
        //stringstream ss(line);
       // string valueString;
        //ss >> valueString;
        //cout<<"valueString ="<<valueString<<endl;
        int count=0;
        if(line.find("sigma")!=string::npos)
        {
            int pos = line.find("=");
            // Copy substring after pos
            string sub = line.substr(pos + 1);
            sigma_nd[count]=stod(sub);
            cout<<"sigmas ="<<stod(sub)<<endl;
        }
        if(line.find("n_d")!=string::npos)
        {
            int pos = line.find("=");
            // Copy substring after pos
            string sub = line.substr(pos + 1);
            nd[count]=stod(sub);
            cout<<"n_ds ="<<stod(sub)<<endl;
        }
        count+=1;
    }
    TCanvas *canv_nd = new TCanvas("canv_ nd", "thickness", 200, 10, 1000, 650);
	canv_nd->SetGrid();
	
	TH2F *hpxnd = new TH2F("hpxnd", "thickness", 20, 0, 300, 100, 0, 1e1);
	
	hpxnd->Draw();
	
	hpxnd->GetYaxis()->SetTitle("time resolution [ns]");
	
	
	hpxnd->GetYaxis()->SetLabelSize(0.05);
	hpxnd->GetXaxis()->SetLabelSize(0.05);
	hpxnd->GetYaxis()->SetTitleSize(0.05);
	hpxnd->GetXaxis()->SetTitleSize(0.05);
	hpxnd->GetXaxis()->SetNoExponent();
	hpxnd->SetTitle("time resolution, sigma_t vs d");
	hpxnd->GetYaxis()->CenterTitle();
	hpxnd->GetXaxis()->SetTitle("d [micro m]");
	hpxnd->GetXaxis()->CenterTitle();
	
	TGraph *gr_nd = new TGraph ( nd, sigma_nd);
	gr_nd ->SetMarkerStyle(20);
	gr_nd ->SetMarkerColor(4);
	gr_nd ->SetMarkerSize(1.0);
	gr_nd ->SetLineWidth(2);
	gr_nd ->SetLineColor(4);
	
	gr_nd->Draw("LP");
}