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
    ifstream outputC("output_C.txt");
    ifstream outputnoise("output_noise.txt");
    ifstream outputthresh("output_threshold.txt");
    ifstream outputpourcent("output_porcentage.txt");
    

    //vector<ifstream> outputs={output, outputT, outputV};
    vector<string> var={"n_d", "T", "V"};
    string line;
    string line1;
    string line2;
    string line3;
    string line4;
    string line5;
    string line6;

    //Capacities=(40 30 20 10 5)
   // Noises=(1 2 2.5 3 3.2)
    
    
    //TVectorD sigma_nd(2);
    //TVectorD nd(2);
    
    TVectorD sigma_area(5);
    TVectorD sigma_T(6);
    TVectorD sigma_C(5);
    TVectorD sigma_noise(5);
    TVectorD sigma_V(5);
    TVectorD sigma_threshold(7);
    TVectorD sigma_pourcent(7);
    TVectorD Thresholds_CFD(7);

    
    
    double arr[] = {110, 200, 430, 500, 600, 300};
    double arr_V[] = { 300, 1000, 1500, 2000, 600};
    double arr_C[] = { 40, 30, 20, 10, 5};
    double arr_noise[]={1, 2, 2.5, 3, 3.2};
    double arr_thresh[]={2, 3, 4, 5, 7, 8, 9};
    double arr_porc[]={0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};


    TVectorD areas(5);
    TVectorD threshold(7, arr_thresh);
    TVectorD pourcentage(7, arr_porc);
    TVectorD Cs(5, arr_C);
    TVectorD Vs(5, arr_V);
    TVectorD Tss(6, arr);
    TVectorD Noises(5, arr_noise);
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
            //cout<<"sigmas area="<<sigma_area[count]<<endl;
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
           //cout<<"sigmas T="<<sigma_T[count]<<endl;
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
            //cout<<"sigmas V="<<sigma_V[count]<<endl;
            //cout<<"count   "<<count<<endl;
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
    count=0;
    while (getline(outputC, line3)) 
    {      
        if(line3.find("sigma")!=string::npos)
        {
            int pos = line3.find("=");
            // Copy substring after pos
            string sub = line3.substr(pos + 1);
            sigma_C[count]=stod(sub);
            cout<<"sigmas C ="<<sigma_C[count]<<endl;
            cout<<"count   "<<count<<endl;
            count+=1;
        }        
    }
    count=0;
    while (getline(outputnoise, line4)) 
    {      
        if(line4.find("sigma")!=string::npos)
        {
            int pos = line4.find("=");
            // Copy substring after pos
            string sub = line4.substr(pos + 1);
            sigma_noise[count]=stod(sub);
            //cout<<"sigmas ="<<sigma_noise[count]<<endl;
            //cout<<"count   "<<count<<endl;
            count+=1;
        }        
    }
    count=0;
    while (getline(outputthresh, line5)) 
    {      
        if(line5.find("sigma")!=string::npos)
        {
            int pos = line5.find("=");
            // Copy substring after pos
            string sub = line5.substr(pos + 1);
            sigma_threshold[count]=stod(sub);
            count+=1;
        }     
    }
    count=0;
    while (getline(outputpourcent, line6)) 
    {      
        if(line6.find("sigma")!=string::npos)
        {
            int pos = line6.find("=");
            // Copy substring after pos
            string sub = line6.substr(pos + 1);
            sigma_pourcent[count]=stod(sub);
        }   
        if(line6.find("threshold CFD=")!=string::npos)
        {
            int pos = line6.find("=");
            // Copy substring after pos
            string sub = line6.substr(pos + 1);
            Thresholds_CFD[count]=stod(sub);
            count+=1;
        }      
    }
 

    gStyle->SetOptStat(00000000);

    TCanvas *canv_area = new TCanvas("canv_area", "area", 200, 10, 1000, 650);
	canv_area->SetGrid();
	
	TH2F *hpxnd = new TH2F("hpxnd", "area", 10000, areas.Min(), areas.Max(), 50, 0, sigma_area.Max());
	
	hpxnd->Draw();
	
	hpxnd->GetYaxis()->SetTitle("Time Resolution [ns]");
	
	
	hpxnd->GetYaxis()->SetLabelSize(0.05);
	hpxnd->GetXaxis()->SetLabelSize(0.05);
	hpxnd->GetYaxis()->SetTitleSize(0.05);
	hpxnd->GetXaxis()->SetTitleSize(0.05);
	hpxnd->GetXaxis()->SetNoExponent(0);
	hpxnd->SetTitle("Time Resolution vs. Area");
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
	
	TH2F *hpxt = new TH2F("hpxt", "temperature", 300, Tss.Min()-1, Tss.Max()+1, 100, 0, sigma_T.Max());
	
	hpxt->Draw();
	
	hpxt->GetYaxis()->SetTitle("Time Resolution [ns]");
	
	
	hpxt->GetYaxis()->SetLabelSize(0.05);
	hpxt->GetXaxis()->SetLabelSize(0.05);
	hpxt->GetYaxis()->SetTitleSize(0.05);
	hpxt->GetXaxis()->SetTitleSize(0.05);
	hpxt->GetXaxis()->SetNoExponent();
	hpxt->SetTitle("Time Resolution vs. Temperature");
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
	cout<<"T min"<<Tss.Min()<<endl;
    cout<<"sigma T max"<<sigma_T.Max()<<endl;

	TH2F *hpxv = new TH2F("hpxv", "Voltage", 20, Vs.Min()-1, Vs.Max()+1, 100, 0, sigma_V.Max());
	
	hpxv->Draw();
	
	hpxv->GetYaxis()->SetTitle("Time Resolution [ns]");
	
	
	hpxv->GetYaxis()->SetLabelSize(0.05);
	hpxv->GetXaxis()->SetLabelSize(0.05);
	hpxv->GetYaxis()->SetTitleSize(0.05);
	hpxv->GetXaxis()->SetTitleSize(0.05);
	hpxv->GetXaxis()->SetNoExponent();
	hpxv->SetTitle("Time Resolution vs. Voltage");
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



    TCanvas *canv_C = new TCanvas("canv_ C", "Capacity", 200, 10, 1000, 650);
	canv_V->SetGrid();
	
	TH2F *hpxc = new TH2F("hpxc", "Capacity", 20, Cs.Min(), Cs.Max(), 100, 0, sigma_C.Max());
	
	hpxc->Draw();
	
	hpxc->GetYaxis()->SetTitle("Time Resolution [ns]");
	
	
	hpxc->GetYaxis()->SetLabelSize(0.05);
	hpxc->GetXaxis()->SetLabelSize(0.05);
	hpxc->GetYaxis()->SetTitleSize(0.05);
	hpxc->GetXaxis()->SetTitleSize(0.05);
	hpxc->GetXaxis()->SetNoExponent();
	hpxc->SetTitle("Time Resolution vs. Capacity");
	hpxc->GetYaxis()->CenterTitle();
	hpxc->GetXaxis()->SetTitle("C [pF]");
	hpxc->GetXaxis()->CenterTitle();
	
	TGraph *gr_C = new TGraph ( Cs, sigma_C);
	gr_C ->SetMarkerStyle(20);
	gr_C ->SetMarkerColor(4);
	gr_C ->SetMarkerSize(1.0);
	gr_C ->SetLineWidth(2);
	gr_C ->SetLineColor(4);
	
	gr_C->Draw("P");


    TCanvas *canv_noise = new TCanvas("canv_ noise", "Noise", 200, 10, 1000, 650);
	canv_noise->SetGrid();
	
	TH2F *hpxnoise = new TH2F("hpxnoise", "Noise", 20, 1, 4, 100, 0, sigma_noise.Max());
	
	hpxnoise->Draw();
	
	hpxnoise->GetYaxis()->SetTitle("Time Resolution [ns]");
	
	
	hpxnoise->GetYaxis()->SetLabelSize(0.05);
	hpxnoise->GetXaxis()->SetLabelSize(0.05);
	hpxnoise->GetYaxis()->SetTitleSize(0.05);
	hpxnoise->GetXaxis()->SetTitleSize(0.05);
	hpxnoise->GetXaxis()->SetNoExponent();
	hpxnoise->SetTitle("Time Resolution vs. Noise");
	hpxnoise->GetYaxis()->CenterTitle();
	hpxnoise->GetXaxis()->SetTitle("I [mA]");
	hpxnoise->GetXaxis()->CenterTitle();
	
	TGraph *gr_noise = new TGraph(Noises, sigma_noise);
	gr_noise ->SetMarkerStyle(20);
	gr_noise ->SetMarkerColor(4);
	gr_noise ->SetMarkerSize(1.0);
	gr_noise ->SetLineWidth(2);
	gr_noise ->SetLineColor(4);
	
	gr_noise->Draw("P");

    TCanvas *canv_thresh = new TCanvas("canv_ thresh", "Threshold", 200, 10, 1000, 650);
	canv_noise->SetGrid();
	
	TH2F *hpxthresh = new TH2F("hpxthresh", "Threshold", 20, threshold.Min(), threshold.Max(), 100, 0, sigma_threshold.Max());
	
	hpxthresh->Draw();
	
	hpxthresh->GetYaxis()->SetTitle("Time Resolution [ns]");
	
	
	hpxthresh->GetYaxis()->SetLabelSize(0.05);
	hpxthresh->GetXaxis()->SetLabelSize(0.05);
	hpxthresh->GetYaxis()->SetTitleSize(0.05);
	hpxthresh->GetXaxis()->SetTitleSize(0.05);
	hpxthresh->GetXaxis()->SetNoExponent();
	hpxthresh->SetTitle("Time Resolution vs. Threshold");
	hpxthresh->GetYaxis()->CenterTitle();
	hpxthresh->GetXaxis()->SetTitle("I [mA]");
	hpxthresh->GetXaxis()->CenterTitle();
	
	TGraph *gr_thresh = new TGraph ( threshold, sigma_threshold);
	gr_thresh ->SetMarkerStyle(20);
	gr_thresh ->SetMarkerColor(4);
	gr_thresh ->SetMarkerSize(1.0);
	gr_thresh ->SetLineWidth(2);
	gr_thresh ->SetLineColor(4);
	
	gr_thresh->Draw("P");

    TCanvas *canv_porc = new TCanvas("canv_ porc", "CFD", 200, 10, 1000, 650);
	canv_porc->SetGrid();
	
	TH2F *hpxtcfd = new TH2F("hpxtcfd", "CFD", 20, Thresholds_CFD.Min(), Thresholds_CFD.Max(), 100, 0, sigma_pourcent.Max());
	
	hpxtcfd->Draw();
	
	hpxtcfd->GetYaxis()->SetTitle("Time Resolution [ns]");
	
	
	hpxtcfd->GetYaxis()->SetLabelSize(0.05);
	hpxtcfd->GetXaxis()->SetLabelSize(0.05);
	hpxtcfd->GetYaxis()->SetTitleSize(0.05);
	hpxtcfd->GetXaxis()->SetTitleSize(0.05);
	hpxtcfd->GetXaxis()->SetNoExponent();
	hpxtcfd->SetTitle("Time Resolution vs. Threshold CFD");
	hpxtcfd->GetYaxis()->CenterTitle();
	hpxtcfd->GetXaxis()->SetTitle("I [mA]");
	hpxtcfd->GetXaxis()->CenterTitle();
	
	TGraph *gr_cfd = new TGraph ( Thresholds_CFD, sigma_pourcent);
	gr_cfd ->SetMarkerStyle(20);
	gr_cfd ->SetMarkerColor(4);
	gr_cfd ->SetMarkerSize(1.0);
	gr_cfd ->SetLineWidth(2);
	gr_cfd ->SetLineColor(4);
	
	gr_cfd->Draw("P");
    
}