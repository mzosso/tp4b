#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <typeinfo>
using namespace std;


	
const int n=1000; //how many pts for d and E
double area=1e5; //m^2
double n_d=100; //thickness, micro
double dx= (double) n_d/n;
int s=100; //how many ns
int n_t=2000; //how many pts for s and I
double dt= (double) s/n_t; //ns
int V_d=300; //depletion voltage
int n_e=100; //nb electrons
double T=298.16; //K
double mu_e=1350; //cm^2/(Vs)
double mu_h=480; //cm^2/(Vs)
double k_B=8.62e-5; //eV/K

static Float_t vel_h_at_position_in_field(float value_of_field, double temp);

Float_t* vel_e( float E_[n], double T_);
double diffusion_coefficient(double T_, double mu, double k_B);
Float_t diffusion_velocity(double mu, int i, Float_t E[n]);

Float_t* vel_h(float E_[n], double T_);

Float_t vel_h_at_position_in_field(float value_of_field, double temp);
vector<vector<double>> tps(float ini[n_e], double n_d,float E[n], double T_); //ini c'est où on injecte // drift time for h and e


void fct_E(float E[n], float d[n],int V_d, int V, bool cst=0); //pour pouvoir définir un champ élec constant facilement
Float_t get_gradient(Float_t E[n], int i);
void add_noise(Float_t I[n_t], double std);
void fct_I(float I[n_t], float t[n_t],double dt,vector<double> tps);

void apply_filter_time_domain(float I[n_t], double n_d);


void tp4b(int V=500)
{
	const int N = 1000;
    float E_[N];
	for(int l = 0; l < N; l++)
	{
		E_[l]=0;
	}
    const double start = 2.0;
    const double end = 5.0;
    const double step = (end - start) / (N - 1);
    for (int i = 0; i < N; i++) {
        E_[i] = pow(10, start + i * step);
		
    }

	float E[n];
	float d[n];
	float t[n_t];
	float I_e[n_t];
	float I_h[n_t];
	float I_tot[n_t];
	float I_tot1[n_t];
	float ini[n_e]; //vecteur d'elec pour pouvoir les placer dans le detecteur
	

	for(int i(0); i<n_t; ++i)
	{
		t[i]=0;
		I_h[i]=0;
		I_e[i]=0;
		I_tot[i]=0;
		I_tot1[i]=0;
	}
	for(int k(0); k<n; ++k)
	{
		E[k]=0;
		d[k]=0;
	}
	for(int c(0); c<n_e;++c)
	{
		ini[c]=0;
	}
	
	
	fct_E(E, d, V_d,V, 1);

	//double mu_e=1350/pow(10,-4);
	//double mu_h=450/pow(10,-4);
	vector<vector<double>> temps=tps(ini, n_d, E,T);

	fct_I(I_e, t, dt, temps[0]);
	fct_I(I_h, t, dt, temps[1]);
	
	for(int i(0); i<n_t;++i)
	{
		I_tot[i]=I_h[i]+I_e[i];
		I_tot1[i]=I_h[i]+I_e[i];
		
	}
	
	while (gPad !=0) gPad ->Close();
	
	gStyle -> SetPadLeftMargin(0.15);
	gStyle -> SetPadRightMargin(0.06);
	gStyle -> SetPadTopMargin(0.12);
	gStyle -> SetPadBottomMargin(0.14);
	gStyle -> SetTitleFontSize(0.06);
	TGaxis::SetMaxDigits(4);
	
	gStyle->SetOptStat(00000000);
	
	TCanvas *canv0 = new TCanvas("canv0", "E", 200, 10, 1000, 650);
	canv0->SetGrid();
	
	TH2F *hpx0 = new TH2F("hpx0", "Electric", 20, 0, 300, 100, 0, 1e1);
	
	hpx0->Draw();
	
	hpx0->GetYaxis()->SetTitle("E [V/micro m]");
	
	
	hpx0->GetYaxis()->SetLabelSize(0.05);
	hpx0->GetXaxis()->SetLabelSize(0.05);
	hpx0->GetYaxis()->SetTitleSize(0.05);
	hpx0->GetXaxis()->SetTitleSize(0.05);
	hpx0->GetXaxis()->SetNoExponent();
	hpx0->SetTitle("Electric field, E vs thickness");
	hpx0->GetYaxis()->CenterTitle();
	hpx0->GetXaxis()->SetTitle("d [micro m]");
	hpx0->GetXaxis()->CenterTitle();
	
	TGraph *gr_E = new TGraph (n, d, E);
	gr_E ->SetMarkerStyle(20);
	gr_E ->SetMarkerColor(4);
	gr_E ->SetMarkerSize(1.0);
	gr_E ->SetLineWidth(2);
	gr_E ->SetLineColor(4);
	
	gr_E->Draw("LP");
	
	
	//Section different signals on same figure
	
	TCanvas *canv1 = new TCanvas("canv1", "I", 200, 10, 1000, 650);
	canv1->SetGrid();
	TH2F *hpx1 = new TH2F("hpx1", "I", 20, 0, 1e8*s, 100, -1, 70);
	hpx1->Draw();
	hpx1->GetYaxis()->SetTitle("I");
	hpx1->GetYaxis()->SetLabelSize(0.05);
	hpx1->GetXaxis()->SetLabelSize(0.05);
	hpx1->GetYaxis()->SetTitleSize(0.05);
	hpx1->GetXaxis()->SetTitleSize(0.05);
	hpx1->GetXaxis()->SetNoExponent();
	hpx1->SetTitle("I, Current vs time");
	hpx1->GetYaxis()->CenterTitle();
	hpx1->GetXaxis()->SetTitle("t [ns]");
	hpx1->GetXaxis()->CenterTitle();
	

	TGraph *gr_I_e = new TGraph (n_t, t, I_e);
	gr_I_e ->SetMarkerStyle(20);
	gr_I_e ->SetMarkerColor(2);
	gr_I_e ->SetMarkerSize(1.0);
	gr_I_e ->SetLineWidth(2);
	gr_I_e ->SetLineColor(2);
	
	
	
	TGraph *gr_I_h = new TGraph (n_t, t, I_h);
	gr_I_h ->SetMarkerStyle(20);
	gr_I_h ->SetMarkerColor(4);
	gr_I_h ->SetMarkerSize(1.0);
	gr_I_h ->SetLineWidth(2);
	gr_I_h ->SetLineColor(4);
	
	
	TGraph *gr_I_tot = new TGraph (n_t, t, I_tot);
	gr_I_tot ->SetMarkerStyle(20);
	gr_I_tot ->SetMarkerColor(1);
	gr_I_tot ->SetMarkerSize(1.0);
	gr_I_tot ->SetLineWidth(2);
	gr_I_tot ->SetLineColor(1);

	auto legend = new TLegend(); //0.2,0.3,0.2,0.3
	vector<TString> mylgd ={"I_{e}","I_{h}", "I_{tot}"};
	legend->AddEntry(gr_I_e,mylgd[0], "l");
	legend->AddEntry(gr_I_h,mylgd[1], "l");
	legend->AddEntry(gr_I_tot,mylgd[2], "l");



	gr_I_e->Draw("L");
	gr_I_h->Draw("L");
	gr_I_tot->Draw("L");
	legend->Draw();
	

	///Noise added after filter
	apply_filter_time_domain(I_tot, n_d);
	add_noise(I_tot, 1.0);
	TCanvas *canvfilter = new TCanvas("canvfilter", "I, noise added after filtering", 200, 10, 1000, 650);
	canvfilter->SetGrid();
	TH2F *hpxfilter = new TH2F("hpxfilter", "I", 20, 0, 1e8*s, 100, -1, 36);
	hpxfilter->Draw();
	hpxfilter->GetYaxis()->SetTitle("I ");
	hpxfilter->GetYaxis()->SetLabelSize(0.05);
	hpxfilter->GetXaxis()->SetLabelSize(0.05);
	hpxfilter->GetYaxis()->SetTitleSize(0.05);
	hpxfilter->GetXaxis()->SetTitleSize(0.05);
	hpxfilter->GetXaxis()->SetNoExponent();
	hpxfilter->SetTitle("I_{tot}, Current vs time, filtered");
	hpxfilter->GetYaxis()->CenterTitle();
	hpxfilter->GetXaxis()->SetTitle("t [ns]");
	hpxfilter->GetXaxis()->CenterTitle();
	TGraph *gr_filter = new TGraph (n_t, t, I_tot);
	gr_filter->SetMarkerStyle(20);
	gr_filter ->SetMarkerColor(2);
	gr_filter ->SetMarkerSize(1.0);
	gr_filter ->SetLineWidth(2);
	gr_filter ->SetLineColor(2);
	
	gr_filter->Draw("L");

///Noise added before filter
	add_noise(I_tot1, 1.0);
	apply_filter_time_domain(I_tot1, n_d);
	TCanvas *canvfilter_ = new TCanvas("canvfilter_", "I, noise added before filtering", 200, 10, 1000, 650);
	canvfilter_->SetGrid();
	TH2F *hpxfilter_ = new TH2F("hpxfilter_", "I", 20, 0, 1e8*s, 100, -1, 36);
	hpxfilter_->Draw();
	hpxfilter_->GetYaxis()->SetTitle("I ");
	hpxfilter_->GetYaxis()->SetLabelSize(0.05);
	hpxfilter_->GetXaxis()->SetLabelSize(0.05);
	hpxfilter_->GetYaxis()->SetTitleSize(0.05);
	hpxfilter_->GetXaxis()->SetTitleSize(0.05);
	hpxfilter_->GetXaxis()->SetNoExponent();
	hpxfilter_->SetTitle("I_{tot}, Current vs time, filtered");
	hpxfilter_->GetYaxis()->CenterTitle();
	hpxfilter_->GetXaxis()->SetTitle("t [ns]");
	hpxfilter_->GetXaxis()->CenterTitle();

	TGraph *gr_filter_ = new TGraph (n_t, t, I_tot1);
	gr_filter_->SetMarkerStyle(20);
	gr_filter_ ->SetMarkerColor(2);
	gr_filter_ ->SetMarkerSize(1.0);
	gr_filter_ ->SetLineWidth(2);
	gr_filter_ ->SetLineColor(2);
	
	gr_filter_->Draw("L");

	
	vector<double> T_={6.0,45.0,110.0,200.0,300.0,430.0};
	size_t nt=T_.size();
	vector<Float_t*> ve(nt);
	vector<Float_t*> vh(nt);
	for(int i;i<nt;++i)
	{
		cout<<T_[i]<<endl;
		ve[i]=vel_e(E_, T_[i]);
		vh[i]=vel_h(E_, T_[i]);
	}
	Float_t* ve_=vel_e(E_, 100);
	Float_t* vh_=vel_h(E_, 100);
	
	TCanvas *canvas2 = new TCanvas("canvas2", "Drift velocity vs E", 200, 10, 1000, 650);
    TH2F *hpx2 = new TH2F("hpx2", "Electric", 100,1e2 ,1e5, n,1e4 , 1e7); //nb binsx, xlow, xup, nbinsy, ydown, yup 
	canvas2->SetGrid();
	hpx2->Draw();
	hpx2->GetYaxis()->SetTitle("v [cm/s]");
	hpx2->GetYaxis()->SetLabelSize(0.05);
	hpx2->GetXaxis()->SetLabelSize(0.05);
	hpx2->GetYaxis()->SetTitleSize(0.05);
	hpx2->GetXaxis()->SetTitleSize(0.05);
	hpx2->GetXaxis()->SetNoExponent();
	hpx2->SetTitle("Drift velocity, v vs E");
	hpx2->GetYaxis()->CenterTitle();
	hpx2->GetXaxis()->SetTitle("E [V/cm]");
	hpx2->GetXaxis()->CenterTitle();
	gPad->SetLogx();
	gPad->SetLogy();
	
	TGraph *graph_Vel = new TGraph(n, E_, ve_);
    graph_Vel->SetTitle("Electron velocity vs E");
    graph_Vel->GetXaxis()->SetTitle("E");
    graph_Vel->GetYaxis()->SetTitle("v [cm/s]");
    graph_Vel ->SetMarkerStyle(20);
	graph_Vel ->SetMarkerColor(2);
	graph_Vel ->SetMarkerSize(1.0);
	graph_Vel ->SetLineWidth(2);
	graph_Vel ->SetLineColor(2);

	TGraph *graph_Velh = new TGraph(n, E_, vh_);
    graph_Velh->SetTitle("Hole velocity vs E");
    graph_Velh->GetXaxis()->SetTitle("E");
    graph_Velh->GetYaxis()->SetTitle("v [cm/s]");
    graph_Velh ->SetMarkerStyle(20);
	graph_Velh ->SetMarkerColor(4);
	graph_Velh ->SetMarkerSize(1.0);
	graph_Velh ->SetLineWidth(2);
	graph_Velh ->SetLineColor(4);

	auto legend1 = new TLegend(); //0.2,0.3,0.2,0.3
	vector<TString> mylgd1 ={"v_{e}","v_{h}"};
	legend1->AddEntry(graph_Vel,mylgd1[0], "l");
	legend1->AddEntry(graph_Velh,mylgd1[1], "l");
	

	graph_Vel->Draw("L");
	graph_Velh->Draw("L");
	legend1->Draw();


	/*// plot different temperatures 
	
    TH2F *hpx_T = new TH2F("hpx_T", "Electric", 100,1e2 ,1e5, n,1e4 , 1e7); //nb binsx, xlow, xup, nbinsy, ydown, yup 
	
	hpx_T->Draw();
	hpx_T->GetYaxis()->SetTitle("v [cm/s]");
	hpx_T->GetYaxis()->SetLabelSize(0.05);
	hpx_T->GetXaxis()->SetLabelSize(0.05);
	hpx_T->GetYaxis()->SetTitleSize(0.05);
	hpx_T->GetXaxis()->SetTitleSize(0.05);
	hpx_T->GetXaxis()->SetNoExponent();
	hpx_T->SetTitle("Drift velocity, v vs E");
	hpx_T->GetYaxis()->CenterTitle();
	hpx_T->GetXaxis()->SetTitle("E [V/cm]");
	hpx_T->GetXaxis()->CenterTitle();
	gPad->SetLogx();
	gPad->SetLogy();

	//vector<TGraph*> graph_T(T_.size());

	auto legend2 = new TLegend(); //0.2,0.3,0.2,0.3
	vector<TString> mylgd2 ={"T=6K","T=45K", "T=110K", "T=200K", "T=300K", "T=430K"};
	vector<int> colors = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange};
	//vector<double> color={1,2,3,4,6,9};
	
	vector<TGraph*> graphs;
for(int i = 0; i < T_.size(); ++i)
{
    auto graph_T = new TGraph(n, E_, ve[i]);
    graph_T->SetTitle("Electron velocity vs E");
    graph_T->GetXaxis()->SetTitle("E");
    graph_T->GetYaxis()->SetTitle("v [cm/s]");
    graph_T->SetLineWidth(2);
    graph_T->SetLineColor(colors[i]);

    legend2->AddEntry(graph_T, mylgd2[i], "l");
    graphs.push_back(graph_T);
}

TCanvas *canvas_T = new TCanvas("canvas", "canvas", 800, 600);
canvas_T->SetGrid();
canvas_T->cd();
graphs[0]->Draw("ALP");
for(int i = 1; i < 2; ++i)
{
    graphs[i]->Draw("LSAME");
}

legend2->Draw();
canvas_T->Update();
}*/

//plot histogram noise
Float_t I_noise[n_t];
for(int i(0); i<n_t; ++i)
{
	I_noise[i]=0;
}
	add_noise(I_noise, 1.0);
	
	TH1F *hist = new TH1F("hist", "I_noise Histogram", 100, -10, 10);

		for(int i(0); i<n_t; ++i)
		{
			hist->Fill(I_noise[i]);
		}
	hist->GetXaxis()->SetLogx(0);
hist->Draw();


}

void fct_E(float E[n], float d[n],int V_d, int V, bool cst=0) //pour pouvoir définir un champ élec constant facilement
{
		
		double a=(double)-V_d/n_d; //pente du champ E, V/micrometre
		for(int i(0); i<n;++i)
		{
			if (cst==0)
				{E[i]=V/n_d;}
			else
			{
				d[i] = (double) abs(n_d/n*i); //micrometre
				E[i]=(double) (a*d[i]+V)/n_d; // V/micrometre
			}    
		}
}
double diffusion_coefficient(double T_, double mu, double k_B) 
{
    return k_B * T_ * mu;
}
Float_t get_gradient(Float_t E[n], int i)
{
	return (E[i+1]-E[i])/dx;
}
Float_t diffusion_velocity(double mu, int i, Float_t E[n], double T_)
{
	return diffusion_coefficient(T_,mu,k_B) / mu * get_gradient(E, i); 
}
	
Float_t* vel_e( float E_[n], double T_)
{
	static Float_t v_e[n];
	for(int i(0); i<n; ++i)
	{	double num=1.42e9*E_[i]*pow(T_,-2.42);
		double den = pow(1.0+pow(E_[i]/1.01/pow(T_,1.55),(2.57e-2*pow(T_,0.66))),(1/(2.57e-2*pow(T_,0.66))));
		v_e[i]=num/den;//without diffusion
		//v_e[i]+=diffusion_velocity(mu_e, i, E_); //with diffusion

    }
	return v_e; //cm/s
}

Float_t* vel_h(float E_[n], double T_)
{
	
	double pot_h=0.46*pow(T_, 0.17);
	static Float_t v_h[n];

	transform(E_, E_+n, v_h, [T_](float E){return vel_h_at_position_in_field(E,T_);});
	return v_h;
}

Float_t vel_h_at_position_in_field(float value_of_field, double temp)
{
	double pot_h=0.46*pow(temp, 0.17);
	double num_=1.31e8*value_of_field*pow(temp,-2.2);
	double den_ = pow(1.0+pow(value_of_field/1.24/pow(temp,1.68),pot_h),(1/pot_h));
		
	return num_/den_;
}

vector<vector<double>> tps(float ini[n_e], double n_d,float E[n], double T_) //ini c'est où on injecte // drift time for h and e
{
	vector<double> v_e(n_e);
	vector<double> v_h(n_e);
	vector<double> dist(n_e);
	vector<int> i_e(n_e);
	vector<int> i_h(n_e);
	vector<double> t_e(n_e);
	vector<double> t_h(n_e);
	vector<vector<double>> t(2, vector<double>(n_e));
	double pot_h=0.46*pow(T_, 0.17);
	
	for(int i(0); i<n_e; ++i)
	{
		ini[i]=(double) n_d/n_e*i+n_d/n_e*0.5; // micrometre
		dist[i]=n_d-ini[i]; //distance à parcourir pour l'électron, micrometre
		i_e[i]=floor(abs(n*dist[i]/n_d)); // avec d[i]=dist=n_d/n*i

		i_h[i]=floor(abs(n*ini[i]/n_d)); // ATTENTION, ROUND DONC PAS E EXACT
		double num=1.42e9*E[i_e[i]]/n_d*1e4*pow(T_,-2.42);
		double den = pow(1.0+pow(E[i_e[i]]/n_d*1e4/1.01/pow(T,1.55),(2.57e-2*pow(T_,0.66))),(1/(2.57e-2*pow(T_,0.66))));
		v_e[i]=num/den;

		double num_=1.31e8*E[i_h[i]]/n_d*1e4*pow(T_,-2.2);
		double den_ = pow(1.0+pow(E[i_h[i]]/n_d*1e4/1.24/pow(T_,1.68),pot_h),(1/pot_h));
		
		v_h[i]=num_/den_;
	
		//v_e[i]=mu_e*E[i_e[i]];
		//v_h[i]=mu_h*E[i_h[i]];
		//t_e[i]=dist[i]/v_e[i];
		//t_h[i]=ini[i]/v_h[i];
		
		t[0][i]=dist[i]*1e4/v_e[i]; //micrometre*1e4/(cm/s)=s
		t[1][i]=ini[i]*1e4/v_h[i];
	}
	return t;
}
void add_noise(Float_t I[n_t], double std)
{
	TRandom *rand = new TRandom();
	vector<double> r(n_t);
		for(int k(0); k<n_t; ++k)
		{
			r[k] = rand->Gaus(0,std);
			I[k]+=r[k];
		}
}
void fct_I(float I[n_t], float t[n_t],double dt,vector<double> tps) 
{
	for(int i(0); i<n_t;++i)
		{
			t[i] = (double) 1e9*s/n_t*i;//ns*1e9=s
		}
		
		
		
		
		vector<int> num_dt(n_e);
		vector<double> height_I(n_e);
		
				for(int i(0); i<n_e;++i)
				{
					//tps in seconds, dt in nanoseconds
					num_dt[i]=floor(abs(tps[i]/(dt)));
					
					height_I[i]=(double) 1.0/tps[i]; // aire du rectangle =1
			
					for(int j(20); j<=num_dt[i]; ++j)
					{
						
						int idx = j;
						if (idx >= n_t) 
						{
							idx = n_t - 1;
						}
						I[idx] += height_I[i];
						
					}
					
				}		
}	

void apply_filter_time_domain(float I[n_t], double n_d)
{
	double R=50;
	double C=(double) area/(n_d*pow(10,-6))*8.854e-12; 
	cout<<C<<endl;
	double alpha=dt / (R*C + dt);
	
	for(int i(1); i<n_t; ++i)
	{
		I[i] = I[i-1] + alpha * (I[i] - I[i-1]); // Esteban: Filtro RC, el valor 0.1 dependera de la capacidad de nuestro detector, ver: https://en.wikipedia.org/wiki/Low-pass_filter
		// en nuestro caso: dt = 50 ps, R = 50 ohm y C = 10 pF (por ejemplo)
		// nos da una frecuencia de corte de 300 MHz, esto cambia con la C del detector.
		//I[i]=I[i]*filter[i];
	}
}