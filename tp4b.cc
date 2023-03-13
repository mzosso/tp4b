#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <typeinfo>
using namespace std;

class Simple
{
	protected:
	double dopage_p;
	double dopage_n;
	double depleted;
	double voltage;
	double e;
	double m;
	
	public:
	Simple(double dopage_p, double dopage_n, double depleted, double voltage, double e, double m) : dopage_p(dopage_p), dopage_n(dopage_n), depleted(depleted), voltage(voltage),
	e(1.6*pow(10,-19)), m(9.1 * pow(10, -31)){}
	double drift_vel_p(double E, double t) const {return e*E*t/m;};
};
class Temperature: public Simple {
	
	private:
	double temperature;
	
	public:
	Temperature(double dopage_p, double dopage_n, double depleted, double voltage,double e, double m, double temperature) : 
	 Simple(dopage_p, dopage_n, depleted, voltage, e, m), temperature(temperature) {}
	 
};
	
int n=1000; //how many pts for d and E
double n_d=300; //thickness
double dx=(double) n_d/n; //dx of E and d
//int V=300; //external
int s=100; //how many ns
int n_t=1000; //how many pts for s and I
double dt= (double) s/n_t;
int V_d=100; //depletion voltage

vector<double> tps(double ini, double n_d, double mu_e,double mu_h, float E[n]) //ini c'est où on injecte // drift time for h and e
{
	vector<double> t(2);
	double dist=n_d-ini; //distance à parcourir pour l'électron
	int i_e=round(n*dist/n_d); // avec d[i]=dist=n_d/n*i
	int i_h=round(n*ini/n_d); // ATTENTION, ROUND DONC PAS E EXACT
	double v_e=mu_e*E[i_e];
	double v_h=mu_h*E[i_h];
	double t_e=dist/v_e;
	double t_h=ini/v_h;
	cout<<"t_e="<<t_e<<endl;
	cout<<"t_h="<<t_h<<endl;
	t[0]=t_e;
	t[1]=t_h;
	return t;
}

void fct_E(float E[n], float d[n],int V_d, int V, bool cst=0) //pour pouvoir définir un champ élec constant facilement
{
		
	
		/*vector<float E[n]> E_sup(nb);
		vector<double> V_sup(nb);*/
		
		double a=(double)-V_d/n_d; //pente du champ E
		for(int i(0); i<n;++i)
		{
			if (cst==0)
				{E[i]=V;}
			else
			{
				/*d.push_back(n_d/n*(n-i));
				E.push_back(V/d[i]);*/
				//d[i]= (double) n_d/n*(n-i);
				d[i] = (double) n_d/n*i;
				if (V_d<=V)
				{
					double dif=(V-V_d);
					E[i]=(double) V_d+dif+a*d[i];
					
				}
				else
				{
					double dif=(V-V_d)/n_d;
					E[i]=(double) (dif+2.0/n_d*V_d*(1-d[i]/double(n_d)));
				}
				/*for (in j(0); j<=nb;++j)
				{
					V_sup[i]=V_d+50-j*nb;
					E_sup[i][j]=(double) V_sup[i]+a*d[i];
				}*/
			}    
	}
}

void fct_I(float I[n], float t[n],double dt,double tps, bool cst=0) 
{
	tps*=1e9;
	for(int i(0); i<n_t;++i)
		{
			//t.push_back(s/n_t*i);
			t[i] = (double) s/n_t*i;
		}
		
	if (cst==0)
	{
		TRandom *rand = new TRandom();
		
		
		double r = 0;
		for(int i(0); i<n_t;++i)
		{
			//I.push_back(0);
			r = rand->Gaus(0,1);
			I[i]=r;
		}
	}
	else
	{
		cout<<"tps"<<tps<<endl;
		cout<<"dt"<<dt<<endl;
		
		int num_dt=round(tps/dt); //nombre de dt
		cout<<"num_dt"<<num_dt<<endl;
		double height_I=(double) 1/tps; // aire du rectangle =1
		cout<<height_I<<endl;
		for(int j(0); j<num_dt; ++j)
		{
			I[j]=height_I;
			
		}
		for(int g(num_dt); g<n_t;++g)
		{
			I[g]=0;
		}
		
	}
}

void tp4b(int V=300)
{
	
	/*vector<double> E(n);
	vector<double> d(n);
	vector<double> t(n_t);
	vector<double> I(n_t);*/
	float E[n];
	float d[n];
	float t[n_t];
	float I_e[n_t];
	float I_h[n_t];
	
	
	
	fct_E(E, d, V_d,V, 1);
	
	double ini=200;
	double mu_e=1350/pow(10,-4);
	double mu_h=450/pow(10,-4);
	vector<double> temps=tps(ini, n_d, mu_e, mu_h, E);
	fct_I(I_e, t, dt, temps[0], 1);
	fct_I(I_h, t, dt, temps[1], 1);
	
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
	
	TH2F *hpx0 = new TH2F("hpx0", "Electric", 20, 0, 300, 100, 0, 1e3);
	
	hpx0->Draw();
	
	hpx0->GetYaxis()->SetTitle("E");
	
	
	hpx0->GetYaxis()->SetLabelSize(0.05);
	hpx0->GetXaxis()->SetLabelSize(0.05);
	hpx0->GetYaxis()->SetTitleSize(0.05);
	hpx0->GetXaxis()->SetTitleSize(0.05);
	hpx0->GetXaxis()->SetNoExponent();
	hpx0->SetTitle("E, d vs E");
	hpx0->GetYaxis()->CenterTitle();
	hpx0->GetXaxis()->SetTitle("d");
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
	
	TMultiGraph *mg = new TMultiGraph();
	mg->SetTitle("Signals");

	
	
	/*TH2F *hpx1 = new TH2F("hpx1", "I", 20, 0, s, 100, -1, 1);
	hpx1->Draw("apl");
	
	
	hpx1->GetYaxis()->SetTitle("I");

	hpx1->GetYaxis()->SetLabelSize(0.05);
	hpx1->GetXaxis()->SetLabelSize(0.05);
	hpx1->GetYaxis()->SetTitleSize(0.05);
	hpx1->GetXaxis()->SetTitleSize(0.05);
	hpx1->GetXaxis()->SetNoExponent();
	hpx1->SetTitle("I, ns vs I");
	hpx1->GetYaxis()->CenterTitle();
	hpx1->GetXaxis()->SetTitle("ns");
	hpx1->GetXaxis()->CenterTitle();*/
	
	
	
	TGraph *gr_I_e = new TGraph (n_t, t, I_e);
	gr_I_e ->SetMarkerStyle(20);
	gr_I_e ->SetMarkerColor(2);
	gr_I_e ->SetMarkerSize(1.0);
	gr_I_e ->SetLineWidth(2);
	gr_I_e ->SetLineColor(2);
	
	
	TGraph *gr_I_h = new TGraph (n_t, t, I_h);
	gr_I_h ->SetMarkerStyle(20);
	gr_I_h ->SetMarkerColor(2);
	gr_I_h ->SetMarkerSize(1.0);
	gr_I_h ->SetLineWidth(2);
	gr_I_h ->SetLineColor(5);
	
	mg->Add(gr_I_h);
	mg->Add(gr_I_e);	
	mg->Draw("AC");
	
	
	
	
}
