#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <typeinfo>
using namespace std;


	
const int n=1000; //how many pts for d and E
double n_d=300; //thickness, micro
double dx=(double) n_d/n; //dx of E and d //micro
//int V=300; //external
int s=100; //how many ns
int n_t=1000; //how many pts for s and I
double dt= (double) s/n_t; //ns
int V_d=100; //depletion voltage
int n_e=10; //nb electrons
double T=298.16; //K

/*const vector<vector<Float_t>> vel(double n_d, float E[n], double T)
{
	vector<double> dist(n);
	vector<double> ini(n);
	vector<int> i_e(n);
	vector<int> i_h(n);
	vector<Float_t> v_e(n);
	vector<Float_t> v_h(n);
	vector<vector<Float_t>> v(2, vector<Float_t>(n));
	for(int i(0); i<n; ++i)
	{
		ini[i]=(double) n_d/n_e*i+n_d/n_e*0.5; //a corregir
		dist[i]=n_d-ini[i]; //distance à parcourir pour l'électron
		i_e[i]=round(n*dist[i]/n_d); // avec d[i]=dist=n_d/n*i
		i_h[i]=round(n*ini[i]/n_d); // ATTENTION, ROUND DONC PAS E EXACT
		v_e[i]=(1.42*pow(10,9)*pow(T, -2.42)*E[i_e[i]])/pow((pow((1+(E[i_e[i]]/1.01)*pow(T,1.55)),2.57*pow(10,-2))*pow(T, 0.66)),(1/2.57*pow(10,-2)*pow(T,0.66)));
		v_h[i]=(1.31*pow(10,8)*pow(T, -2.2)*E[i_h[i]])/pow((pow((1+(E[i_h[i]]/1.24)*pow(T,1.68)),0.46)*pow(T, 0.17)),(1/0.46*pow(T,0.17)));
		v[0][i]=v_e[i];
		v[1][i]=v_h[i];
	}
	
	return v;
}
*/

Float_t* vel_e( float E_[n], double T)
{
	static Float_t v_e[n];
	for(int i(0); i<n; ++i)
	{
		v_e[i]=(1.42*pow(10,9)*pow(T, -2.42)*E_[i])/pow((pow((1+(E_[i]/1.01)*pow(T,1.55)),2.57*pow(10,-2))*pow(T, 0.66)),(1/(2.57*pow(10,-2)*pow(T,0.66))));
	}
	
	return v_e; //cm/s
}

/*Float_t* vel_h(double n_d, float v_h[n], float E[n], double T)
{
	vector<double> dist(n);
	vector<double> ini(n);
	vector<int> i_h(n);
	for(int i(0); i<n; ++i)
	{
		ini[i]=(double) n_d/n_e*i+n_d/n_e*0.5; //micro
		dist[i]=n_d-ini[i]; //distance à parcourir pour l'électron, micro
		i_h[i]=round(n*ini[i]/n_d); // ATTENTION, ROUND DONC PAS E EXACT
		v_h[i]=(1.31*pow(10,8)*pow(T, -2.2)*E[i_h[i]])/pow((pow((1+(E[i_h[i]]/1.24)*pow(T,1.68)),0.46)*pow(T, 0.17)),(1/0.46*pow(T,0.17)));
	}
	
	return v_h;
}*/

vector<vector<double>> tps(float ini[n_e], double n_d,float E[n], double T) //ini c'est où on injecte // drift time for h and e
{
	vector<double> v_e(n_e);
	vector<double> v_h(n_e);
	vector<double> dist(n_e);
	vector<int> i_e(n_e);
	vector<int> i_h(n_e);
	vector<double> t_e(n_e);
	vector<double> t_h(n_e);
	vector<vector<double>> t(2, vector<double>(n_e));

	
	for(int i(0); i<n_e; ++i)
	{
		ini[i]=(double) n_d/n_e*i+n_d/n_e*0.5; // micrometre
		dist[i]=n_d-ini[i]; //distance à parcourir pour l'électron, micrometre
		i_e[i]=round(n*dist[i]/n_d); // avec d[i]=dist=n_d/n*i
		//cout<<"hello"<<i<<endl;
		//cout<<"dist"<<dist[i]<<endl;
		//cout<<"i_e"<<i_e[i]<<endl;
		i_h[i]=round(n*ini[i]/n_d); // ATTENTION, ROUND DONC PAS E EXACT
		v_e[i]=(1.42*pow(10,9)*pow(T, -2.42)*E[i_e[i]])/pow((pow((1+(E[i_e[i]]/1.01)*pow(T,1.55)),2.57*pow(10,-2))*pow(T, 0.66)),(1/2.57*pow(10,-2)*pow(T,0.66))); //cm /s
		//cout<<"hola"<<i<<endl;
		//cout<<"ini"<<ini[i]<<endl;
		v_h[i]=(1.31*pow(10,8)*pow(T, -2.2)*E[i_h[i]])/pow((pow((1+(E[i_h[i]]/1.24)*pow(T,1.68)),0.46)*pow(T, 0.17)),(1/0.46*pow(T,0.17))); //cm/s
		//v_e[i]=mu_e*E[i_e[i]];
		//v_h[i]=mu_h*E[i_h[i]];
		//t_e[i]=dist[i]/v_e[i];
		//t_h[i]=ini[i]/v_h[i];
		
		t[0][i]=dist[i]*1e4/v_e[i]; //micrometre*1e4/(cm/s)=s
		t[1][i]=ini[i]*1e4/v_h[i];
	}
	cout<<dist[0];
	cout<<ini[0];
	
	//cout<<"t_e="<<t_e<<endl;
	//cout<<"t_h="<<t_h<<endl;
	
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
				{E[i]=V/n_d;}
			else
			{
				/*d.push_back(n_d/n*(n-i));
				E.push_back(V/d[i]);*/
				//d[i]= (double) n_d/n*(n-i);
				d[i] = (double) n_d/n*i; //micrometre
				if (V_d<=V)
				{
					double dif=(V-V_d)/n_d;//volt/micrometre
					E[i]=(double) 1e4*(V_d/n_d+dif+a*d[i]/n_d);// volt/micrometre*1e4=volt/cm
					
				}
				else
				{
					double dif=(V-V_d)/n_d;// volt/micrometre
					E[i]=(double) 1e4*(dif+2.0/n_d*V_d*(1-d[i]/double(n_d)));//I transform it to volt/cm
				}
				/*for (in j(0); j<=nb;++j)
				{
					V_sup[i]=V_d+50-j*nb;
					E_sup[i][j]=(double) V_sup[i]+a*d[i];
				}*/
			}    
	}
}


void fct_I(float I[n_t], float t[n_t],double dt,vector<double> tps, bool cst=0) 
{
	//vector<vector<double>> I_nv(n_t,vector<double>(n_e));
	for(int i(0); i<n_t;++i)
		{
			//t.push_back(s/n_t*i);
			t[i] = (double)1e9* s/n_t*i;//ns*1e9=s
		}
		
	if (cst==0)
	{
		TRandom *rand = new TRandom();
		
		
		vector<double> r(n_t);
		for(int i(0); i<n_t;++i)
		{
			//I.push_back(0);
			r[i] = rand->Gaus(0,1);
			I[i]=r[i];
		}
	}
	else
	{
		//cout<<"tps"<<tps<<endl;
		//cout<<"dt"<<dt<<endl;
		vector<int> num_dt(n_e);
		vector<double> height_I(n_e);
		
				for(int i(0); i<n_e;++i)
				{
					//cout<<"i"<<i;
					//tps[i]*=1e9;
					num_dt[i]=round(tps[i]/(dt*1e9));//nombre de dt
					
					cout<<i<<endl<<"temps  "<<tps[i]<<endl;
					//cout<<"dt  "<<dt<<endl;
					height_I[i]=(double) 1.0/tps[i]; // aire du rectangle =1
					cout<<"height   "<<height_I[i]<<endl;
			
					for(int j(0); j<num_dt[i]; ++j)
					{
						//cout<<"gutenTAg"<<j<<endl;
						I[j]+=height_I[i];
						//cout<<"size I_nv_ligne"<<I_nv[i].size()<<endl<<"size I_nv colonne"<<I_nv.size()<<j<<endl;		
					}
				}
					/*for(int g(num_dt[i]); g<n_t;++g)
					{
						I[g]=0;
						//cout<<"hello"<<g<<endl;
					}*/
					
		}	
		
}
		




void tp4b(int V=300)
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
		cout<<"E_  "<<E_[i]<<endl;
    }
	/*vector<double> E(n);
	vector<double> d(n);
	vector<double> t(n_t);
	vector<double> I(n_t);*/
	/*float v_e[n];
	float v_h[n];*/
	float E[n];
	float d[n];
	float t[n_t];
	float I_e[n_t];
	float I_h[n_t];
	float I_tot[n_t];
	float ini[n_e]; //vecteur d'elec pour pouvoir les placer dans le detecteur
	
	for(int i(0); i<n_t; ++i)
	{
		t[i]=0;
		I_h[i]=0;
		I_e[i]=0;
		I_tot[i]=0;
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

	fct_I(I_e, t, dt, temps[0], 1);
	fct_I(I_h, t, dt, temps[1], 1);
	
	for(int i(0); i<n_t;++i)
	{
		I_tot[i]=I_h[i]+I_e[i];
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
	/*TMultiGraph *mg = new TMultiGraph();
	mg->SetTitle("Signals");*/
	TH2F *hpx1 = new TH2F("hpx1", "I", 20, 0, s, 100, -1, 100);
	hpx1->Draw();
	hpx1->GetYaxis()->SetTitle("I");
	hpx1->GetYaxis()->SetLabelSize(0.05);
	hpx1->GetXaxis()->SetLabelSize(0.05);
	hpx1->GetYaxis()->SetTitleSize(0.05);
	hpx1->GetXaxis()->SetTitleSize(0.05);
	hpx1->GetXaxis()->SetNoExponent();
	hpx1->SetTitle("I, ns vs I");
	hpx1->GetYaxis()->CenterTitle();
	hpx1->GetXaxis()->SetTitle("ns");
	hpx1->GetXaxis()->CenterTitle();
	

	TGraph *gr_I_e = new TGraph (n_t, t, I_e);
	gr_I_e ->SetMarkerStyle(20);
	gr_I_e ->SetMarkerColor(2);
	gr_I_e ->SetMarkerSize(1.0);
	gr_I_e ->SetLineWidth(2);
	gr_I_e ->SetLineColor(2);
	
	
	
	TGraph *gr_I_h = new TGraph (n_t, t, I_h);
	gr_I_h ->SetMarkerStyle(20);
	gr_I_h ->SetMarkerColor(3);
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
	


	Float_t* ve=vel_e(E_, 300);
	for(int j(0); j<n; ++j)
	{
		cout<<"ve "<<j<< ":  "<<ve[j]<<endl;
	}
	TCanvas *canvas2 = new TCanvas("canvas", "Electron velocity vs E", 200, 10, 1000, 650);
    TH2F *hpx2 = new TH2F("hpx2", "Electric", 100, 1e5,n, 1e6, 1e5, 1e6);
	
	hpx2->Draw();
	hpx2->GetYaxis()->SetTitle("vel");
	hpx2->GetYaxis()->SetLabelSize(0.05);
	hpx2->GetXaxis()->SetLabelSize(0.05);
	hpx2->GetYaxis()->SetTitleSize(0.05);
	hpx2->GetXaxis()->SetTitleSize(0.05);
	hpx2->GetXaxis()->SetNoExponent();
	hpx2->SetTitle("vel, v vs E");
	hpx2->GetYaxis()->CenterTitle();
	hpx2->GetXaxis()->SetTitle("E");
	hpx2->GetXaxis()->CenterTitle();

	
	TGraph *graph_Vel = new TGraph(n, E_, ve);
    graph_Vel->SetTitle("Electron velocity vs E");
    graph_Vel->GetXaxis()->SetTitle("E");
    graph_Vel->GetYaxis()->SetTitle("Velocity");
    graph_Vel ->SetMarkerStyle(20);
	graph_Vel ->SetMarkerColor(2);
	graph_Vel ->SetMarkerSize(1.0);
	graph_Vel ->SetLineWidth(2);
	graph_Vel ->SetLineColor(2);
	
	graph_Vel->Draw("L");
	
	
	//mg->Add(gr_I_h);
	//mg->Add(gr_I_e);	
	//mg->Draw("APLC");
	
	
}