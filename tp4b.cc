#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <typeinfo>
using namespace std;


	
const int n=1000; //how many pts for d and E
double n_d=100; //thickness, micro
double dx=(double) n_d/n; //dx of E and d //micro
//int V=300; //external
int s=100; //how many ns
int n_t=1000; //how many pts for s and I
double dt= (double) s/n_t; //ns
int V_d=300; //depletion voltage
int n_e=10; //nb electrons
double T=298.16; //K
double cutoff=2; //GHz

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
static Float_t vel_h_at_position_in_field(float value_of_field, double temp);
Float_t* vel_e( float E_[n], double T)
{
	static Float_t v_e[n];
	for(int i(0); i<n; ++i)
	{	double num=1.42e9*E_[i]*pow(T,-2.42);
		double den = pow(1.0+pow(E_[i]/1.01/pow(T,1.55),(2.57e-2*pow(T,0.66))),(1/(2.57e-2*pow(T,0.66))));
		v_e[i]=num/den;
		//v_e[i]=(1.42*pow(10,9)*pow(T, -2.42)*E_[i])/pow((pow((1+(E_[i]/1.01)*pow(T,1.55)),2.57*pow(10,-2))*pow(T, 0.66)),(1/(2.57*pow(10,-2)*pow(T,0.66))));
	}
	
	return v_e; //cm/s
}

Float_t* vel_h(float E_[n], double T)
{
	
	double pot_h=0.46*pow(T, 0.17);
	static Float_t v_h[n];

	transform(E_, E_+n, v_h, [T](float E){return vel_h_at_position_in_field(E,T);});
	return v_h;
}

Float_t vel_h_at_position_in_field(float value_of_field, double temp)
{
	double pot_h=0.46*pow(temp, 0.17);
	double num_=1.31e8*value_of_field*pow(temp,-2.2);
	double den_ = pow(1.0+pow(value_of_field/1.24/pow(temp,1.68),pot_h),(1/pot_h));
		
	return num_/den_;
}

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
	double pot_h=0.46*pow(T, 0.17);
	
	for(int i(0); i<n_e; ++i)
	{

		
		ini[i]=(double) n_d/n_e*i+n_d/n_e*0.5; // micrometre
		dist[i]=n_d-ini[i]; //distance à parcourir pour l'électron, micrometre
		i_e[i]=floor(abs(n*dist[i]/n_d)); // avec d[i]=dist=n_d/n*i
		//cout<<"hello"<<i<<endl;
		//cout<<"dist"<<dist[i]<<endl;
		//cout<<"i_e"<<i_e[i]<<endl;
		i_h[i]=floor(abs(n*ini[i]/n_d)); // ATTENTION, ROUND DONC PAS E EXACT
		//the function E returns us actually a voltage, which means I have to derive it by the thickness (in micro) and then
		//multiply by 1e4
		//E[i_h[i]]=E[i_h[i]]/n_d*1e4;
		//E[i_e[i]]=E[i_e[i]]/n_d*1e4;
		double num=1.42e9*E[i_e[i]]/n_d*1e4*pow(T,-2.42);
		double den = pow(1.0+pow(E[i_e[i]]/n_d*1e4/1.01/pow(T,1.55),(2.57e-2*pow(T,0.66))),(1/(2.57e-2*pow(T,0.66))));
		v_e[i]=num/den;
		//v_e[i]=(1.42*pow(10,9)*pow(T, -2.42)*E[i_e[i]]/n_d*1e4)/pow((1+pow(E[i_e[i]]/n_d*1e4/(1.01*pow(T,1.55)),2.57*pow(10,-2))*pow(T, 0.66)),(1/(2.57*pow(10,-2)*pow(T,0.66)))); //cm /s
		//cout<<"hola"<<i<<endl;
		//cout<<"ini"<<ini[i]<<endl;
		double num_=1.31e8*E[i_h[i]]/n_d*1e4*pow(T,-2.2);
		double den_ = pow(1.0+pow(E[i_h[i]]/n_d*1e4/1.24/pow(T,1.68),pot_h),(1/pot_h));
		
		v_h[i]=num_/den_;
		//v_h[i]=(1.31*pow(10,8)*pow(T, -2.2)*E[i_h[i]]/n_d*1e4)/pow((1+pow(E[i_h[i]]/n_d*1e4/(1.24*pow(T,1.68)),0.46)*pow(T, 0.17)),(1/(0.46*pow(T,0.17)))); //cm/s
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
		
	
		
		//This is actually the Voltage as a function of the thickness, so the electric field would be the V/d
		double a=(double)-V_d/n_d; //pente du champ E, V/micrometre
		for(int i(0); i<n;++i)
		{
			if (cst==0)
				{E[i]=V;}
			else
			{
				/*d.push_back(n_d/n*(n-i));
				E.push_back(V/d[i]);*/
				//d[i]= (double) n_d/n*(n-i);
				d[i] = (double) abs(n_d/n*i); //micrometre
				//if (V_d<=V)
				//{
					//double dif=(V-V_d)/n_d;//volt/micrometre
					//E[i]=(double) 1e-4*(V_d/n_d+dif+a*d[i]/n_d);// volt/micrometre*1e4=volt/cm
					//E[i]=(double) (a*d[i]+V)/n_d/1e4; // /n_d/1e4//en V/micro/1e4=V/cm
					E[i]=(double) a*d[i]+V;
					//cout<<"E1   "<<E[i]<<endl;
					
				//}
				/*else
				{
					double dif=(V-V_d)/n_d;// volt/micrometre
					E[i]=(double) abs(1e-4*(dif+2.0/n_d*V_d*(1-d[i]/double(n_d))));//I transform it to volt/cm
					cout<<"E2   "<<E[i]<<endl;
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
					//tps in seconds, dt in nanoseconds
					//num_dt[i]=floor(abs(tps[i]/1e9/(dt)));//nombre de dt (arbitrarily put *1e1)
					num_dt[i]=floor(abs(1e1*tps[i]/(dt)));
					//cout<<i<<endl<<"temps  "<<tps[i]<<endl;
					//cout<<"num_dt  "<<num_dt[i]<<endl;
					height_I[i]=(double) 1.0/tps[i]; // aire du rectangle =1
					//cout<<"height   "<<height_I[i]<<endl;
			
					for(int j(0); j<=num_dt[i]; ++j)
					{

						int idx = j;
						if (idx >= n_t) {
							idx = n_t - 1;
						}
						I[idx] += height_I[i];
						//cout<<"gutenTAg"<<j<<endl;
						//cout<<"I   "<<I[j]<<endl;
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
void I_to_frequence(float I[n_t], float freq[n_t]) //∫f(x)e(-i2 pi k x)dx=f(k), x runs from 0 to n_t
{
	//normalising factors??
	float normalization=1.0/n_t;
   for(int k(0); k<n_t; ++k) //fill freq vector
	{	
		freq[k] = 0;
		for (int i = 0; i < n_t; ++i) //integral
		{
			float Re = cos(2*M_PI*k*i/n_t);
			//cout<<"Re  "<<Re<<endl;
			float Im = -sin(2*M_PI*k*i/n_t);
			//cout<<"Im  "<<Im<<endl;
			//cout<<"I_e"<<I[i]<<endl;
			freq[k] += I[i] * (Re + Im);
		}
		cout<<"unfiltered, freq  "<<freq[k]<<endl;
		//freq[k] *= normalization; // added normalization factor
	}
	/*double op=tan(deg*M_PI/180)*1;//tan takes in radians so transform deg to rad
	int nb_dt_0=ceil(abs(1/(dt*(cutoff+op))))+1; //when filter=0
	int nb_dt=ceil(abs(1/(dt*cutoff)));//when not 1 anymore*/
}


/*void filter_(double cutoff,Float_t filter[n_t],int f[n_t],float t[n_t], double deg) //cutoff in Ghz, deg in degrees
{
	//static Float_t f[n_t];

	for(int i(0); i<n_t; ++i)
	{
		f[n_t-i-1]=floor(abs(1e9*1/(t[i]+1))); //in order for f to increase t is in seconds, we want f in GHz, so *1e9
		cout<<"f   "<<f[n_t-i-1]<<endl;
	}
	int j=0;
	while ((j<n_t) and (f[j]<=cutoff)){
		filter[f[j]]=1;
		j+=1;
	} //filters nothing until cutoff
	
	int j=0;

	double a=tan(deg);
	//double b=1-a*f[cut];
	b=1+tan(deg)*cutoff;
	//int max=n_t-j;
	int l=0;
	while((j<n_t) and (l<max) and (a*l+b>=0)) {
		filter[j]=-a*l+b;
		j+=1;
		l+=0;
	}
	
	if(a*l+b<=0)
	{
		while(j<n_t) {
			filter[j]=0;
			j+=1;
		}
	}


}*/


void construct_filter_(float filter[n_t], double deg)
{
	//filter=0 when f in (cut+op,inf)
	//filter =a*l+b when f in (cut, cut+op)
	//filter=1 when f in (0,cut)
	double op=tan(deg)*1;
	int nb_dt_0=ceil(abs(1/(dt*(cutoff+op))));
	int nb_dt=ceil(abs(1/(dt*cutoff)));
	
	cout <<"nb_dt  "<<nb_dt<<endl;

	for(int i(0); i<=nb_dt_0; ++i)
	{
		filter[i]=0;
	}
	
	double a=-tan(deg*M_PI/180);
	cout<<"a   "<<a<<endl;
	double b=1+tan(deg*M_PI/180)*cutoff;
	//double b=1;
	cout<<"b   "<<b<<endl;

	//int l=0;
	cout <<"filter(nb_dt_0)  "<<a*1.0/((nb_dt_0+1)*dt)+b<<endl;
	for(int i(nb_dt_0+1); i<=nb_dt; ++i)
	{
		//cout<<"i  "<<i<<endl;
		//cout<<"filter    "<<(a*l+b)<<endl;
		cout<<"filter    "<<(a*1.0/(i*dt)+b)<<endl;
		filter[i]=(a*1/(i*dt)+b);
	}
	for(int i(nb_dt+1); i<n_t; ++i)
	{
		filter[i]=1;
	}
}

void apply_filter(float freq[n_t], float filter[n_t])
{
	for(int i(0); i<n_t; ++i)
	{
		freq[i]=freq[i]*filter[i];
		cout<<"filtered in freq domain   "<<freq[i]<<endl;
	}
}

void inverse_fourier(float freq[n_t], float new_I[n_t])
{
	//normalising factors?? 
	float normalization_factor = 1.0 / n_t;
	//i is time and k frequence
   for(int i(0); i<n_t; ++i) //fill new_I vector
	{	
		new_I[i] = 0;
		for (int k = 0; i < n_t; ++k) //integral over frequence
		{
			float Re = cos(2*M_PI*k*i/n_t);
			float Im = sin(2*M_PI*k*i/n_t);
			new_I[i] += freq[k] * (Re + Im);
			//cout<<"filtered in time domain   "<<new_I[i]<<endl;
		}
		//new_I[k] *= normalization_factor;
		cout<<"new_I   "<<new_I[k]<<endl;
		
	}

}

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
		//cout<<"E_  "<<E_[i]<<endl;
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
	Float_t filter[n_t];
	Float_t freq[n_t];
	Float_t new_I[n_t];

	for(int i(0); i<n_t; ++i)
	{
		t[i]=0;
		I_h[i]=0;
		I_e[i]=0;
		I_tot[i]=0;
		filter[i]=0;
		freq[i]=0;
		new_I[i]=0;
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
	
	/*for(int i(0); i<n_t;++i)
	{
		I_tot[i]=I_h[i]+I_e[i];
		cout<<"I_tot   "<<I_tot[i]<<endl;
	}*/
	
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
	
	hpx0->GetYaxis()->SetTitle("U [V]");
	
	
	hpx0->GetYaxis()->SetLabelSize(0.05);
	hpx0->GetXaxis()->SetLabelSize(0.05);
	hpx0->GetYaxis()->SetTitleSize(0.05);
	hpx0->GetXaxis()->SetTitleSize(0.05);
	hpx0->GetXaxis()->SetNoExponent();
	hpx0->SetTitle("Voltage, U vs thickness");
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
	/*TMultiGraph *mg = new TMultiGraph();
	mg->SetTitle("Signals");*/
	TH2F *hpx1 = new TH2F("hpx1", "I", 20, 0, 1e8*s, 100, -1, 1000);
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
	
	//filter_(cutoff, &filter[n_t], &f[n_t], &t[n_t],1.0);
	
	I_to_frequence(I_e, freq);
	construct_filter_(filter, 1.0);
	apply_filter(freq, filter);
	inverse_fourier(freq, new_I);
	TCanvas *canvfilter = new TCanvas("canvfilter", "I", 200, 10, 1000, 650);
	canvfilter->SetGrid();
	/*TMultiGraph *mg = new TMultiGraph();
	mg->SetTitle("Signals");*/
	TH2F *hpxfilter = new TH2F("hpxfilter", "I", 20, 0, 1e8*s, 100, -1, 1000);
	hpxfilter->Draw();
	hpxfilter->GetYaxis()->SetTitle("I");
	hpxfilter->GetYaxis()->SetLabelSize(0.05);
	hpxfilter->GetXaxis()->SetLabelSize(0.05);
	hpxfilter->GetYaxis()->SetTitleSize(0.05);
	hpxfilter->GetXaxis()->SetTitleSize(0.05);
	hpxfilter->GetXaxis()->SetNoExponent();
	hpxfilter->SetTitle("I, Current vs time");
	hpxfilter->GetYaxis()->CenterTitle();
	hpxfilter->GetXaxis()->SetTitle("t [ns]");
	hpxfilter->GetXaxis()->CenterTitle();
	

	TGraph *gr_filter = new TGraph (n_t, t, new_I);
	gr_filter->SetMarkerStyle(20);
	gr_filter ->SetMarkerColor(2);
	gr_filter ->SetMarkerSize(1.0);
	gr_filter ->SetLineWidth(2);
	gr_filter ->SetLineColor(2);
	
	


	Float_t* ve=vel_e(E_, 300);
	/*for(int j(0); j<n; ++j)
	{
		cout<<"ve "<<j<< ":  "<<ve[j]<<endl;
	}*/
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
	
	TGraph *graph_Vel = new TGraph(n, E_, ve);
    graph_Vel->SetTitle("Electron velocity vs E");
    graph_Vel->GetXaxis()->SetTitle("E");
    graph_Vel->GetYaxis()->SetTitle("v [cm/s]");
    graph_Vel ->SetMarkerStyle(20);
	graph_Vel ->SetMarkerColor(2);
	graph_Vel ->SetMarkerSize(1.0);
	graph_Vel ->SetLineWidth(2);
	graph_Vel ->SetLineColor(2);
	
	

	Float_t* vh=vel_h(E_, 300);
	/*TCanvas *canvas3 = new TCanvas("canvas3", "Hole velocity vs E", 200, 10, 1000, 650);
    TH2F *hpx3 = new TH2F("hpx3", "Electric", 100,1e2 ,1e5, n,1e4 , 1e7); //nb binsx, xlow, xup, nbinsy, ydown, yup 
	
	hpx3->Draw();
	hpx3->GetYaxis()->SetTitle("vel");
	hpx3->GetYaxis()->SetLabelSize(0.05);
	hpx3->GetXaxis()->SetLabelSize(0.05);
	hpx3->GetYaxis()->SetTitleSize(0.05);
	hpx3->GetXaxis()->SetTitleSize(0.05);
	hpx3->GetXaxis()->SetNoExponent();
	hpx3->SetTitle("vel_h, v vs E");
	hpx3->GetYaxis()->CenterTitle();
	hpx3->GetXaxis()->SetTitle("E");
	hpx3->GetXaxis()->CenterTitle();
	gPad->SetLogx();
	gPad->SetLogy();*/
	
	TGraph *graph_Velh = new TGraph(n, E_, vh);
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
	
	//mg->Add(gr_I_h);
	//mg->Add(gr_I_e);	
	//mg->Draw("APLC");
}