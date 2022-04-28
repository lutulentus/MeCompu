//
//  main.cpp
//  SchroedingerTD
//
//  da completare
//

#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>

using namespace std;

void solve_tridiagonal(int n, complex<double> *d, complex<double> *u, complex<double> *l,complex<double> *b, complex<double> *a ){
  //b sono i dati e a l'incognita; b, a, d,u,l  
    complex<double>* alfa = new complex<double>[n];
    complex<double>* beta = new complex<double>[n];
    
    alfa[0]=-d[0]/u[0];//dj=Mj,i uj=Mj,i+1
    beta[0]=b[0]/u[0];//b=bi
    
    for(int i=1;i < n; i++){
        alfa[i]=(-l[i-1]/(u[i]*alfa[i-1])-d[i]/u[i]);//lj=Mj,i-1
        beta[i]=b[i]/u[i]+l[i-1]*beta[i-1]/(u[i]*alfa[i-1]);//
                 
    }
    
    
    a[n-1]=(b[n-1]/l[n-2]+beta[n-2]/alfa[n-2])/(1./alfa[n-2]+d[n-1]/l[n-2]);
    
    for(int i=n-1;i>0;i--){
        a[i-1]=a[i]/alfa[i-1]-beta[i-1]/alfa[i-1];
    }
    
    
    free(alfa);
    free(beta);
    
    
}


int main(int argc, const char * argv[])
{
    double L,x0,q,sigma;//pacchetto
    double a,b, V0;//potenziale
    double dt;
    int Nx,Nsteps,Nprint;//step fisici psi(Nx[0])=psi(Nx[Nx])=0
    complex<double> pic (0.,1.0);
    
    cout << "Lunghezza L\n";
    cin >> L;
    cout << "Posizione x0\n";
    cin >> x0;
    cout << "Momento q\n";
    cin >> q;
    cout << "Sigma \n";
    cin >> sigma;
    cout << "Limite potenziale a\n";
    cin >> a;
    cout << "Limite potenziale b\n";
    cin >> b;
    cout << "Valore potenziale V0\n";
    cin >> V0;
    cout << "Numero punti asse x Nx\n";
    cin >> Nx;
    cout << "Time step dt \n";
    cin >> dt;
    cout << "Numero time steps\n";
    cin >> Nsteps;
    cout << "Scrivi ogni numero passi\n";
    cin >> Nprint;

    
    double h=L/(Nx-1);//se ho 3 gridstep ho solo 2 intervalli
    int Nmat=Nx-2;//formula 221 dispense
    double norm;
    
    complex<double> *psi0 = new complex<double>[Nmat];
    complex<double> *psi1 = new complex<double>[Nmat];
    complex<double> *d = new complex<double>[Nmat];
    complex<double> *u = new complex<double>[Nmat];
    complex<double> *l = new complex<double>[Nmat];
    complex<double> *f = new complex<double>[Nmat];
//SETTA PSI0 e NORMALIZZA
    
//NMAT E' IL NUMERO DI GRID STEPS INTERNI DOVE LA FUNZIONE D'ONDA PUO'
//ESSERE DIVERSA DA ZERO-> per condiz su D(H) ψ(a)=ψ(b)=0
//Nx E' IL NUMERO TOTALE DI GRID STEPS COMPRESI GLI ESTREMI

    norm=0.;//normalizza psi
    //int psi² dx = (sum psi[i]²)·dx
    //ciclo for per sum psi²
    for(int i=1;i<Nx-1;i++){//c.cont. => psi[0]=psi[N]=0
        double x=h*i;
    //formula lab7.pdf: psi(x,t0)=e^(iqx)·e^(-(x-x0)²/2σ)    
        psi0[i-1]=exp(pic*q*x)*exp(-pow(x-x0,2)/(2*pow(sigma,2)));
        norm+=pow(abs(psi0[i-1]),2);
        
    }

    norm=norm*L/Nx;//L/Nx=dx e cosi ho la norma²
    norm=sqrt(norm);
    cout << "Norma" <<norm <<'\n';
    for(int i=0;i<Nmat;i++){
        psi0[i]=psi0[i]/norm;
    }
    
    //SETTA MATRICE M
    
    for(int i=0;i<Nmat;i++){
        u[i]=1.;
        l[i]=1.;
    }
    
   
// QUI CODARE LA PARTE DIAGONALE DELLA MATRICE
// STARE ATTENTI AGLI ESTREMI: nel punto precedente all'inizio/dopo la fine psi=0
    //el diagonali i=j di (221)

 
    ofstream fileg;
    string nomeg;
    nomeg = (string) "psi_tutta.dat";
    fileg.open(nomeg,ios::out);
    fileg.precision(10);

    
    //LOOP SU PASSI
    for(int n=0;n<Nsteps;n++){
    

// QUI CODARE L'UPDATE DELLA FUNZIONE D'ONDA   
//CHE VIENE MESSA IN psi1
        //chiamare solve tridiagonal
        if(n%Nprint==0){
            for(int i=0;i<Nmat;i++){
                double x=h*(i+1);
                fileg << x << "  " << pow(abs(psi0[i]),2) << '\n';
               
            }

            fileg << '\n';

        }
        
        //a ogni passo calcola norma; ci aspettiamo che si conservi con Crank-Nicolson
        norm=0.;
        for(int i=0;i<Nmat;i++){
            psi0[i]=psi1[i];
            norm+=pow(abs(psi0[i]),2);
        }
        norm=norm*L/Nx;
        cout << "Passo " << n << "Norma " << norm << '\n';

    }
    
    fileg.close();
    free(psi0);
    free(psi1);
    free(d);
    free(u);
    free(l);
    free(f);
    
    return 0;
}

