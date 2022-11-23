#define _USE_MATH_DEFINES // for C++
#include <cmath>

#include <iostream>
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace std;

const double t_0 = 0;                       //tiempo de inicio de la simulación
const double r_0 = 6671;                    //radio incial //=r_t
const double v_r_0 = 0;                     // vel radial inicial 
const double theta_0 = 0;                   // angulo inicial
const double v_theta_0 = 1.163552835E-3;    // velocidad angular inicial
const double t_f = 40000;                   //tiempo maximo de la simumlación
const double dt = 0.1;                      // paso de tiempo

const double G = 6.673E-20;                 //constante de la gravedad en km^3
const double m_T = 5.972E24;                //masa de la tierra en kg
const double P_T = 100;                      //potencia del motor en newtons
const double D = 0.75;                      //coheficiente de friccion del cohete 
const double ro_0 = 1.225;                  //dencidad de la atmosfera en la superficie
const double r_T = 6371;                    //densidad del aire en la superficie
const double rendimiento = 1;            //gasto de combustible por segundo a maxima potencia
double const m_0 = 10000;                   //masa inicial del cohete 3 000 000
const double combustible = 500;             //masa del combustible (masa utilizable)

bool comb_disp = true;                      //variable que permite saber si hay conbustible para realizar una maniobra
double potencia (int i, double t) {         //retorna el porcentaje de potencia del motor
    if (comb_disp == true){
        if(5400<=t && t<=5700){
            return 1;
        }
        else{
            return 0;
        }
    }
    else {
        return 0;
    }
}

double phi (int i, double t) {              // retorna el rumbo el cohete 
    double theta_0 = 0;
    return (-M_PI/2);//placeholder
}

double a_r (int i, double m, double t, double r,double v_r, double theta, double v_theta) { //aceleracion radial


    return r*v_theta*v_theta - G*m_T/(r*r) + potencia(i, t)*(P_T/m)*cos(theta - phi(i, t)) - (D/m)*ro_0*exp(-(r-r_T)/10.4)*v_r*v_r;
}

double a_theta (int i, double m, double t, double r,double v_r, double theta, double v_theta) {
    return -2*(v_r*v_theta)/r + potencia(i, t)*(P_T/(m*r))*sin(theta - phi(i, t)) - (D/(m*r))*ro_0*exp(-(r-r_T)/10.4)*r*r*v_theta*v_theta;
}


void runge_kutta_4 (double t_0, double r_0, double v_r0, double theta_0, double v_theta_0, double t_f, double dt, string nombre_archivo) {
    ofstream datos;
    datos.open (nombre_archivo);
    
    double m = m_0;

    // L depende de velocidad r
    // M depende de aceleracion r
    // H depende de vel angular
    // Q depende de accel angua+lar 
    // el sufijo _n signfica que representa el paso anterior de tiempo 
    // el sufijo _np'''''''''''''''''''''   el paso que estamos calculando

    double r_n = r_0;           double v_r_n = v_r0;                double theta_n = theta_0;           double v_theta_n = v_theta_0;
    double r_np;                double v_r_np;                      double theta_np;                    double v_theta_np;
    double L1;                  double M1;                          double H1;                          double Q1;
    double L2;                  double M2;                          double H2;                          double Q2;
    double L3;                  double M3;                          double H3;                          double Q3;
    double L4;                  double M4;                          double H4;                          double Q4;
    double suma_r;              double suma_v_r;                    double suma_theta;                  double suma_v_theta;
    int i = 0; //el paso en el que vamos 
     
    for (double t = t_0 ; t <= t_f; t = t + dt) {
        datos<<t<<","<<r_n<<","<<theta_n<<","<<m<<"\n";
        L1 = v_r_n;                                                 M1 = a_r (i, m, t, r_n, v_r_n, theta_n, v_theta_n);
        H1 = v_theta_n;                                             Q1 = a_theta (i, m, t, r_n, v_r_n, theta_n, v_theta_n);
        L2 = v_r_n + (0.5*dt*M1);                                   M2 = a_r (i, m, t + 0.5*dt, r_n + (0.5*dt*L1), v_r_n + (0.5*dt*M1), theta_n + (0.5*dt*H1), v_theta_n + (0.5*dt*Q1));
        H2 = v_theta_n + (0.5*dt*Q1);                               Q2 = a_theta (i, m, t + 0.5*dt, r_n + (0.5*dt*L1), v_r_n + (0.5*dt*M1), theta_n + (0.5*dt*H1), v_theta_n + (0.5*dt*Q1));
        L3 = v_r_n + (0.5*dt*M2);                                   M3 = a_r (i, m, t + 0.5*dt, r_n + (0.5*dt*L2), v_r_n + (0.5*dt*M2), theta_n + (0.5*dt*H2), v_theta_n + (0.5*dt*Q2));
        H3 = v_theta_n + (0.5*dt*Q2);                               Q3 = a_theta (i, m, t + 0.5*dt, r_n + (0.5*dt*L2), v_r_n + (0.5*dt*M2), theta_n + (0.5*dt*H2), v_theta_n + (0.5*dt*Q2));
        L4 = v_r_n + dt*M3;                                         M4 = a_r (i, m, t + dt, r_n + dt*L3, v_r_n + dt*M3, theta_n + dt*H3, v_theta_n + dt*Q3);
        H4 = v_theta_n + dt*Q3;                                     Q4 = a_theta (i, m, t + dt, r_n + dt*L3, v_r_n + dt*M3, theta_n + dt*H3, v_theta_n + dt*Q3);
        r_np = r_n + (dt*(L1 + (2*L2) + (2*L3) + L4))/6;            v_r_np = v_r_n + (dt*(M1 + (2*M2) + (2*M3) + M4))/6;
        theta_np = theta_n + (dt*(H1 + (2*H2) + (2*H3) + H4))/6;    v_theta_np = v_theta_n + (dt*(Q1 + (2*Q2) + (2*Q3) + Q4))/6;
        m = m - rendimiento*potencia(i,t);
        r_n = r_np;                                                                                 v_r_n = v_r_np;
        theta_n = theta_np;                                                                         v_theta_n = v_theta_np;
       
        i = i + 1;
        if (r_n <= r_T) {
            break;
        }
        if (m <= (m_0 - combustible)) {
            comb_disp = false;
        }
    }
    datos.close();
}

int main() {
    runge_kutta_4 (t_0, r_0, v_r_0, theta_0, v_theta_0, t_f, dt,"datos_runge_kutta.txt"); //ejecutar la sim
    return 0;
}