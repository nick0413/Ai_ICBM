#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace std;

const double G = 6.673E-20;                 //constante de la gravedad en km^3
const double m_T = 5.972E24;                //masa de la tierra en kg
const double P_T = 260;                      //potencia del motor en KiloNewtons
const double D = 0.75;                      //coheficiente de friccion del cohete
const double ro_0 = 1.225;                  //densidad de la atmosfera en la superficie kg/m^3
const double r_T = 6371;                    //radio de la tierra en km
const double rendimiento = 83.49;           //gasto de combustible en kg  por segundo a maxima potencia
const double m_0 = 18000;                   //masa inicial del cohete en kg
const double combustible = 17000;             //masa del combustible en kg(masa utilizable)


const int num_cer = 50;                      //numero de cerebros que seran simulados         
const int num_dat = 2000;                  //numero de datos 
double daticos[2*num_cer][num_dat];            //arreglo para almacenar los datos de cada cerebro
double condiciones[1][2];                   //arreglo para guardar las condiciones de simulacion
int i_cerebro;                              //indice del cerebro que estamos simulando

ofstream fitness;                           //declara la variable para el ouput de datos de fitness


void leer_archivo_datos (string nombre_archivo){//extrae los datos de vuelo de un archivo txt y los organiza en un arreglo
    ifstream archivo_datos;
    archivo_datos.open(nombre_archivo);
    
    for (int i=0; i<2*num_cer; i++){
        for (int j = 0; j<num_dat; j++){
            archivo_datos>> daticos[i][j];
        }
    }
    archivo_datos.close();
}

void leer_archivo_condiciones (string nombre_archivo){//extrae los datos de un archivo txt y los guarda en un arreglo
    ifstream archivo_condiciones;
    archivo_condiciones.open(nombre_archivo);

    for (int j = 0; j < 2; j++){
        archivo_condiciones >> condiciones[0][j];
    }
    
    archivo_condiciones.close();

}
const double t_0 = 0;                       //tiempo de inicio de la simulación
const double r_0 = r_T;                    //radio incial //=r_t
const double v_r_0 = 0;                     // vel radial inicial 
const double theta_0 = condiciones[0][0];    // angulo inicial
double theta_objetivo = condiciones[0][1]; //angulo objetivo de los cerebros
double heading;                             //variable orientacion del cohete
const double v_theta_0 = 0;                 // velocidad angular inicial
const double t_f = 20000;                   //tiempo maximo de la simumlación
const double dt = 0.1;                      // paso de tiempo

bool comb_disp = true;                      //variable que permite saber si hay conbustible para realizar una maniobra


double potencia (int i, double t) {         //retorna el porcentaje de potencia del motor
    //if (comb_disp==true){
        //if (0<=t && t<=1000){
            //return 1;
        //}
        //else{
            //return 0;                     //BORRAME
        //}
    //}
    //else {
        //return 0;
    //} 
    if (comb_disp == true){                 //retorna el valor del impulso en cada instante de tiempo
        return 10*abs(daticos[2*i_cerebro][i]); 
    }
    else{
        return 0;
    }

}

double phi (int i, double t) {              // retorna el rumbo el cohete
    //if (0<=t && t <= 60){
        //return theta_0;
    //}
    //else{                     //BORRAME
        //return theta_0 + M_PI/6;
    //}

    heading = heading +0.05*daticos[2*i_cerebro+1][i]/8;//retornar el valor de la direccion en cada instante
    //cout <<heading<<" "<<data[2*i_cerebro+1][i]<<"\n";
    return heading;
}

double a_r (int i, double m, double t, double r,double v_r, double theta, double v_theta) { //aceleracion radial
    return r*v_theta*v_theta - G*m_T/(r*r) + potencia(i, t)*(P_T/m)*cos(theta + phi(i, t)) - (D/m)*ro_0*exp(-(r-r_T)/10.4)*v_r*v_r;
}

double a_theta (int i, double m, double t, double r,double v_r, double theta, double v_theta) {//acerelacion angular
    return -2*(v_r*v_theta)/r + potencia(i, t)*(P_T/(m*r))*sin(theta + phi(i, t)) - (D/(m*r))*ro_0*exp(-(r-r_T)/10.4)*r*r*v_theta*v_theta;
}


void runge_kutta_4 (double t_0, double r_0, double v_r0, double theta_0, double v_theta_0, double t_f, double dt, string nombre_archivo) {
    ofstream datos;
    datos.open (nombre_archivo);
    
    double m = m_0;         //resetea las siguientes variables para cada simulacion
    heading = theta_0;
    comb_disp = true;

    double thrust_prom=0;
    double max_height=0;
    double prom_thetav=0;
    // L depende de velocidad r
    // M depende de aceleracion r
    // H depende de vel angular
    // Q depende de accel angular 
    // el sufijo _n representa el paso anterior de tiempo 
    // el sufijo _np'''''''''' el paso que estamos calculando

    double r_n = r_0;           double v_r_n = v_r0;                double theta_n = theta_0;           double v_theta_n = v_theta_0;
    double r_np;                double v_r_np;                      double theta_np;                    double v_theta_np;
    double L1;                  double M1;                          double H1;                          double Q1;
    double L2;                  double M2;                          double H2;                          double Q2;
    double L3;                  double M3;                          double H3;                          double Q3;
    double L4;                  double M4;                          double H4;                          double Q4;
    double suma_r;              double suma_v_r;                    double suma_theta;                  double suma_v_theta;
    int i = 0; //el paso en el que vamos 
     
    for (double t = t_0 ; t <= t_f; t = t + dt) {
        if (r_n < r_T) {
            //fitness<< ((m_0-m)*(m_0-m))*1E-6+ (theta_objetivo-theta_n)*(theta_objetivo-theta_n)*r_n*r_n<< "\n";
            thrust_prom=thrust_prom/(t/dt);
            max_height=max_height/r_T;
            
            //cout<<thrust_prom<<"\n";
            fitness<< setprecision(10) <<(theta_objetivo-theta_n)<<"\n";
            //fitness<< setprecision(10) << pow(prom_thetav+1,40)*abs(pow(thrust_prom,10)*pow(max_height,20)) << "\n";
            //fitness<< setprecision(10) << pow(prom_thetav+1,40)*prom_thetav << "\n";
            //cout<<"fitness_writen"<<"\n";
            //if(i_cerebro==18){fitness<<(m_0)<<"-----"<<(theta_objetivo)<<"\n";}
            break;
        }
        datos<<t<<","<<r_n<<","<<theta_n<<","<<daticos[2*i_cerebro+1][i]<<","<<potencia(i,t)<<"\n"; //DELETEME
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

        m = m - rendimiento*potencia(i,t)*dt;        //actualiza el combustible consumido en el interbalo dt
        r_n = r_np;                                                                                 v_r_n = v_r_np;         //actualiza las 4 variables para el siguiente paso
        theta_n = theta_np;                                                                         v_theta_n = v_theta_np;

        if(r_n>max_height){max_height=r_n;};

        i = i + 1;
        thrust_prom+=potencia(i,t)*dt;
        prom_thetav+=v_theta_np;
        if (m <= (m_0 - combustible)) {             //actualiza el valor de la variable que determina si hay combustible 
            comb_disp = false;
        }
        
    }
    datos.close();
}

int main() {


    leer_archivo_datos("datos_gen1.txt");

    leer_archivo_condiciones("condiciones_gen1.txt");
    fitness.open("fitness.txt");

    //cout<<condiciones[0][0]<<" "<<condiciones[0][1]<<"\n";
    theta_objetivo = condiciones[0][1];
    cout<<theta_objetivo<<"\n";
    for (i_cerebro = 0; i_cerebro < num_cer; i_cerebro ++){
        //resetear el valor de heading antes de la siguiente simulacion 
        runge_kutta_4 (t_0, r_0+10, v_r_0, theta_0, v_theta_0, t_f, dt,"datos_runge_kutta"+to_string(i_cerebro)+".txt"); //ejecutar la sim
    }
    fitness.close();

    return 0;


}