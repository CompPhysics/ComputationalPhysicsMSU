#include "solver.h"
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <ctype.h>
using namespace std;

/**
 * Default initializaiton is euler
 */
solver::solver(){
    this->method="Euler";
    this->origin="sun";
    this->num_planets=this->planets.size();
}

/**
 * Creates a solver with the specified method
 * @param method ("Euler" or "VV")
 */
solver::solver(string method, string origin){
    this->method=method;
    this->origin=origin;
    this->num_planets=this->planets.size();
}

/**
 * Adds planet to solver array
 * SUN IS ASSUMED TO BE PLANET[0]
 * @param p (planet)
 */
void solver::add_planet(planet p){
    if(this->method=="Euler" && this->planets.size()==2){
        cout<<"Error, Euler method only written for binary system"<<endl;
    }
    else{
        this->planets.push_back(p);
        this->num_planets=this->planets.size();
        //create ofstream map
        string filename=this->method;
        filename+="_";
        filename+=p.name;
        string file_x=filename+"_x.txt";
        string file_y=filename+"_y.txt";
        string file_z=filename+"_z.txt";
        string file_xy=filename+"_xy.txt";
        string file_dist=filename+"_dist.txt";
        vector<shared_ptr<ofstream>>v;
        v.push_back(make_shared<ofstream>(file_x));
        v.push_back(make_shared<ofstream>(file_y));
        v.push_back(make_shared<ofstream>(file_z));
        v.push_back(make_shared<ofstream>(file_xy));
        v.push_back(make_shared<ofstream>(file_dist));
        this->output[p.name]=v;
    }
    
}

/**
 * Calculates acceleration/r for the given planet p
 * @param p planet to calculate acceleration of
 * @return acceleration/r
 */
double solver::calc_acc(planet p){
    double pi= acos(-1);
    double min;
    int num_planets=this->num_planets;
    
    double acc=1;
    if(p.name!="sun"){
        acc=1/pow(p.distance(this->planets[0]),3);
    }
    for(int p1=1;p1<num_planets;p1++){
        //don't add force due to itself
        if(this->planets[p1].name==p.name){
            continue;
        }
        //acc+= (Mp1/Ms)*(1/rpp1)^3 (rpp1 is distance between p and p1)
        acc+=this->planets[p1].get_mass()/(this->planets[0].get_mass()*pow(p.distance(this->planets[p1]),3));
    }
    
    acc=4*pi*pi*acc;
    return acc;
}

/**
 * 
 * @return angular momentum of system
 */
double solver::calc_ang_momentum(){
    double p[3];
    double L[3];
    double v_vec[3];
    double mass;
    L[0]=0;
    L[1]=0;
    L[2]=0;
    for(int i=0;i<this->num_planets;i++){
        mass=this->planets[i].get_mass();
        v_vec[0]=this->planets[i].get_velocity()[0];
        p[0]=mass*v_vec[0];
        v_vec[1]=this->planets[i].get_velocity()[1];
        p[1]=mass*v_vec[1];
        v_vec[2]=this->planets[i].get_velocity()[2];
        p[2]=mass*v_vec[2];
        L[0]+=(this->planets[i].x[1]*p[2]-this->planets[i].x[2]*p[1]);
        L[1]-=(this->planets[i].x[0]*p[2]-this->planets[i].x[2]*p[0]);
        L[2]+=(this->planets[i].x[0]*p[1]-this->planets[i].x[1]*p[0]);
    }    
    double ang_mom=sqrt(L[0]*L[0]+L[1]*L[1]+L[2]*L[2]);
    return ang_mom;
}

/**
 * 
 * @return total energy of system
 */
double solver::calc_energy(){
    int num_planets = this->num_planets;
    double kinetic=0;
    double mass;
    double velocity2; //v^2
    double v_vec[3];
    
    //calculate kinetic energy
    for(int p=0;p<num_planets;p++){
        mass=this->planets[p].get_mass();
        v_vec[0]=this->planets[p].get_velocity()[0];
        v_vec[1]=this->planets[p].get_velocity()[1];
        v_vec[2]=this->planets[p].get_velocity()[2];
        velocity2=v_vec[0]*v_vec[0]+v_vec[1]*v_vec[1]+v_vec[2]*v_vec[2];
        kinetic+=0.5*mass*velocity2;
    }
    
    //calculate potential energy
    double pi= acos(-1);
    double potential=0;
    for(int p=0;p<num_planets-1;p++){
        for(int p2=p+1;p2<num_planets;p2++){
            //U=G*Ms*(m1*m2)/(r12*Ms)
            potential+=this->planets[p].get_mass()*this->planets[p2].get_mass()/
                    (this->planets[0].get_mass()*this->planets[p].distance(this->planets[p2]));
        }
    }
    potential=potential*4*pi*pi;
    
    return (kinetic+potential);
}

/**
 * Writes x,y,z,xy,and dist of planet p at time t to files
 * @param t time
 * @param p planet to write
 */
void solver::write(double t, planet p){
    auto& out_x=this->output[p.name][0];
    *out_x<<t<<"   "<<p.x[0]<<endl;
    auto& out_y=this->output[p.name][1];
    *out_y<<t<<"   "<<p.x[1]<<endl;
    auto& out_z=this->output[p.name][2];
    *out_z<<t<<"   "<<p.x[2]<<endl;
    auto& out_xy=this->output[p.name][3];
    *out_xy<<p.x[0]<<"   "<<p.x[1]<<endl;
    auto& dist=this->output[p.name][4];
    double distance=sqrt(p.x[0]*p.x[0]+p.x[1]*p.x[1]+p.x[2]*p.x[2]);
    *dist<<t<<"   "<<distance<<endl;
}

/**
 * Solve differential equation using method assigned to class
 * @param n number of integration points
 * @param h step size
 */
void solver::solve(int n, double h){
    double old_x[3];
    double old_v[3];
    double old_dist;
    double E_i=calc_energy();
    double L_i=calc_ang_momentum();
    int min;
    if(this->origin=="sun"){
        min=1;
    }
    else if(this->origin=="cm"){
        min=0;
    }
    //loop through all planets in system and solve for each
     //Euler only used for 2 planets
     if (this->method=="Euler"){
            double pi=acos(-1);
            double p=1; //only need to calculate for earth which is planets[1]]
            //write initial position to files
            this->write(0,this->planets[1]);
            for (int i=1; i<=n;i++){
                //get distance from last iteration
                old_dist=this->planets[p].distance(this->planets[0]);
                //get values from last iteration and calculate new vals
                for(int dim=0;dim<3;dim++){
                    old_x[dim]=this->planets[p].x[dim];
                    old_v[dim]=this->planets[p].v[dim];
                    this->planets[p].x[dim]=old_x[dim]+h*old_v[dim];
                    this->planets[p].v[dim]=old_v[dim]-h*4*pi*pi*old_x[dim]/pow(old_dist,3);
                }
                //write positions to files
                this->write(i*h,this->planets[1]);
            }  
     }

        if (this->method=="VV"){
            double old_a[3];
            double a[3];
            //write initial positions to files
            for(int p=min;p<num_planets;p++){
                this->write(0,this->planets[p]);
            }
            //iterate over time
            for (int i=1; i<=n;i++){
                //store acc(ti) in planet acc member
                for(int p=min;p<num_planets;p++){
                    this->planets[p].acc=calc_acc(this->planets[p]);
                }
                for(int p=min;p<num_planets;p++){
                    //get values from last iteration and calculate new position
                    for(int dim=0;dim<3;dim++){
                        old_x[dim]=this->planets[p].x[dim];
                        old_v[dim]=this->planets[p].v[dim];
                        old_a[dim]=-this->planets[p].acc*old_x[dim];
                        this->planets[p].x[dim]=old_x[dim]+h*old_v[dim]+h*h*old_a[dim]/2;
                        //calculate updated acceleration
                        //acc=a/r for x component of a need acc*x=a*x/r (same for y,z))
                        a[dim]=-this->planets[p].acc*this->planets[p].x[dim];
                        this->planets[p].v[dim]=old_v[dim]+h*(a[dim]+old_a[dim])/2;
                        this->write(i*h,this->planets[p]);
                    }
                }
            }
        }
    double E_f=calc_energy();
    cout<<"Energy change (%): "<<100*(E_i-E_f)/E_i<<endl;
    double L_f=calc_ang_momentum();
    cout<<"Angular momentum change (%): "<<100*(L_i-L_f)/L_i<<endl;
}

solver::~solver() {
}

