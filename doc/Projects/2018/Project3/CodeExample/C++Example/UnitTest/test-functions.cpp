#include "catch.hpp"
#include "planet.h"
#include "solver.h"
#include <fstream>
#include <iostream>
using namespace std;


TEST_CASE("Testing Binary Euler"){
    double h = 0.01;
    double n = 500;
    planet earth("earth",0.000003,1.,0.0,0.0,0.0,6.3,0.);
    planet sun ("sun",1.,0.,0.,0.,0.,0.,0.);
    
    solver euler("Euler","sun");
    euler.add_planet(sun);
    euler.add_planet(earth);
    euler.solve(n,h);
    
    //these are files output by code
    ifstream euler_x("Euler_earth_x.txt");
    ifstream euler_y("Euler_earth_y.txt");
    ifstream euler_z("Euler_earth_z.txt");
    ifstream euler_dist("Euler_earth_dist.txt");
    
    //test last value in each file
    double t,x,y,z,dist;
    while(euler_x>>t>>x){};
    REQUIRE(x==Approx(1.71371));
    while(euler_y>>t>>y){};
    REQUIRE(y==Approx(1.12968));
    while(euler_z>>t>>z){};
    REQUIRE(z==Approx(0.00));
    while(euler_dist>>t>>dist){};
    REQUIRE(dist==Approx(2.05256));
    
    euler_x.close();
    euler_y.close();
    euler_z.close();
    euler_dist.close();
}

TEST_CASE("Testing Binary VV"){
    double h = 0.01;
    double n = 500;
    
    planet earth("earth",0.000003,1.,0.0,0.0,0.0,6.3,0.);
    planet sun ("sun",1.,0.,0.,0.,0.,0.,0.);
    
    solver vv("VV","sun");
    vv.add_planet(sun);
    vv.add_planet(earth);
    vv.solve(n,h);
    
    //these are files output by code
    ifstream vv_x("VV_earth_x.txt");
    ifstream vv_y("VV_earth_y.txt");
    ifstream vv_z("VV_earth_z.txt");
    ifstream vv_dist("VV_earth_dist.txt");
    
    //test last value in each file
    double t,x,y,z,dist;
    while(vv_x>>t>>x){};
    REQUIRE(x==Approx(0.962506));
    while(vv_y>>t>>y){};
    REQUIRE(y==Approx(-0.288839));
    while(vv_z>>t>>z){};
    REQUIRE(z==Approx(0.00));
    while(vv_dist>>t>>dist){};
    REQUIRE(dist==Approx(1.00491));
    
    vv_x.close();
    vv_y.close();
    vv_z.close();
    vv_dist.close();
}

TEST_CASE("Testing 3-planet Acceleration"){
    double pi=acos(-1);
    planet sun ("sun",1.,0.,0.,0.,0.,0.,0.);
    planet earth("earth",0.000003,-0.9922,-5.2975E-2,-1.444E-4,0.25679,-6.292089,-0.0005);
    planet jupiter ("jupiter",0.00095,-5.21645,-1.58016,0.12322,0.76674,-2.505676,-0.0067);
    
    solver vv("VV","sun");
    vv.add_planet(sun);
    vv.add_planet(earth);
    vv.add_planet(jupiter);

    //Test earth acceleration in all 3 directions
    double a=4*pi*pi*(1/pow(earth.distance(sun),3)+jupiter.get_mass()/pow(earth.distance(jupiter),3));
    double acc_x;
    acc_x=earth.x[0]*vv.calc_acc(earth);
    double ax=a*earth.x[0];
    REQUIRE(acc_x==Approx(ax).epsilon(0.001));
    
    double acc_y;
    acc_y=earth.x[1]*vv.calc_acc(earth);
    double ay=a*earth.x[1];
    REQUIRE(acc_y==Approx(ay).epsilon(0.001));
    
    double acc_z;
    acc_z=earth.x[2]*vv.calc_acc(earth);
    double az=a*earth.x[2];
    REQUIRE(acc_z==Approx(az).epsilon(0.001));
    
    //Test jupiter acceleration in all 3 directions
    a=4*pi*pi*(1/pow(jupiter.distance(sun),3)+earth.get_mass()/pow(earth.distance(jupiter),3));
    acc_x=jupiter.x[0]*vv.calc_acc(jupiter);
    ax=a*jupiter.x[0];
    REQUIRE(acc_x==Approx(ax).epsilon(0.001));
    
    acc_y=jupiter.x[1]*vv.calc_acc(jupiter);
    ay=a*jupiter.x[1];
    REQUIRE(acc_y==Approx(ay).epsilon(0.001));
    
    acc_z=jupiter.x[2]*vv.calc_acc(jupiter);
    az=a*jupiter.x[2];
    REQUIRE(acc_z==Approx(az).epsilon(0.001));
}

TEST_CASE("Testing 3-planet VV"){
    double h = 0.01;
    double n = 25000;
    planet sun ("sun",1.,0.,0.,0.,0.,0.,0.);
    planet earth("earth",3E-6,
            -9.921659045836120E-01,-5.297464688174806E-02,-1.444074140071765E-04,
            0.256790253188838299,-6.29207466280881475,0.000500268555676713355);
    planet jupiter("jupiter",9.5E-4,
            -5.216449805254690E+00,-1.580162951996290E+00,1.232243010564932E-01,
            0.76674377434577883,-2.5056764652495431,-0.0067372069838469784);
    solver vv("VV","sun");
    vv.add_planet(sun);
    vv.add_planet(earth);
    vv.add_planet(jupiter);
    vv.solve(n,h);
    
    //these are earth files output by code
    ifstream vv_earth_x("VV_earth_x.txt");
    ifstream vv_earth_y("VV_earth_y.txt");
    ifstream vv_earth_z("VV_earth_z.txt");
    ifstream vv_earth_dist("VV_earth_dist.txt");
    
    //test last value in each file
    double t,x,y,z,dist;
    while(vv_earth_x>>t>>x){};
    REQUIRE(x==Approx(-0.778202));
    while(vv_earth_y>>t>>y){};
    REQUIRE(y==Approx(0.616003));
    while(vv_earth_z>>t>>z){};
    REQUIRE(z==Approx(-0.000161535));
    while(vv_earth_dist>>t>>dist){};
    REQUIRE(dist==Approx(0.992501));
    
    vv_earth_x.close();
    vv_earth_y.close();
    vv_earth_z.close();
    vv_earth_dist.close();
    
    //these are jupiter files output by code
    ifstream vv_jupiter_x("VV_jupiter_x.txt");
    ifstream vv_jupiter_y("VV_jupiter_y.txt");
    ifstream vv_jupiter_z("VV_jupiter_z.txt");
    ifstream vv_jupiter_dist("VV_jupiter_dist.txt");
    
    //test last value in each file
    while(vv_jupiter_x>>t>>x){};
    REQUIRE(x==Approx(0.14351));
    while(vv_jupiter_y>>t>>y){};
    REQUIRE(y==Approx(-5.20162));
    while(vv_jupiter_z>>t>>z){};
    REQUIRE(z==Approx(0.0184012));
    while(vv_jupiter_dist>>t>>dist){};
    REQUIRE(dist==Approx(5.20364));
    
    vv_jupiter_x.close();
    vv_jupiter_y.close();
    vv_jupiter_z.close();
    vv_jupiter_dist.close();
    
}

TEST_CASE("Testing angular momentum"){
    planet sun ("sun",1.,0.,0.,0.,0.,0.,0.);
    planet earth("earth",0.000003,-0.9922,-5.2975E-2,-1.444E-4,0.25679,-6.292089,-0.0005);
    planet jupiter ("jupiter",0.00095,-5.21645,-1.58016,0.12322,0.76674,-2.505676,-0.0067);
    
    solver vv("VV","sun");
    vv.add_planet(sun);
    vv.add_planet(earth);
    vv.add_planet(jupiter);
    double L=vv.calc_ang_momentum();
    
    double px=earth.get_mass()*earth.v[0];
    double py=earth.get_mass()*earth.v[1];
    double pz=earth.get_mass()*earth.v[2];
    double Lx_e=earth.x[1]*pz-earth.x[2]*py;
    double Ly_e=-(earth.x[0]*pz-earth.x[2]*px);
    double Lz_e=earth.x[0]*py-earth.x[1]*px;
    px=jupiter.get_mass()*jupiter.v[0];
    py=jupiter.get_mass()*jupiter.v[1];
    pz=jupiter.get_mass()*jupiter.v[2];
    double Lx_j=jupiter.x[1]*pz-jupiter.x[2]*py;
    double Ly_j=-(jupiter.x[0]*pz-jupiter.x[2]*px);
    double Lz_j=jupiter.x[0]*py-jupiter.x[1]*px;
    double L_tot=sqrt((Lx_e+Lx_j)*(Lx_e+Lx_j)+(Ly_e+Ly_j)*(Ly_e+Ly_j)+(Lz_e+Lz_j)*(Lz_e+Lz_j));
    
    REQUIRE(L==Approx(L_tot));
}


