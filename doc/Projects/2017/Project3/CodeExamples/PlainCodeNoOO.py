def computerdoitforme(t, planet1, planet2, planet3, planet4, planet5, planet6, planet7, planet8, planet9):
#These are the Verlet equations for making the position value update
    def coorx(xi, h,vxi, axi):
        return xi + h*vxi+h**2*axi/2
    def coory(yi, h,vyi, ayi):
        return yi + h*vyi+h**2*ayi/2
    def coorz(zi, h,vzi, azi):
        return zi + h*vzi+h**2*azi/2
#These are the Verlet equations for updating velocity of the planet
    def velx(vxi, h, ax_i_1,ax_i):
        return vxi + (h/2)*(ax_i_1+ax_i)
    def vely(vyi, h, ay_i_1,ay_i):
        return vyi + (h/2)*(ay_i_1+ay_i)
    def velz(vzi, h, az_i_1,az_i):
        return vzi + (h/2)*(az_i_1+az_i)
#These are the acceleration equation when we include all of the planets and the separations of the planet with the 8 others
    def accx(mass1, mass2, mass3, mass4, mass5, mass6, mass7, mass8, coor1,coor2, coor3, coor4, coor5, coor6, coor7, coor8, coor9,  dist, sepdist1, sepdist2, sepdist3, sepdist4, sepdist5, sepdist6, sepdist7, sepdist8):
          return -4*math.pi**2*coor1/(dist**3)-4*math.pi**2*(coor1-coor2)*mass1/(sepdist1**3)-4*math.pi**2*(coor1-coor3)*mass2/(sepdist2**3)-4*math.pi**2*(coor1-coor4)*mass3/(sepdist3**3)-4*math.pi**2*(coor1-coor5)*mass4/(sepdist4**3)-4*math.pi**2*(coor1-coor6)*mass5/(sepdist5**3)-4*math.pi**2*(coor1-coor7)*mass6/(sepdist6**3)-4*math.pi**2*(coor1-coor8)*mass7/(sepdist7**3)-4*math.pi**2*(coor1-coor9)*mass8/(sepdist8**3)
    def accy(mass1, mass2, mass3, mass4, mass5, mass6, mass7, mass8, coor1,coor2, coor3, coor4, coor5, coor6, coor7, coor8, coor9,  dist, sepdist1, sepdist2, sepdist3, sepdist4, sepdist5, sepdist6, sepdist7, sepdist8):
          return -4*math.pi**2*coor1/(dist**3)-4*math.pi**2*(coor1-coor2)*mass1/(sepdist1**3)-4*math.pi**2*(coor1-coor3)*mass2/(sepdist2**3)-4*math.pi**2*(coor1-coor4)*mass3/(sepdist3**3)-4*math.pi**2*(coor1-coor5)*mass4/(sepdist4**3)-4*math.pi**2*(coor1-coor6)*mass5/(sepdist5**3)-4*math.pi**2*(coor1-coor7)*mass6/(sepdist6**3)-4*math.pi**2*(coor1-coor8)*mass7/(sepdist7**3)-4*math.pi**2*(coor1-coor9)*mass8/(sepdist8**3)
    def accz(mass1, mass2, mass3, mass4, mass5, mass6, mass7, mass8, coor1,coor2, coor3, coor4, coor5, coor6, coor7, coor8, coor9,  dist, sepdist1, sepdist2, sepdist3, sepdist4, sepdist5, sepdist6, sepdist7, sepdist8):
          return -4*math.pi**2*coor1/(dist**3)-4*math.pi**2*(coor1-coor2)*mass1/(sepdist1**3)-4*math.pi**2*(coor1-coor3)*mass2/(sepdist2**3)-4*math.pi**2*(coor1-coor4)*mass3/(sepdist3**3)-4*math.pi**2*(coor1-coor5)*mass4/(sepdist4**3)-4*math.pi**2*(coor1-coor6)*mass5/(sepdist5**3)-4*math.pi**2*(coor1-coor7)*mass6/(sepdist6**3)-4*math.pi**2*(coor1-coor8)*mass7/(sepdist7**3)-4*math.pi**2*(coor1-coor9)*mass8/(sepdist8**3)

#This defines the number of years, integration points, and step size
    time = t #The number of years we want to loop over
    h = 1/365 #The step size, defined as one day
    n = int(t/h) #The total numbers of iterations

    #This sets up the zero vectors for all of the quantities that will update over time
    coordinatesx1 = np.zeros(n+1)
    coordinatesy1 = np.zeros(n+1)
    coordinatesz1 = np.zeros(n+1)
    velocitiesx1 = np.zeros(n+1)
    velocitiesy1 = np.zeros(n+1)
    velocitiesz1 = np.zeros(n+1)

    coordinatesx2 = np.zeros(n+1)
    coordinatesy2 = np.zeros(n+1)
    coordinatesz2 = np.zeros(n+1)
    velocitiesx2 = np.zeros(n+1)
    velocitiesy2 = np.zeros(n+1)
    velocitiesz2 = np.zeros(n+1)

    coordinatesx3 = np.zeros(n+1)
    coordinatesy3 = np.zeros(n+1)
    coordinatesz3 = np.zeros(n+1)
    velocitiesx3 = np.zeros(n+1)
    velocitiesy3 = np.zeros(n+1)
    velocitiesz3 = np.zeros(n+1)

    coordinatesx4 = np.zeros(n+1)
    coordinatesy4 = np.zeros(n+1)
    coordinatesz4 = np.zeros(n+1)
    velocitiesx4 = np.zeros(n+1)
    velocitiesy4 = np.zeros(n+1)
    velocitiesz4 = np.zeros(n+1)

    coordinatesx5 = np.zeros(n+1)
    coordinatesy5 = np.zeros(n+1)
    coordinatesz5 = np.zeros(n+1)
    velocitiesx5 = np.zeros(n+1)
    velocitiesy5 = np.zeros(n+1)
    velocitiesz5 = np.zeros(n+1)

    coordinatesx6 = np.zeros(n+1)
    coordinatesy6 = np.zeros(n+1)
    coordinatesz6 = np.zeros(n+1)
    velocitiesx6 = np.zeros(n+1)
    velocitiesy6 = np.zeros(n+1)
    velocitiesz6 = np.zeros(n+1)

    coordinatesx7 = np.zeros(n+1)
    coordinatesy7 = np.zeros(n+1)
    coordinatesz7 = np.zeros(n+1)
    velocitiesx7 = np.zeros(n+1)
    velocitiesy7 = np.zeros(n+1)
    velocitiesz7 = np.zeros(n+1)

    coordinatesx8 = np.zeros(n+1)
    coordinatesy8 = np.zeros(n+1)
    coordinatesz8 = np.zeros(n+1)
    velocitiesx8 = np.zeros(n+1)
    velocitiesy8 = np.zeros(n+1)
    velocitiesz8 = np.zeros(n+1)

    coordinatesx9 = np.zeros(n+1)
    coordinatesy9 = np.zeros(n+1)
    coordinatesz9 = np.zeros(n+1)
    velocitiesx9 = np.zeros(n+1)
    velocitiesy9 = np.zeros(n+1)
    velocitiesz9 = np.zeros(n+1)


    #Set up initial values for the positions and the velocities
    coordinatesx1[0] = planet1.x
    coordinatesy1[0] = planet1.y
    coordinatesz1[0] = planet1.z
    velocitiesx1[0] = planet1.vx
    velocitiesy1[0] = planet1.vy
    velocitiesz1[0] = planet1.vz

    coordinatesx2[0] = planet2.x
    coordinatesy2[0] = planet2.y
    coordinatesz2[0] = planet2.z
    velocitiesx2[0] = planet2.vx
    velocitiesy2[0] = planet2.vy
    velocitiesz2[0] = planet2.vz

    coordinatesx3[0] = planet3.x
    coordinatesy3[0] = planet3.y
    coordinatesz3[0] = planet3.z
    velocitiesx3[0] = planet3.vx
    velocitiesy3[0] = planet3.vy
    velocitiesz3[0] = planet3.vz

    coordinatesx4[0] = planet4.x
    coordinatesy4[0] = planet4.y
    coordinatesz4[0] = planet4.z
    velocitiesx4[0] = planet4.vx
    velocitiesy4[0] = planet4.vy
    velocitiesz4[0] = planet4.vz

    coordinatesx5[0] = planet5.x
    coordinatesy5[0] = planet5.y
    coordinatesz5[0] = planet5.z
    velocitiesx5[0] = planet5.vx
    velocitiesy5[0] = planet5.vy
    velocitiesz5[0] = planet5.vz

    coordinatesx6[0] = planet6.x
    coordinatesy6[0] = planet6.y
    coordinatesz6[0] = planet6.z
    velocitiesx6[0] = planet6.vx
    velocitiesy6[0] = planet6.vy
    velocitiesz6[0] = planet6.vz

    coordinatesx7[0] = planet7.x
    coordinatesy7[0] = planet7.y
    coordinatesz7[0] = planet7.z
    velocitiesx7[0] = planet7.vx
    velocitiesy7[0] = planet7.vy
    velocitiesz7[0] = planet7.vz

    coordinatesx8[0] = planet8.x
    coordinatesy8[0] = planet8.y
    coordinatesz8[0] = planet8.z
    velocitiesx8[0] = planet8.vx
    velocitiesy8[0] = planet8.vy
    velocitiesz8[0] = planet8.vz

    coordinatesx9[0] = planet9.x
    coordinatesy9[0] = planet9.y
    coordinatesz9[0] = planet9.z
    velocitiesx9[0] = planet9.vx
    velocitiesy9[0] = planet9.vy
    velocitiesz9[0] = planet9.vz

#Setting up the radii for the planets defined by an earlier defined function as well as the initial values I programmed in in the beginning
    rad1 = r(planet1)
    rad2 = r(planet2)
    rad3 = r(planet3)
    rad4 = r(planet4)
    rad5 = r(planet5)
    rad6 = r(planet6)
    rad7 = r(planet7)
    rad8 = r(planet8)
    rad9 = r(planet9)

#This does the stuff. All of the stuff.
    for i in range(n):
        #Define x coordinates and velocities for all 9 planets
        x1_i = coordinatesx1[i]
        x2_i = coordinatesx2[i]
        x3_i = coordinatesx3[i]
        x4_i = coordinatesx4[i]
        x5_i = coordinatesx5[i]
        x6_i = coordinatesx6[i]
        x7_i = coordinatesx7[i]
        x8_i = coordinatesx8[i]
        x9_i = coordinatesx9[i]

        vx1_i = velocitiesx1[i]
        vx2_i = velocitiesx2[i]
        vx3_i = velocitiesx3[i]
        vx4_i = velocitiesx4[i]
        vx5_i = velocitiesx5[i]
        vx6_i = velocitiesx6[i]
        vx7_i = velocitiesx7[i]
        vx8_i = velocitiesx8[i]
        vx9_i = velocitiesx9[i]

        #Define y coordinates and velocities for all 9 planets
        y1_i = coordinatesy1[i]
        y2_i = coordinatesy2[i]
        y3_i = coordinatesy3[i]
        y4_i = coordinatesy4[i]
        y5_i = coordinatesy5[i]
        y6_i = coordinatesy6[i]
        y7_i = coordinatesy7[i]
        y8_i = coordinatesy8[i]
        y9_i = coordinatesy9[i]

        vy1_i = velocitiesy1[i]
        vy2_i = velocitiesy2[i]
        vy3_i = velocitiesy3[i]
        vy4_i = velocitiesy4[i]
        vy5_i = velocitiesy5[i]
        vy6_i = velocitiesy6[i]
        vy7_i = velocitiesy7[i]
        vy8_i = velocitiesy8[i]
        vy9_i = velocitiesy9[i]

        #Define z coordinates and velocities for all 9 planets
        z1_i = coordinatesz1[i]
        z2_i = coordinatesz2[i]
        z3_i = coordinatesz3[i]
        z4_i = coordinatesz4[i]
        z5_i = coordinatesz5[i]
        z6_i = coordinatesz6[i]
        z7_i = coordinatesz7[i]
        z8_i = coordinatesz8[i]
        z9_i = coordinatesz9[i]

        vz1_i = velocitiesz1[i]
        vz2_i = velocitiesz2[i]
        vz3_i = velocitiesz3[i]
        vz4_i = velocitiesz4[i]
        vz5_i = velocitiesz5[i]
        vz6_i = velocitiesz6[i]
        vz7_i = velocitiesz7[i]
        vz8_i = velocitiesz8[i]
        vz9_i = velocitiesz9[i]

        #Distance between the planets
        #planet1 to something else radii
        rsep12 = ((x1_i-x2_i)**2+(y1_i-y2_i)**2+(z1_i-z2_i)**2)**(1/2)
        rsep13 = ((x1_i-x3_i)**2+(y1_i-y3_i)**2+(z1_i-z3_i)**2)**(1/2)
        rsep14 = ((x1_i-x4_i)**2+(y1_i-y4_i)**2+(z1_i-z4_i)**2)**(1/2)
        rsep15 = ((x1_i-x5_i)**2+(y1_i-y5_i)**2+(z1_i-z5_i)**2)**(1/2)
        rsep16 = ((x1_i-x6_i)**2+(y1_i-y6_i)**2+(z1_i-z6_i)**2)**(1/2)
        rsep17 = ((x1_i-x7_i)**2+(y1_i-y7_i)**2+(z1_i-z7_i)**2)**(1/2)
        rsep18 = ((x1_i-x8_i)**2+(y1_i-y8_i)**2+(z1_i-z8_i)**2)**(1/2)
        rsep19 = ((x1_i-x9_i)**2+(y1_i-y9_i)**2+(z1_i-z9_i)**2)**(1/2)

        #planet2 to something else radii
        rsep21 = ((x2_i-x1_i)**2+(y2_i-y1_i)**2+(z2_i-z1_i)**2)**(1/2)
        rsep23 = ((x2_i-x3_i)**2+(y2_i-y3_i)**2+(z2_i-z3_i)**2)**(1/2)
        rsep24 = ((x2_i-x4_i)**2+(y2_i-y4_i)**2+(z2_i-z4_i)**2)**(1/2)
        rsep25 = ((x2_i-x5_i)**2+(y2_i-y5_i)**2+(z2_i-z5_i)**2)**(1/2)
        rsep26 = ((x2_i-x6_i)**2+(y2_i-y6_i)**2+(z2_i-z6_i)**2)**(1/2)
        rsep27 = ((x2_i-x7_i)**2+(y2_i-y7_i)**2+(z2_i-z7_i)**2)**(1/2)
        rsep28 = ((x2_i-x8_i)**2+(y2_i-y8_i)**2+(z2_i-z8_i)**2)**(1/2)
        rsep29 = ((x2_i-x9_i)**2+(y2_i-y9_i)**2+(z2_i-z9_i)**2)**(1/2)

        #planet3 to something else radii
        rsep31 = ((x3_i-x1_i)**2+(y3_i-y1_i)**2+(z3_i-z1_i)**2)**(1/2)
        rsep32 = ((x3_i-x2_i)**2+(y3_i-y2_i)**2+(z3_i-z2_i)**2)**(1/2)
        rsep34 = ((x3_i-x4_i)**2+(y3_i-y4_i)**2+(z3_i-z4_i)**2)**(1/2)
        rsep35 = ((x3_i-x5_i)**2+(y3_i-y5_i)**2+(z3_i-z5_i)**2)**(1/2)
        rsep36 = ((x3_i-x6_i)**2+(y3_i-y6_i)**2+(z3_i-z6_i)**2)**(1/2)
        rsep37 = ((x3_i-x7_i)**2+(y3_i-y7_i)**2+(z3_i-z7_i)**2)**(1/2)
        rsep38 = ((x3_i-x8_i)**2+(y3_i-y8_i)**2+(z3_i-z8_i)**2)**(1/2)
        rsep39 = ((x3_i-x9_i)**2+(y3_i-y9_i)**2+(z3_i-z9_i)**2)**(1/2)

        #planet4 to something else radii
        rsep41 = ((x4_i-x1_i)**2+(y4_i-y1_i)**2+(z4_i-z1_i)**2)**(1/2)
        rsep42 = ((x4_i-x2_i)**2+(y4_i-y2_i)**2+(z4_i-z2_i)**2)**(1/2)
        rsep43 = ((x4_i-x3_i)**2+(y4_i-y3_i)**2+(z4_i-z3_i)**2)**(1/2)
        rsep45 = ((x4_i-x5_i)**2+(y4_i-y5_i)**2+(z4_i-z5_i)**2)**(1/2)
        rsep46 = ((x4_i-x6_i)**2+(y4_i-y6_i)**2+(z4_i-z6_i)**2)**(1/2)
        rsep47 = ((x4_i-x7_i)**2+(y4_i-y7_i)**2+(z4_i-z7_i)**2)**(1/2)
        rsep48 = ((x4_i-x8_i)**2+(y4_i-y8_i)**2+(z4_i-z8_i)**2)**(1/2)
        rsep49 = ((x4_i-x9_i)**2+(y4_i-y9_i)**2+(z4_i-z9_i)**2)**(1/2)

        #planet5 to something else radii
        rsep51 = ((x5_i-x1_i)**2+(y5_i-y1_i)**2+(z5_i-z1_i)**2)**(1/2)
        rsep52 = ((x5_i-x2_i)**2+(y5_i-y2_i)**2+(z5_i-z2_i)**2)**(1/2)
        rsep53 = ((x5_i-x3_i)**2+(y5_i-y3_i)**2+(z5_i-z3_i)**2)**(1/2)
        rsep54 = ((x5_i-x4_i)**2+(y5_i-y4_i)**2+(z5_i-z4_i)**2)**(1/2)
        rsep56 = ((x5_i-x6_i)**2+(y5_i-y6_i)**2+(z5_i-z6_i)**2)**(1/2)
        rsep57 = ((x5_i-x7_i)**2+(y5_i-y7_i)**2+(z5_i-z7_i)**2)**(1/2)
        rsep58 = ((x5_i-x8_i)**2+(y5_i-y8_i)**2+(z5_i-z8_i)**2)**(1/2)
        rsep59 = ((x5_i-x9_i)**2+(y5_i-y9_i)**2+(z5_i-z9_i)**2)**(1/2)

        #planet6 to something else radii
        rsep61 = ((x6_i-x1_i)**2+(y6_i-y1_i)**2+(z6_i-z1_i)**2)**(1/2)
        rsep62 = ((x6_i-x2_i)**2+(y6_i-y2_i)**2+(z6_i-z2_i)**2)**(1/2)
        rsep63 = ((x6_i-x3_i)**2+(y6_i-y3_i)**2+(z6_i-z3_i)**2)**(1/2)
        rsep64 = ((x6_i-x4_i)**2+(y6_i-y4_i)**2+(z6_i-z4_i)**2)**(1/2)
        rsep65 = ((x6_i-x5_i)**2+(y6_i-y5_i)**2+(z6_i-z5_i)**2)**(1/2)
        rsep67 = ((x6_i-x7_i)**2+(y6_i-y7_i)**2+(z6_i-z7_i)**2)**(1/2)
        rsep68 = ((x6_i-x8_i)**2+(y6_i-y8_i)**2+(z6_i-z8_i)**2)**(1/2)
        rsep69 = ((x6_i-x9_i)**2+(y6_i-y9_i)**2+(z6_i-z9_i)**2)**(1/2)

        #planet7 to something else radii
        rsep71 = ((x7_i-x1_i)**2+(y7_i-y1_i)**2+(z7_i-z1_i)**2)**(1/2)
        rsep72 = ((x7_i-x2_i)**2+(y7_i-y2_i)**2+(z7_i-z2_i)**2)**(1/2)
        rsep73 = ((x7_i-x3_i)**2+(y7_i-y3_i)**2+(z7_i-z3_i)**2)**(1/2)
        rsep74 = ((x7_i-x4_i)**2+(y7_i-y4_i)**2+(z7_i-z4_i)**2)**(1/2)
        rsep75 = ((x7_i-x5_i)**2+(y7_i-y5_i)**2+(z7_i-z5_i)**2)**(1/2)
        rsep76 = ((x7_i-x6_i)**2+(y7_i-y6_i)**2+(z7_i-z6_i)**2)**(1/2)
        rsep78 = ((x7_i-x8_i)**2+(y7_i-y8_i)**2+(z7_i-z8_i)**2)**(1/2)
        rsep79 = ((x7_i-x9_i)**2+(y7_i-y9_i)**2+(z7_i-z9_i)**2)**(1/2)

        #planet8 to something else radii
        rsep81 = ((x8_i-x1_i)**2+(y8_i-y1_i)**2+(z8_i-z1_i)**2)**(1/2)
        rsep82 = ((x8_i-x2_i)**2+(y8_i-y2_i)**2+(z8_i-z2_i)**2)**(1/2)
        rsep83 = ((x8_i-x3_i)**2+(y8_i-y3_i)**2+(z8_i-z3_i)**2)**(1/2)
        rsep84 = ((x8_i-x4_i)**2+(y8_i-y4_i)**2+(z8_i-z4_i)**2)**(1/2)
        rsep85 = ((x8_i-x5_i)**2+(y8_i-y5_i)**2+(z8_i-z5_i)**2)**(1/2)
        rsep86 = ((x8_i-x6_i)**2+(y8_i-y6_i)**2+(z8_i-z6_i)**2)**(1/2)
        rsep87 = ((x8_i-x7_i)**2+(y8_i-y7_i)**2+(z8_i-z7_i)**2)**(1/2)
        rsep89 = ((x8_i-x9_i)**2+(y8_i-y9_i)**2+(z8_i-z9_i)**2)**(1/2)

        #planet9 to something else radii
        rsep91 = ((x9_i-x1_i)**2+(y9_i-y1_i)**2+(z9_i-z1_i)**2)**(1/2)
        rsep92 = ((x9_i-x2_i)**2+(y9_i-y2_i)**2+(z9_i-z2_i)**2)**(1/2)
        rsep93 = ((x9_i-x3_i)**2+(y9_i-y3_i)**2+(z9_i-z3_i)**2)**(1/2)
        rsep94 = ((x9_i-x4_i)**2+(y9_i-y4_i)**2+(z9_i-z4_i)**2)**(1/2)
        rsep95 = ((x9_i-x5_i)**2+(y9_i-y5_i)**2+(z9_i-z5_i)**2)**(1/2)
        rsep96 = ((x9_i-x6_i)**2+(y9_i-y6_i)**2+(z9_i-z6_i)**2)**(1/2)
        rsep97 = ((x9_i-x7_i)**2+(y9_i-y7_i)**2+(z9_i-z7_i)**2)**(1/2)
        rsep98 = ((x9_i-x8_i)**2+(y9_i-y8_i)**2+(z9_i-z8_i)**2)**(1/2)

        #x accelerations including interactions from the other planets
        ax1_i = accx(planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, x1_i, x2_i, x3_i, x4_i, x5_i, x6_i, x7_i, x8_i, x9_i, rad1, rsep12, rsep13, rsep14, rsep15, rsep16, rsep17, rsep18, rsep19)
        ax2_i = accx(planet1.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, x2_i, x1_i, x3_i, x4_i, x5_i, x6_i, x7_i, x8_i, x9_i, rad2, rsep21, rsep23, rsep24, rsep25, rsep26, rsep27, rsep28, rsep29)
        ax3_i = accx(planet1.mass, planet2.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, x3_i, x1_i, x2_i, x4_i, x5_i, x6_i, x7_i, x8_i, x9_i, rad3, rsep31, rsep32, rsep34, rsep35, rsep36, rsep37, rsep38, rsep39)
        ax4_i = accx(planet1.mass, planet2.mass, planet3.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, x4_i, x1_i, x2_i, x3_i, x5_i, x6_i, x7_i, x8_i, x9_i, rad4, rsep41, rsep42, rsep43, rsep45, rsep46, rsep47, rsep48, rsep49)
        ax5_i = accx(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, x5_i, x1_i, x2_i, x3_i, x4_i, x6_i, x7_i, x8_i, x9_i, rad5, rsep51, rsep52, rsep53, rsep54, rsep56, rsep57, rsep58, rsep59)
        ax6_i = accx(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet7.mass, planet8.mass, planet9.mass, x6_i, x1_i, x2_i, x3_i, x4_i, x5_i, x7_i, x8_i, x9_i, rad6, rsep61, rsep62, rsep63, rsep64, rsep65, rsep67, rsep68, rsep69)
        ax7_i = accx(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet8.mass, planet9.mass, x7_i, x1_i, x2_i, x3_i, x4_i, x5_i, x6_i, x8_i, x9_i, rad7, rsep71, rsep72, rsep73, rsep74, rsep75, rsep76, rsep78, rsep79)
        ax8_i = accx(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet9.mass, x8_i, x1_i, x2_i, x3_i, x4_i, x5_i, x6_i, x7_i, x9_i, rad8, rsep81, rsep82, rsep83, rsep84, rsep85, rsep86, rsep87, rsep89)
        ax9_i = accx(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, x9_i, x1_i, x2_i, x3_i, x4_i, x5_i, x6_i, x7_i, x8_i, rad9, rsep91, rsep92, rsep93, rsep94, rsep95, rsep96, rsep97, rsep98)

        #y accelerations including interactions from the other planets
        ay1_i = accy(planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, y1_i, y2_i, y3_i, y4_i, y5_i, y6_i, y7_i, y8_i, y9_i, rad1, rsep12, rsep13, rsep14, rsep15, rsep16, rsep17, rsep18, rsep19)
        ay2_i = accy(planet1.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, y2_i, y1_i, y3_i, y4_i, y5_i, y6_i, y7_i, y8_i, y9_i, rad2, rsep21, rsep23, rsep24, rsep25, rsep26, rsep27, rsep28, rsep29)
        ay3_i = accy(planet1.mass, planet2.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, y3_i, y1_i, y2_i, y4_i, y5_i, y6_i, y7_i, y8_i, y9_i, rad3, rsep31, rsep32, rsep34, rsep35, rsep36, rsep37, rsep38, rsep39)
        ay4_i = accy(planet1.mass, planet2.mass, planet3.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, y4_i, y1_i, y2_i, y3_i, y5_i, y6_i, y7_i, y8_i, y9_i, rad4, rsep41, rsep42, rsep43, rsep45, rsep46, rsep47, rsep48, rsep49)
        ay5_i = accy(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, y5_i, y1_i, y2_i, y3_i, y4_i, y6_i, y7_i, y8_i, y9_i, rad5, rsep51, rsep52, rsep53, rsep54, rsep56, rsep57, rsep58, rsep59)
        ay6_i = accy(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet7.mass, planet8.mass, planet9.mass, y6_i, y1_i, y2_i, y3_i, y4_i, y5_i, y7_i, y8_i, y9_i, rad6, rsep61, rsep62, rsep63, rsep64, rsep65, rsep67, rsep68, rsep69)
        ay7_i = accy(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet8.mass, planet9.mass, y7_i, y1_i, y2_i, y3_i, y4_i, y5_i, y6_i, y8_i, y9_i, rad7, rsep71, rsep72, rsep73, rsep74, rsep75, rsep76, rsep78, rsep79)
        ay8_i = accy(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet9.mass, y8_i, y1_i, y2_i, y3_i, y4_i, y5_i, y6_i, y7_i, y9_i, rad8, rsep81, rsep82, rsep83, rsep84, rsep85, rsep86, rsep87, rsep89)
        ay9_i = accy(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, y9_i, y1_i, y2_i, y3_i, y4_i, y5_i, y6_i, y7_i, y8_i, rad9, rsep91, rsep92, rsep93, rsep94, rsep95, rsep96, rsep97, rsep98)

        #z accelerations including interactions from the other planets
        az1_i = accz(planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, z1_i, z2_i, z3_i, z4_i, z5_i, z6_i, z7_i, z8_i, z9_i, rad1, rsep12, rsep13, rsep14, rsep15, rsep16, rsep17, rsep18, rsep19)
        az2_i = accz(planet1.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, z2_i, z1_i, z3_i, z4_i, z5_i, z6_i, z7_i, z8_i, z9_i, rad2, rsep21, rsep23, rsep24, rsep25, rsep26, rsep27, rsep28, rsep29)
        az3_i = accz(planet1.mass, planet2.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, z3_i, z1_i, z2_i, z4_i, z5_i, z6_i, z7_i, z8_i, z9_i, rad3, rsep31, rsep32, rsep34, rsep35, rsep36, rsep37, rsep38, rsep39)
        az4_i = accz(planet1.mass, planet2.mass, planet3.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, z4_i, z1_i, z2_i, z3_i, z5_i, z6_i, z7_i, z8_i, z9_i, rad4, rsep41, rsep42, rsep43, rsep45, rsep46, rsep47, rsep48, rsep49)
        az5_i = accz(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, z5_i, z1_i, z2_i, z3_i, z4_i, z6_i, z7_i, z8_i, z9_i, rad5, rsep51, rsep52, rsep53, rsep54, rsep56, rsep57, rsep58, rsep59)
        az6_i = accz(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet7.mass, planet8.mass, planet9.mass, z6_i, z1_i, z2_i, z3_i, z4_i, z5_i, z7_i, z8_i, z9_i, rad6, rsep61, rsep62, rsep63, rsep64, rsep65, rsep67, rsep68, rsep69)
        az7_i = accz(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet8.mass, planet9.mass, z7_i, z1_i, z2_i, z3_i, z4_i, z5_i, z6_i, z8_i, z9_i, rad7, rsep71, rsep72, rsep73, rsep74, rsep75, rsep76, rsep78, rsep79)
        az8_i = accz(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet9.mass, z8_i, z1_i, z2_i, z3_i, z4_i, z5_i, z6_i, z7_i, z9_i, rad8, rsep81, rsep82, rsep83, rsep84, rsep85, rsep86, rsep87, rsep89)
        az9_i = accz(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, z9_i, z1_i, z2_i, z3_i, z4_i, z5_i, z6_i, z7_i, z8_i, rad9, rsep91, rsep92, rsep93, rsep94, rsep95, rsep96, rsep97, rsep98)

        #next x coordinate in the series
        x1_i_1 = coorx(x1_i, h,vx1_i,ax1_i)
        x2_i_1 = coorx(x2_i, h,vx2_i,ax2_i)
        x3_i_1 = coorx(x3_i, h,vx3_i,ax3_i)
        x4_i_1 = coorx(x4_i, h,vx4_i,ax4_i)
        x5_i_1 = coorx(x5_i, h,vx5_i,ax5_i)
        x6_i_1 = coorx(x6_i, h,vx6_i,ax6_i)
        x7_i_1 = coorx(x7_i, h,vx7_i,ax7_i)
        x8_i_1 = coorx(x8_i, h,vx8_i,ax8_i)
        x9_i_1 = coorx(x9_i, h,vx9_i,ax9_i)
        #next y coordinate in the series
        y1_i_1 = coory(y1_i, h,vy1_i,ay1_i)
        y2_i_1 = coory(y2_i, h,vy2_i,ay2_i)
        y3_i_1 = coory(y3_i, h,vy3_i,ay3_i)
        y4_i_1 = coory(y4_i, h,vy4_i,ay4_i)
        y5_i_1 = coory(y5_i, h,vy5_i,ay5_i)
        y6_i_1 = coory(y6_i, h,vy6_i,ay6_i)
        y7_i_1 = coory(y7_i, h,vy7_i,ay7_i)
        y8_i_1 = coory(y8_i, h,vy8_i,ay8_i)
        y9_i_1 = coory(y9_i, h,vy9_i,ay9_i)
        #next z coordinate in the series
        z1_i_1 = coorz(z1_i, h,vz1_i,az1_i)
        z2_i_1 = coorz(z2_i, h,vz2_i,az2_i)
        z3_i_1 = coorz(z3_i, h,vz3_i,az3_i)
        z4_i_1 = coorz(z4_i, h,vz4_i,az4_i)
        z5_i_1 = coorz(z5_i, h,vz5_i,az5_i)
        z6_i_1 = coorz(z6_i, h,vz6_i,az6_i)
        z7_i_1 = coorz(z7_i, h,vz7_i,az7_i)
        z8_i_1 = coorz(z8_i, h,vz8_i,az8_i)
        z9_i_1 = coorz(z9_i, h,vz9_i,az9_i)

        #updates the lists that I made previously
        coordinatesx1[i+1] = x1_i_1
        coordinatesx2[i+1] = x2_i_1
        coordinatesx3[i+1] = x3_i_1
        coordinatesx4[i+1] = x4_i_1
        coordinatesx5[i+1] = x5_i_1
        coordinatesx6[i+1] = x6_i_1
        coordinatesx7[i+1] = x7_i_1
        coordinatesx8[i+1] = x8_i_1
        coordinatesx9[i+1] = x9_i_1

        coordinatesy1[i+1] = y1_i_1
        coordinatesy2[i+1] = y2_i_1
        coordinatesy3[i+1] = y3_i_1
        coordinatesy4[i+1] = y4_i_1
        coordinatesy5[i+1] = y5_i_1
        coordinatesy6[i+1] = y6_i_1
        coordinatesy7[i+1] = y7_i_1
        coordinatesy8[i+1] = y8_i_1
        coordinatesy9[i+1] = y9_i_1

        coordinatesz1[i+1] = z1_i_1
        coordinatesz2[i+1] = z2_i_1
        coordinatesz3[i+1] = z3_i_1
        coordinatesz4[i+1] = z4_i_1
        coordinatesz5[i+1] = z5_i_1
        coordinatesz6[i+1] = z6_i_1
        coordinatesz7[i+1] = z7_i_1
        coordinatesz8[i+1] = z8_i_1
        coordinatesz9[i+1] = z9_i_1

        #accelaration updates for the next time step
        ax1_i_1 = accx(planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, x1_i_1, x2_i_1, x3_i_1, x4_i_1, x5_i_1, x6_i_1, x7_i_1, x8_i_1, x9_i_1, rad1, rsep12, rsep13, rsep14, rsep15, rsep16, rsep17, rsep18, rsep19)
        ax2_i_1 = accx(planet1.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, x2_i_1, x1_i_1, x3_i_1, x4_i_1, x5_i_1, x6_i_1, x7_i_1, x8_i_1, x9_i_1, rad2, rsep12, rsep23, rsep24, rsep25, rsep26, rsep27, rsep28, rsep29)
        ax3_i_1 = accx(planet1.mass, planet2.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, x3_i_1, x1_i_1, x2_i_1, x4_i_1, x5_i_1, x6_i_1, x7_i_1, x8_i_1, x9_i_1, rad3, rsep13, rsep23, rsep34, rsep35, rsep36, rsep37, rsep38, rsep39)
        ax4_i_1 = accx(planet1.mass, planet2.mass, planet3.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, x4_i_1, x1_i_1, x2_i_1, x3_i_1, x5_i_1, x6_i_1, x7_i_1, x8_i_1, x9_i_1, rad4, rsep14, rsep24, rsep34, rsep45, rsep46, rsep47, rsep48, rsep49)
        ax5_i_1 = accx(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, x5_i_1, x1_i_1, x2_i_1, x3_i_1, x4_i_1, x6_i_1, x7_i_1, x8_i_1, x9_i_1, rad5, rsep15, rsep25, rsep35, rsep45, rsep56, rsep57, rsep58, rsep59)
        ax6_i_1 = accx(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet7.mass, planet8.mass, planet9.mass, x6_i_1, x1_i_1, x2_i_1, x3_i_1, x4_i_1, x5_i_1, x7_i_1, x8_i_1, x9_i_1, rad6, rsep16, rsep26, rsep36, rsep46, rsep56, rsep67, rsep68, rsep69)
        ax7_i_1 = accx(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet8.mass, planet9.mass, x7_i_1, x1_i_1, x2_i_1, x3_i_1, x4_i_1, x5_i_1, x6_i_1, x8_i_1, x9_i_1, rad7, rsep17, rsep27, rsep37, rsep47, rsep57, rsep67, rsep78, rsep79)
        ax8_i_1 = accx(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet9.mass, x8_i_1, x1_i_1, x2_i_1, x3_i_1, x4_i_1, x5_i_1, x6_i_1, x7_i_1, x9_i_1, rad8, rsep18, rsep28, rsep38, rsep48, rsep58, rsep68, rsep78, rsep89)
        ax9_i_1 = accx(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, x9_i_1, x1_i_1, x2_i_1, x3_i_1, x4_i_1, x5_i_1, x6_i_1, x7_i_1, x8_i_1, rad9, rsep19, rsep29, rsep39, rsep49, rsep59, rsep69, rsep79, rsep89)

        ay1_i_1 = accy(planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, y1_i_1, y2_i_1, y3_i_1, y4_i_1, y5_i_1, y6_i_1, y7_i_1, y8_i_1, y9_i_1, rad1, rsep12, rsep13, rsep14, rsep15, rsep16, rsep17, rsep18, rsep19)
        ay2_i_1 = accy(planet1.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, y2_i_1, y1_i_1, y3_i_1, y4_i_1, y5_i_1, y6_i_1, y7_i_1, y8_i_1, y9_i_1, rad2, rsep12, rsep23, rsep24, rsep25, rsep26, rsep27, rsep28, rsep29)
        ay3_i_1 = accy(planet1.mass, planet2.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, y3_i_1, y1_i_1, y2_i_1, y4_i_1, y5_i_1, y6_i_1, y7_i_1, y8_i_1, y9_i_1, rad3, rsep13, rsep23, rsep34, rsep35, rsep36, rsep37, rsep38, rsep39)
        ay4_i_1 = accy(planet1.mass, planet2.mass, planet3.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, y4_i_1, y1_i_1, y2_i_1, y3_i_1, y5_i_1, y6_i_1, y7_i_1, y8_i_1, y9_i_1, rad4, rsep14, rsep24, rsep34, rsep45, rsep46, rsep47, rsep48, rsep49)
        ay5_i_1 = accy(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, y5_i_1, y1_i_1, y2_i_1, y3_i_1, y4_i_1, y6_i_1, y7_i_1, y8_i_1, y9_i_1, rad5, rsep15, rsep25, rsep35, rsep45, rsep56, rsep57, rsep58, rsep59)
        ay6_i_1 = accy(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet7.mass, planet8.mass, planet9.mass, y6_i_1, y1_i_1, y2_i_1, y3_i_1, y4_i_1, y5_i_1, y7_i_1, y8_i_1, y9_i_1, rad6, rsep16, rsep26, rsep36, rsep46, rsep56, rsep67, rsep68, rsep69)
        ay7_i_1 = accy(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet8.mass, planet9.mass, y7_i_1, y1_i_1, y2_i_1, y3_i_1, y4_i_1, y5_i_1, y6_i_1, y8_i_1, y9_i_1, rad7, rsep17, rsep27, rsep37, rsep47, rsep57, rsep67, rsep78, rsep79)
        ay8_i_1 = accy(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet9.mass, y8_i_1, y1_i_1, y2_i_1, y3_i_1, y4_i_1, y5_i_1, y6_i_1, y7_i_1, y9_i_1, rad8, rsep18, rsep28, rsep38, rsep48, rsep58, rsep68, rsep78, rsep89)
        ay9_i_1 = accy(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, y9_i_1, y1_i_1, y2_i_1, y3_i_1, y4_i_1, y5_i_1, y6_i_1, y7_i_1, y8_i_1, rad9, rsep19, rsep29, rsep39, rsep49, rsep59, rsep69, rsep79, rsep89)

        az1_i_1 = accz(planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, z1_i_1, z2_i_1, z3_i_1, z4_i_1, z5_i_1, z6_i_1, z7_i_1, z8_i_1, z9_i_1, rad1, rsep12, rsep13, rsep14, rsep15, rsep16, rsep17, rsep18, rsep19)
        az2_i_1 = accz(planet1.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, z2_i_1, z1_i_1, z3_i_1, z4_i_1, z5_i_1, z6_i_1, z7_i_1, z8_i_1, z9_i_1, rad2, rsep12, rsep23, rsep24, rsep25, rsep26, rsep27, rsep28, rsep29)
        az3_i_1 = accz(planet1.mass, planet2.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, z3_i_1, z1_i_1, z2_i_1, z4_i_1, z5_i_1, z6_i_1, z7_i_1, z8_i_1, z9_i_1, rad3, rsep13, rsep23, rsep34, rsep35, rsep36, rsep37, rsep38, rsep39)
        az4_i_1 = accz(planet1.mass, planet2.mass, planet3.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, z4_i_1, z1_i_1, z2_i_1, z3_i_1, z5_i_1, z6_i_1, z7_i_1, z8_i_1, z9_i_1, rad4, rsep14, rsep24, rsep34, rsep45, rsep46, rsep47, rsep48, rsep49)
        az5_i_1 = accz(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet6.mass, planet7.mass, planet8.mass, planet9.mass, z5_i_1, z1_i_1, z2_i_1, z3_i_1, z4_i_1, z6_i_1, z7_i_1, z8_i_1, z9_i_1, rad5, rsep15, rsep25, rsep35, rsep45, rsep56, rsep57, rsep58, rsep59)
        az6_i_1 = accz(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet7.mass, planet8.mass, planet9.mass, z6_i_1, z1_i_1, z2_i_1, z3_i_1, z4_i_1, z5_i_1, z7_i_1, z8_i_1, z9_i_1, rad6, rsep16, rsep26, rsep36, rsep46, rsep56, rsep67, rsep68, rsep69)
        az7_i_1 = accz(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet8.mass, planet9.mass, z7_i_1, z1_i_1, z2_i_1, z3_i_1, z4_i_1, z5_i_1, z6_i_1, z8_i_1, z9_i_1, rad7, rsep17, rsep27, rsep37, rsep47, rsep57, rsep67, rsep78, rsep79)
        az8_i_1 = accz(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet9.mass, z8_i_1, z1_i_1, z2_i_1, z3_i_1, z4_i_1, z5_i_1, z6_i_1, z7_i_1, z9_i_1, rad8, rsep18, rsep28, rsep38, rsep48, rsep58, rsep68, rsep78, rsep89)
        az9_i_1 = accz(planet1.mass, planet2.mass, planet3.mass, planet4.mass, planet5.mass, planet6.mass, planet7.mass, planet8.mass, z9_i_1, z1_i_1, z2_i_1, z3_i_1, z4_i_1, z5_i_1, z6_i_1, z7_i_1, z8_i_1, rad9, rsep19, rsep29, rsep39, rsep49, rsep59, rsep69, rsep79, rsep89)

        #update velocities for the next time step
        vx1_i_1 = velx(vx1_i,h,ax1_i_1,ax1_i)
        vx2_i_1 = velx(vx2_i,h,ax2_i_1,ax2_i)
        vx3_i_1 = velx(vx3_i,h,ax3_i_1,ax3_i)
        vx4_i_1 = velx(vx4_i,h,ax4_i_1,ax4_i)
        vx5_i_1 = velx(vx5_i,h,ax5_i_1,ax5_i)
        vx6_i_1 = velx(vx6_i,h,ax6_i_1,ax6_i)
        vx7_i_1 = velx(vx7_i,h,ax7_i_1,ax7_i)
        vx8_i_1 = velx(vx8_i,h,ax8_i_1,ax8_i)
        vx9_i_1 = velx(vx9_i,h,ax9_i_1,ax9_i)

        vy1_i_1 = vely(vy1_i,h,ay1_i_1,ay1_i)
        vy2_i_1 = vely(vy2_i,h,ay2_i_1,ay2_i)
        vy3_i_1 = vely(vy3_i,h,ay3_i_1,ay3_i)
        vy4_i_1 = vely(vy4_i,h,ay4_i_1,ay4_i)
        vy5_i_1 = vely(vy5_i,h,ay5_i_1,ay5_i)
        vy6_i_1 = vely(vy6_i,h,ay6_i_1,ay6_i)
        vy7_i_1 = vely(vy7_i,h,ay7_i_1,ay7_i)
        vy8_i_1 = vely(vy8_i,h,ay8_i_1,ay8_i)
        vy9_i_1 = vely(vy9_i,h,ay9_i_1,ay9_i)

        vz1_i_1 = velz(vz1_i,h,az1_i_1,az1_i)
        vz2_i_1 = velz(vz2_i,h,az2_i_1,az2_i)
        vz3_i_1 = velz(vz3_i,h,az3_i_1,az3_i)
        vz4_i_1 = velz(vz4_i,h,az4_i_1,az4_i)
        vz5_i_1 = velz(vz5_i,h,az5_i_1,az5_i)
        vz6_i_1 = velz(vz6_i,h,az6_i_1,az6_i)
        vz7_i_1 = velz(vz7_i,h,az7_i_1,az7_i)
        vz8_i_1 = velz(vz8_i,h,az8_i_1,az8_i)
        vz9_i_1 = velz(vz9_i,h,az9_i_1,az9_i)

        #updating the lists made earlier
        velocitiesx1[i+1] = vx1_i_1
        velocitiesx2[i+1] = vx2_i_1
        velocitiesx3[i+1] = vx3_i_1
        velocitiesx4[i+1] = vx4_i_1
        velocitiesx5[i+1] = vx5_i_1
        velocitiesx6[i+1] = vx6_i_1
        velocitiesx7[i+1] = vx7_i_1
        velocitiesx8[i+1] = vx8_i_1
        velocitiesx9[i+1] = vx9_i_1

        velocitiesy1[i+1] = vy1_i_1
        velocitiesy2[i+1] = vy2_i_1
        velocitiesy3[i+1] = vy3_i_1
        velocitiesy4[i+1] = vy4_i_1
        velocitiesy5[i+1] = vy5_i_1
        velocitiesy6[i+1] = vy6_i_1
        velocitiesy7[i+1] = vy7_i_1
        velocitiesy8[i+1] = vy8_i_1
        velocitiesy9[i+1] = vy9_i_1

        velocitiesz1[i+1] = vz1_i_1
        velocitiesz2[i+1] = vz2_i_1
        velocitiesz3[i+1] = vz3_i_1
        velocitiesz4[i+1] = vz4_i_1
        velocitiesz5[i+1] = vz5_i_1
        velocitiesz6[i+1] = vz6_i_1
        velocitiesz7[i+1] = vz7_i_1
        velocitiesz8[i+1] = vz8_i_1
        velocitiesz9[i+1] = vz9_i_1

    return (coordinatesx1, coordinatesx2, coordinatesx3, coordinatesx4, coordinatesx5, coordinatesx6, coordinatesx7, coordinatesx8, coordinatesx9), (coordinatesy1, coordinatesy2, coordinatesy3, coordinatesy4, coordinatesy5, coordinatesy6, coordinatesy7, coordinatesy8, coordinatesy9), (coordinatesz1, coordinatesz2, coordinatesz3, coordinatesz4, coordinatesz5, coordinatesz6, coordinatesz7, coordinatesz8, coordinatesz9)
