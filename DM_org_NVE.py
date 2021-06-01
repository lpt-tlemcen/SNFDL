import math

N = 500                   # Nombre de particules
Ro = 0.85                 # Densité réduite
Tr = 2                    # Température réduite
ncfc=round((N/4.)**(1/3))   # Nombre de petites boîtes CFC
L = (N/Ro)**(1/3)         # Longueur réduite
Ls2 = L/2.                # Moitié de la longueur de la boîtes256256

Niter = 10               # Nombre d'itération

# Paramètre du Potentiel Lennard Jones Cubic Spline LJCS , Réf. :
# Bjørn Hafskjold, Karl Patrick Travis, Amanda Bailey Hass, Morten Hammer, Ailo Aasen & Øivind Wilhelmsen (2019)
# Thermodynamic properties of the 3D Lennard-Jones/spline model, Molecular Physics, DOI: 10.1080/00268976.2019.1664780

Aljs = -24192/3211       # Paramètre du Potentiel Lennard Jones Cubic Spline LJCS
Bljs = -387072/61009     # Paramètre du Potentiel Lennard Jones Cubic Spline LJCS
rsljs = (26/7)**(1/6)    # Paramètre du Potentiel Lennard Jones Cubic Spline LJCS
rcljs = 67/48*rsljs      # Paramètre du Potentiel Lennard Jones Cubic Spline LJCS

dt = 0.001              # Delta(t) Pas de temps en unité réduite
dts2 = dt/2             # Delta(t)/2 Moitié du pas de temps en unité réduite

print("N=     ", N)
print("Ro=    ", Ro)
print("Tr=    ", Tr)
print("ncfc = ", ncfc)
print("L = "   , L)

print("Aljs=     ", Aljs)
print("Aljs=    ", Bljs)
print("rsljs=    ", rsljs)
print("rcljs=    ", rcljs)

#pour la premier boite (initialisation)

# itinialisation des tableaux des positions 

x=[]
y=[]
z=[]

#les positions initiales :

# on va déplacer notre petite boite CFC avec ses 04 atomes dans les 03 axes X,Y,Z

for i in range (ncfc): #déplacement de la petite boite CFC selon X
    for j in range(ncfc): #déplacement de la petite boite CFC selon Y
        for k in range(ncfc): #déplacement de la petite boite CFC selon Z
            x.append(0+i)  #x atome 1 de la petite boite CFC
            y.append(0+j)  #y atome 1 de la petite boite CFC
            z.append(0+k)  #z atome 1 de la petite boite CFC
            x.append(1/2+i) #x atome 2 de la petite boite CFC
            y.append(1/2+j) #y atome 2 de la petite boite CFC
            z.append(0+k)   #z atome 2 de la petite boite CFC
            x.append(1/2+i) #x atome 3 de la petite boite CFC
            y.append(0+j)   #y atome 3 de la petite boite CFC
            z.append(1/2+k) #z atome 3 de la petite boite CFC
            x.append(0+i)   #x atome 4 de la petite boite CFC
            y.append(1/2+j) #y atome 4 de la petite boite CFC
            z.append(1/2+k) #z atome 4 de la petite boite CFC

# Un rescaling spatial (*L/ncfc) pour s'assurer que ce soient les coordonnées "réelles" des particules à l'intérieur de la boîte cubique de longeur L

for i in range (N) :
    x[i] = x[i]*L/ncfc
    y[i] = y[i]*L/ncfc
    z[i] = z[i]*L/ncfc

#les vitesses initiales :

# appeler le générateur de nombres pseudo-aléatoires
from random import seed
from random import random

# initier le  générateur de nombre pseudo-aléatoire
seed(1)

# itinialisation des tableaux des positions 

vx=[]
vy=[]
vz=[]

# on va générer les vitesses de nos particules extraite d'une distribution de Maxwell-Boltzmann à T = Tr pour les 03 composantes vx, vy & vz

for i in range (N):
    Ax = random() #nombre pseudo aléatoire
    vxi = math.sqrt(-2 * Tr * math.log(Ax)) #composante x de la vitesse extraite de la distribution de Maxwell-Boltzmann
    vx.append(vxi)
    Ay=random() #nombre pseudo aléatoire
    vyi=math.sqrt(-2*Tr*math.log(Ay)) #composante y de la vitesse extraite de la distribution de Maxwell-Boltzmann
    vy.append(vyi)
    Az=random() #nombre pseudo aléatoire
    vzi=math.sqrt(-2*Tr*math.log(Az)) #composante z de la vitesse extraite de la distribution de Maxwell-Boltzmann
    vz.append(vzi)

# On va s'assurer que la quantité totale de mouvement soit nulle 

vcmx = sum(vx)/N
vcmy = sum(vy)/N
vcmz = sum(vz)/N

S = 0

for i in range (N):
    vx[i] = vx[i] - vcmx
    vy[i] = vy[i] - vcmy
    vz[i] = vz[i] - vcmz
    S = S + vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]

# On va opérer un rescaling des vitesses initiales à la température T = Tr pour accélerer la phase d'équilibration

ls = math.sqrt(Tr*3*N/S)

for i in range (N):
    vx[i] = vx[i]*ls
    vy[i] = vy[i]*ls
    vz[i] = vz[i]*ls

vcmx = sum(vx)/N
vcmy = sum(vy)/N
vcmz = sum(vz)/N

Ec = 0 # On initialise l'Energie Cinétique Ec

# On va s'assurer que la quantité totale de mouvement reste nulle et calculer l'Energie Cinétique réduite Ec

for i in range (N):
    vx[i] = vx[i] - vcmx
    vy[i] = vy[i] - vcmy
    vz[i] = vz[i] - vcmz
    Ec = Ec + vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]

# On calcul la température instantanée Tinst

Tinst = Ec/(3*N)

print ("La température Instantanée est = ",Tinst)

Ep = 0               # On initialise l'énergie potentielle

Virial = 0           # On initialise le viriel

# On commence les itérations

f = open ("result_themod_NVE.txt", "w+") #Ouverture du fichier "result_themod_NVE.txt" pour l'écriture des données
                                         #pour visualisations

for io in range (Niter):

# Calcul des Forces

# On inititalise les forces (qui rappelons le, dépendent des positions)
    Fx = [] 
    Fy = []
    Fz = []

    for i in range(N):
    
        Fx.append(0)
        Fy.append(0)
        Fz.append(0)

# On commence à calculer les forces
# Rappel : les forces sont centrales : ne dépendent que de rij la distance entre la particule i et la particule j
# Rappel : les forces sont dérivées d'un potentiel radiale (ne dépend que de rij) Lennard Jones Cubic Spline LJCS

    for i in range(N-1):
  
        for j in range(i+1,N):
                     
            xij = x[i] - x[j]
            yij = y[i] - y[j]
            zij = z[i] - z[j]

            # Conditions de l'Image Minimum pour xij

            if (math.fabs(xij) > Ls2) :
                if (xij > 0) :
                   xijn = xij - L
                elif (xij < 0) :
                   xijn = xij + L
            else : 
                   xijn = xij

            # Conditions de l'Image Minimum pour yij

            if (math.fabs(yij) > Ls2) :
                if (yij > 0) :
                   yijn = yij - L
                elif (yij < 0) :
                   yijn = yij + L
            else : 
                   yijn = yij

            # Conditions de l'Image Minimum pour zij

            if (math.fabs(zij) > Ls2) :
                if (zij > 0) :
                   zijn = zij - L
                elif (zij < 0) :
                   zijn = zij + L
            else : 
                   zijn = zij

            # Maintenant que les conditions de l'Image Minimum sont remplies :
            # On calcul rij la distance entre la particule i et la particule j

            rij = math.sqrt(xijn*xijn + yijn*yijn + zijn*zijn)

            # Maintenant on procède au calcul de la dérivée du potentiel LJCS

            r14 = rij**14

            r8  = rij**8

            r12 = rij**12

            r6  = rij**6

            if (rij < rsljs) :
                dF  = 4*((12/r14)-(6/r8))
                Phi = 4*((1/r12) - (1/r6))
            elif (rij < rcljs) :
                dF = -1/rij*(2*Aljs*(rij-rcljs)/rsljs**2 + 3*Bljs*(rij-rcljs)**2/rsljs**3)
                Phi = Aljs * (rij - rcljs) ** 2 / (rsljs ** 2) + Bljs * (rij - rcljs) ** 3 / (rsljs ** 3)
            else :
                dF = 0
                Phi = 0

            # Maintenant on procède au calcul des forces d'après la dérivée du potentiel

            Fxij = dF * xijn 
            Fyij = dF * yijn
            Fzij = dF * zijn

            Fx[i] = Fx[i] + Fxij
            Fy[i] = Fy[i] +  Fyij
            Fz[i] = Fz[i] +  Fzij
         
            Fx[j] = Fx[j] - Fxij
            Fy[j] = Fy[j] - Fyij
            Fz[j] = Fz[j] - Fzij

            Ep = Ep + Phi

            Virial = Virial + dF * rij

# Fin du Calcul des Forces

# Calcul des nouvelles vitesses selon l'Algorithme de Verlet

    for i in range(N):
        vx[i] = vx[i] + dts2*Fx[i]
        vy[i] = vy[i] + dts2*Fy[i]
        vz[i] = vz[i] + dts2*Fz[i]

    vcmx = sum(vx)/N
    vcmy = sum(vy)/N
    vcmz = sum(vz)/N

# On va s'assurer que la quantité totale de mouvement soit nulle

    for i in range(N):
        vx[i] = vx[i] - vcmx
        vy[i] = vy[i] - vcmy
        vz[i] = vz[i] - vcmz

# Calcul des nouvelles positions selon l'Algorithme de Verlet

    for i in range(N):
        
        xn = x[i] + dt*vx[i]

        # Conditions aux Limites Periodiques pour x

        if xn > L :
           x[i] = xn - L
        elif xn < 0 :
           x[i] = xn + L
        else : 
           x[i] = xn
      
        yn = y[i] + dt*vy[i]

        # Conditions aux Limites Periodiques pour y

        if yn > L :
           y[i] = yn - L
        elif yn < 0 :
           y[i] = yn + L
        else : 
           y[i] = yn
        
        zn = z[i] + dt*vz[i]    

        # Conditions aux Limites Periodiques pour z

        if zn > L :
           z[i] = zn - L
        elif zn < 0 :
           z[i] = zn + L
        else : 
           z[i] = zn

# Les positons on changées dont on doit recalculer les Forces

# On re-inititalise les forces (les positions ont changée)
# On inititalise les forces (qui rappelons le, dépendent des positions)

    Fx = []
    Fy = []
    Fz = []

    for i in range(N):
    
        Fx.append(0)
        Fy.append(0)
        Fz.append(0)

    Ep = 0               # On initialise l'énergie potentielle

    Virial = 0           # On initialise le viriel

    for i in range(N-1):
  
        for j in range(i+1,N):
                     
            xij = x[i] - x[j]
            yij = y[i] - y[j]
            zij = z[i] - z[j]

            # Conditions de l'Image Minimum pour xij

            if (math.fabs(xij) > Ls2) :
                if (xij > 0) :
                   xijn = xij - L
                elif (xij < 0) :
                   xijn = xij + L
            else : 
                   xijn = xij

            # Conditions de l'Image Minimum pour yij

            if (math.fabs(yij) > Ls2) :
                if (yij > 0) :
                   yijn = yij - L
                elif (yij < 0) :
                   yijn = yij + L
            else : 
                   yijn = yij

            # Conditions de l'Image Minimum pour zij

            if (math.fabs(zij) > Ls2) :
                if (zij > 0) :
                   zijn = zij - L
                elif (zij < 0) :
                   zijn = zij + L
            else : 
                   zijn = zij

            # Maintenant que les conditions de l'Image Minimum sont remplies :
            # On calcul rij la distance entre la particule i et la particule j

            rij = math.sqrt(xijn*xijn + yijn*yijn + zijn*zijn)    

            # Maintenant on procède au calcul de la dérivée du potentiel LJCS

            r14 = rij ** 14

            r12 = rij ** 12

            r8 = rij ** 8

            r6 = rij ** 6

            if (rij < rsljs):
                dF  = 4 * ((12 / r14) - (6 / r8))
                Phi = 4 * (1  / r12 - 1 / r6)
            elif (rij < rcljs):
                dF = -1 / rij * (2 * Aljs * (rij - rcljs) / rsljs ** 2 + 3 * Bljs * (rij - rcljs) ** 2 / rsljs ** 3)
                Phi = Aljs * (rij - rcljs) ** 2 / (rsljs ** 2) + Bljs * (rij - rcljs) ** 3 / (rsljs ** 3)
            else:
                dF = 0
                Phi = 0

            # Maintenant on procède au calcul des forces d'après la dérivée du potentiel

            Fxij = dF * xijn
            Fyij = dF * yijn
            Fzij = dF * zijn

            Fx[i] = Fx[i] + Fxij
            Fy[i] = Fy[i] + Fyij
            Fz[i] = Fz[i] + Fzij

            Fx[j] = Fx[j] - Fxij
            Fy[j] = Fy[j] - Fyij
            Fz[j] = Fz[j] - Fzij

            Ep = Ep + Phi

            Virial = Virial + dF * rij

# Fin du Calcul des Forces

   # Calcul des nouvelles vitesses selon l'Algorithme de Verlet

    for i in range(N):

        vx[i] = vx[i] + dts2*Fx[i]
        vy[i] = vy[i] + dts2*Fy[i]
        vz[i] = vz[i] + dts2*Fz[i]

    vcmx = sum(vx)/float(N)
    vcmy = sum(vy)/float(N)
    vcmz = sum(vz)/float(N)
   
    Ec = 0 # On initialise l'Energie Cinétique Ec

    # On va s'assurer que la quantité totale de mouvement reste nulle et calculer l'Energie Cinétique réduite Ec

    for i in range(N):

        vx[i] = vx[i] - vcmx
        vy[i] = vy[i] - vcmy
        vz[i] = vz[i] - vcmz
        
        Ec = Ec + vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]

# On calcul la température instantanée Tinst

    Tinst = Ec/(3*N)

# On calcul la pression instantanée Pinst
    
    Pinst = Ro*(Ec + Virial)/(3*N)

# On calcul l'énergie totale

    ET = Ec + Ep

# On affiche  les résultats

    print (io,Tinst,Pinst,ET)

    f.write("%s %s %s %s\n" % (io,Tinst,Pinst,ET))

f.close()
