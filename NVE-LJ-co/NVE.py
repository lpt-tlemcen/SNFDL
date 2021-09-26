import math

read_input = open ("NVE_input_data.txt", "r") # Ouverture du fichier "NVE_input_data.txt" pour la lecture des
                                              # données de simulation

with open('NVE_input_data.txt') as f:
     for line in f:
          a = line.split()
          N = int(a[0])  # Nombre de particules
          print("N =", N)
          Ro = float(a[1])  # Densité réduite
          print("Ro =", Ro)
          Tr = float(a[2])  # Température réduite
          print("Tr =", Tr)
          Niter = int(a[3])  # Nombre d'itération
          print("Niter =", Niter)
          dt = float(a[4])  # Delta(t) Pas de temps en unité réduite
          print("dt =", dt)
          Niter_equilibrium = int(a[5]) # Nombre d'itération pour équilibrer le système
          print("Niter_equilibrium =", Niter_equilibrium)
          Niter_cor = int(a[6]) # Nombre d'itération inclut dans le calcul du RMS
          print("Niter_cor =", Niter_cor)
          Conf_Type = int(a[7]) # Nombre d'itération inclut dans le calcul du RMS
          print("Conf_Type =", Conf_Type)


ncfc = round((N / 4.) ** (1 / 3))  # Nombre de petites boîtes CFC
L = (N / Ro) ** (1 / 3)  # Longueur réduite
Ls2 = L / 2.  # Moitié de la longueur de la boîtes

print("ncfc = ", ncfc)
print("L = ", L)

# Paramètre du Potentiel Lennard Jones Cubic Spline LJCS , Réf. :
# Bjørn Hafskjold, Karl Patrick Travis, Amanda Bailey Hass, Morten Hammer, Ailo Aasen & Øivind Wilhelmsen (2019)
# Thermodynamic properties of the 3D Lennard-Jones/spline model, Molecular Physics, DOI: 10.1080/00268976.2019.1664780

VLJtr  = - 0.0163    # valeur du potentiel de L_J pour rc=2.5
rcLJtr = 2.5       #rayon de trancature potentiel LJ

dts2 = dt / 2  # Delta(t)/2 Moitié du pas de temps en unité réduite

Delta_r = 0.01
Maxbin = round((Ls2 - Delta_r) / Delta_r)
number_pi = 3.14659824235543546
gr_const = (4. / 3.) * number_pi / (N / Ro)

if (Conf_Type  == 1):
    # pour la premier boite (initialisation)

    # itinialisation des tableaux des positions

    x = []
    y = []
    z = []

    # les positions initiales :

    # on va déplacer notre petite boite CFC avec ses 04 atomes dans les 03 axes X,Y,Z

    for i in range(ncfc):  # déplacement de la petite boite CFC selon X
        for j in range(ncfc):  # déplacement de la petite boite CFC selon Y
            for k in range(ncfc):  # déplacement de la petite boite CFC selon Z
                x.append(0 + i)  # x atome 1 de la petite boite CFC
                y.append(0 + j)  # y atome 1 de la petite boite CFC
                z.append(0 + k)  # z atome 1 de la petite boite CFC
                x.append(1 / 2 + i)  # x atome 2 de la petite boite CFC
                y.append(1 / 2 + j)  # y atome 2 de la petite boite CFC
                z.append(0 + k)  # z atome 2 de la petite boite CFC
                x.append(1 / 2 + i)  # x atome 3 de la petite boite CFC
                y.append(0 + j)  # y atome 3 de la petite boite CFC
                z.append(1 / 2 + k)  # z atome 3 de la petite boite CFC
                x.append(0 + i)  # x atome 4 de la petite boite CFC
                y.append(1 / 2 + j)  # y atome 4 de la petite boite CFC
                z.append(1 / 2 + k)  # z atome 4 de la petite boite CFC

    # Un rescaling spatial (*L/ncfc) pour s'assurer que ce soient les coordonnées "réelles" des particules à l'intérieur de la boîte cubique de longeur L

    for i in range(N):
        x[i] = x[i] * L / ncfc
        y[i] = y[i] * L / ncfc
        z[i] = z[i] * L / ncfc

    # les vitesses initiales :

    # appeler le générateur de nombres pseudo-aléatoires
    from random import seed
    from random import random

    # initier le  générateur de nombre pseudo-aléatoire
    seed(1)

    # itinialisation des tableaux des positions

    vx = []
    vy = []
    vz = []

    # on va générer les vitesses de nos particules extraite d'une distribution de Maxwell-Boltzmann à T = Tr pour les 03 composantes vx, vy & vz

    for i in range(N):
        Ax = random()  # nombre pseudo aléatoire
        vxi = math.sqrt(
            -2 * Tr * math.log(Ax))  # composante x de la vitesse extraite de la distribution de Maxwell-Boltzmann
        vx.append(vxi)
        Ay = random()  # nombre pseudo aléatoire
        vyi = math.sqrt(
            -2 * Tr * math.log(Ay))  # composante y de la vitesse extraite de la distribution de Maxwell-Boltzmann
        vy.append(vyi)
        Az = random()  # nombre pseudo aléatoire
        vzi = math.sqrt(
            -2 * Tr * math.log(Az))  # composante z de la vitesse extraite de la distribution de Maxwell-Boltzmann
        vz.append(vzi)

    # On va s'assurer que la quantité totale de mouvement soit nulle

    vcmx = sum(vx) / N
    vcmy = sum(vy) / N
    vcmz = sum(vz) / N

    S = 0

    for i in range(N):
        vx[i] = vx[i] - vcmx
        vy[i] = vy[i] - vcmy
        vz[i] = vz[i] - vcmz
        S = S + vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]

    # On va opérer un rescaling des vitesses initiales à la température T = Tr pour accélerer la phase d'équilibration

    ls = math.sqrt(Tr * 3 * N / S)

    for i in range(N):
        vx[i] = vx[i] * ls
        vy[i] = vy[i] * ls
        vz[i] = vz[i] * ls

    vcmx = sum(vx) / N
    vcmy = sum(vy) / N
    vcmz = sum(vz) / N

else :

    # Lecture de la configuration position et vitesse

    x = []
    y = []
    z = []

    vx = []
    vy = []
    vz = []

    fconf = open("Config_output_NVE.txt", "r")  # Ouverture du fichier "result_config_NVE_org" pour la lecture de la
    # configuration de sortie coordonnées et vitesses

    with open('Config_output_NVE.txt') as f:
        for line in f:
            a = line.split()
            indice = int(a[0])
            xpart = float(a[1])
            ypart = float(a[2])
            zpart = float(a[3])
            vxpart = float(a[4])
            vypart = float(a[5])
            vzpart = float(a[6])
            x.append(xpart)
            y.append(ypart)
            z.append(zpart)
            vx.append(vxpart)
            vy.append(vypart)
            vz.append(vzpart)


Ec = 0  # On initialise l'Energie Cinétique Ec

# On va s'assurer que la quantité totale de mouvement reste nulle et calculer l'Energie Cinétique réduite Ec

vcmx = sum(vx)/N
vcmy = sum(vy)/N
vcmz = sum(vz)/N

for i in range(N):
    vx[i] = vx[i] - vcmx
    vy[i] = vy[i] - vcmy
    vz[i] = vz[i] - vcmz
    Ec = Ec + vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]

# On calcul la température instantanée Tinst

Tinst = Ec / (3 * N)

print("La température Instantanée est = ", Tinst)

Ep = 0  # On initialise l'énergie potentielle

Virial = 0  # On initialise le viriel

hist = []

for ib in range(10000):
    hist.append(0)

RMS = []
VAF = []

for ic in range(Niter_cor):
    RMS.append(0)
    VAF.append(0)

x_sclp = []
y_sclp = []
z_sclp = []

for ib in range(N):
    x_sclp.append(0)
    y_sclp.append(0)
    z_sclp.append(0)

xx = []
yy = []
zz = []

xx = [[0 for j in range(0,N)] for i in range(0,Niter_cor)]
yy = [[0 for j in range(0,N)] for i in range(0,Niter_cor)]
zz = [[0 for j in range(0,N)] for i in range(0,Niter_cor)]

vxx = []
vyy = []
vzz = []

vxx = [[0 for j in range(0,N)] for i in range(0,Niter_cor)]
vyy = [[0 for j in range(0,N)] for i in range(0,Niter_cor)]
vzz = [[0 for j in range(0,N)] for i in range(0,Niter_cor)]

# On commence les itérations

f = open ("result_themod_NVE.txt", "w+") #Ouverture du fichier "result_themod_NVE.txt" pour l'écriture des données
                                         # pour visualisations

frms = open("result_RMS_NVE.txt", "w+")  # Ouverture du fichier "result_RMS_NVE.txt" pour l'écriture des données
                                         #pour visualisations

fvaf = open("result_VAF_NVE.txt", "w+")  # Ouverture du fichier "result_VAF_NVE.txt" pour l'écriture des données
                                         #pour visualisations

N_indic_bloc = 1

# On commence à itérer

for io in range(Niter):

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

    for i in range(N - 1):

        for j in range(i + 1, N):

            xij = x[i] - x[j]
            yij = y[i] - y[j]
            zij = z[i] - z[j]

            # Conditions de l'Image Minimum pour xij

            if (math.fabs(xij) > Ls2):
                if (xij > 0):
                    xijn = xij - L
                elif (xij < 0):
                    xijn = xij + L
            else:
                xijn = xij

            # Conditions de l'Image Minimum pour yij

            if (math.fabs(yij) > Ls2):
                if (yij > 0):
                    yijn = yij - L
                elif (yij < 0):
                    yijn = yij + L
            else:
                yijn = yij

            # Conditions de l'Image Minimum pour zij

            if (math.fabs(zij) > Ls2):
                if (zij > 0):
                    zijn = zij - L
                elif (zij < 0):
                    zijn = zij + L
            else:
                zijn = zij

            # Maintenant que les conditions de l'Image Minimum sont remplies :
            # On calcul rij la distance entre la particule i et la particule j

            rij = math.sqrt(xijn * xijn + yijn * yijn + zijn * zijn)

            ibin = round(rij / Delta_r) + 1

            if (ibin <= Maxbin):
                hist[ibin] = hist[ibin] + 2
            else:
                hist[ibin] = hist[ibin] + 2

            # Maintenant on procède au calcul de la dérivée du potentiel LJCS

            r14 = rij ** 14

            r8 = rij ** 8

            r12 = rij ** 12

            r6 = rij ** 6

            if (rij < rcLJtr):
                dF = 4 * ((12 / r14) - (6 / r8))
                Phi = 4 * ((1 / r12) - (1 / r6)) - VLJtr

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

    for i in range(N):  # iL1 propagateur partiel
        vx[i] = vx[i] + dts2 * Fx[i]
        vy[i] = vy[i] + dts2 * Fy[i]
        vz[i] = vz[i] + dts2 * Fz[i]

    vcmx = sum(vx) / N
    vcmy = sum(vy) / N
    vcmz = sum(vz) / N

    # On va s'assurer que la quantité totale de mouvement soit nulle

    for i in range(N):
        vx[i] = vx[i] - vcmx
        vy[i] = vy[i] - vcmy
        vz[i] = vz[i] - vcmz

    for i in range(N):  # iL2 propagateur partiel (central)

        xn = x[i] + dt * vx[i]

        x_sclp[i] = x_sclp[i] + ( xn - x[i] )

        # Conditions aux Limites Periodiques pour x

        if xn > L:
            x[i] = xn - L
        elif xn < 0:
            x[i] = xn + L
        else:
            x[i] = xn

        yn = y[i] + dt * vy[i]

        y_sclp[i] = y_sclp[i] + ( yn - y[i] )

        # Conditions aux Limites Periodiques pour y

        if yn > L:
            y[i] = yn - L
        elif yn < 0:
            y[i] = yn + L
        else:
            y[i] = yn

        zn = z[i] + dt * vz[i]

        z_sclp[i] = z_sclp[i] + ( zn - z[i] )

        # Conditions aux Limites Periodiques pour z

        if zn > L:
            z[i] = zn - L
        elif zn < 0:
            z[i] = zn + L
        else:
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

    Ep = 0  # On initialise l'énergie potentielle

    Virial = 0  # On initialise le viriel

    for i in range(N - 1):

        for j in range(i + 1, N):

            xij = x[i] - x[j]
            yij = y[i] - y[j]
            zij = z[i] - z[j]

            # Conditions de l'Image Minimum pour xij

            if (math.fabs(xij) > Ls2):
                if (xij > 0):
                    xijn = xij - L
                elif (xij < 0):
                    xijn = xij + L
            else:
                xijn = xij

            # Conditions de l'Image Minimum pour yij

            if (math.fabs(yij) > Ls2):
                if (yij > 0):
                    yijn = yij - L
                elif (yij < 0):
                    yijn = yij + L
            else:
                yijn = yij

            # Conditions de l'Image Minimum pour zij

            if (math.fabs(zij) > Ls2):
                if (zij > 0):
                    zijn = zij - L
                elif (zij < 0):
                    zijn = zij + L
            else:
                zijn = zij

            # Maintenant que les conditions de l'Image Minimum sont remplies :
            # On calcul rij la distance entre la particule i et la particule j

            rij = math.sqrt(xijn * xijn + yijn * yijn + zijn * zijn)

            ibin = round(rij / Delta_r) + 1

            if (ibin <= Maxbin):
                hist[ibin] = hist[ibin] + 2
            else:
                hist[ibin] = hist[ibin] + 2

            # Maintenant on procède au calcul de la dérivée du potentiel LJCS

            r14 = rij ** 14

            r8 = rij ** 8

            r12 = rij ** 12

            r6 = rij ** 6

            if (rij < rcLJtr):
                dF = 4 * ((12 / r14) - (6 / r8))
                Phi = 4 * ((1 / r12) - (1 / r6)) - VLJtr

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

    for i in range(N):  # iL1 propagateur partiel
        vx[i] = vx[i] + dts2 * Fx[i]
        vy[i] = vy[i] + dts2 * Fy[i]
        vz[i] = vz[i] + dts2 * Fz[i]

    vcmx = sum(vx) / N
    vcmy = sum(vy) / N
    vcmz = sum(vz) / N

    # On va s'assurer que la quantité totale de mouvement soit nulle

    for i in range(N):
        vx[i] = vx[i] - vcmx
        vy[i] = vy[i] - vcmy
        vz[i] = vz[i] - vcmz

    Ec = 0  # On initialise l'Energie Cinétique Ec

    for i in range(N):
        Ec = Ec + vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]

    # On calcul la température instantanée Tinst

    Tinst = Ec / (3 * N)

    # On calcul la pression instantanée Pinst

    Pinst = Ro * (Ec + Virial) / (3 * N)

    # On calcul l'énergie totale

    ET = Ec + Ep

    # On affiche  les résultats

    print(io, Tinst, Pinst, ET/N)

    if (io >= Niter_equilibrium):

        f.write("%s %20.10f %s20.10f %20.10f\n" % (io, Tinst, Pinst, ET/N))

        N_indic = (io - Niter_equilibrium) % Niter_cor

        for ipc in range(0, N):
            xx[N_indic][ipc] = x_sclp[ipc]
            yy[N_indic][ipc] = y_sclp[ipc]
            zz[N_indic][ipc] = z_sclp[ipc]
            vxx[N_indic][ipc] = vx[ipc]
            vyy[N_indic][ipc] = vy[ipc]
            vzz[N_indic][ipc] = vz[ipc]

        if (N_indic == (Niter_cor - 1)):
            for i in range(Niter_cor):
                for j in range(N):
                    RMS[i] = RMS[i] + (xx[i][j] - xx[0][j]) ** 2 + (yy[i][j] - yy[0][j]) ** 2 + (
                                zz[i][j] - zz[0][j]) ** 2
                    VAF[i] = VAF[i] + vxx[i][j]*vxx[0][j]  + vyy[i][j]*vyy[0][j] + vzz[i][j]*vzz[0][j]

                RMS[i] = RMS[i] / N
                VAF[i] = VAF[i] / N
                frms.write("%3d %20.10f %20.10f \n" % (N_indic_bloc, i * dt, RMS[i]))
                fvaf.write("%3d %20.10f %20.10f \n" % (N_indic_bloc, i * dt, VAF[i]))
            N_indic_bloc = N_indic_bloc + 1
            # print("N_indic_bloc",N_indic_bloc)

f.close()

frms.close()

fvaf.close()

fgr = open("result_Gr_NVE.txt", "w+")  # Ouverture du fichier "result_Gr_NVE.txt" pour l'écriture des données
                                       # pour visualisations

for ib in range(Maxbin):
    r_lower = (ib - 1) * Delta_r
    r_uper = r_lower + Delta_r
    rr = r_uper ** 3. - r_lower ** 3.
    En_ideal = gr_const * rr
    Gr = hist[ib] / (Niter * 1.) / (N ** 2.) / En_ideal
    fgr.write("%s %20.10f \n" % (ib, Gr))

fgr.close()

fconf = open ("Config_output_NVE.txt", "w+") #Ouverture du fichier "Config_output_NVE.txt" pour l'écriture de la
                                             # configuration de sortie coordonnées et vitesses
for i in range(N):
    fconf.write("%s %s %s %s %s %s %s \n" % (i,x[i],y[i],z[i],vx[i],vy[i],vz[i]))

fconf.close()
