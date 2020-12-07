# Programme de simulation informatique produisant différents graphiques
# Groupe 11.27
# Thibault, Clara, Nicolas, Romain, Hugolin et Safiya
 
import math
import matplotlib.pyplot as plt
import numpy as np
 

### Constantes
 
g = 9.81             # Constante de gravitation [m/s^2] 
 
### Paramètres du système
 
D = 10                  # coef d'amortissement
I = 6.596               # Moment d'inertie [kg.m^2]
mtot = 5                # Masse totale [kg]
m = .4                  # Masse déplacée [kg]
m1 = 3                  # Masse de la barge [kg]
m2 = 1.5                # Masse de la grue [kg]
d = 0.85                # Distance [m]
L = 0.6                 # Largeur de la base carrée [m]
h1 = 0.15               # Hauteur de la base [m]
h2 = 0.5                # Hauteur du chargement [m]
h3 = 0.1                # Hauteur de la masse portée [m]
hc = mtot/(L**2*997)    # 997 étant la masse volumique de l'eau [m]
 
### Paramètres de la simulation
 
step = 0.01     # Pas (dt) [s]
end = 20        # Durée [s]
x_0 = 0         # Angle initial [rad]
x1_0 = 0        # Vitesse angulaire initiale [rad/s]
x2_0 = 0        # Accélération angulaire initiale [rad/s**2]
sum_c_0 = 0     # Somme des couples initiaux [N.m]
Ca_0 = -m*g*d   # Couple initiale de chavirement [N.m]
Cr_0 = 0        # Couple de redressement initial [N.m]
dist_0 = 0      # Distance entre le centre de gravité et le centre de poussée selon X initiale [m]
 
### Matrices de stockage des données
 
t = np.arange(0, end, step)  # Durée [s]
x = np.empty_like(t)         # Angle [rad]
x1 = np.empty_like(t)        # Vitesse angulaire [rad/s]
x2 = np.empty_like(t)        # Accélération angulaire [rad/s**2]
sum_c = np.empty_like(t)     # Somme des couples [N.m]
Ca = np.empty_like(t)        # Couple de chavirement [N.m]
Cr = np.empty_like(t)        # Couple de redressement [N.m]
dist = np.empty_like(t)      # Distance entre la barge et la charge [m]
xc = np.empty_like(t)        # Centre de poussée, coordonnée du x [m]
xg = np.empty_like(t)        # Centre de gravité, coordonnée du x [m]
distance = np.empty_like(t)  # Distance en xc et xg [m]
 
 
def simulation(m):
    
    # Equation: I*ang'' + D*ang' = C
 
    ### Paramètres initiaux de l'équation
    
    x[0] = x_0      
    x1[0] = x1_0    
    x2[0] = x2_0    
    sum_c[0] = sum_c_0
    Ca[0] = Ca_0
    Cr[0] = Cr_0
    dist[0] = dist_0
 
    dt = step
    P = mtot*g
    Pa = m*g
    
    # calcul pas-à-pas de l'EDO de second ordre. On stocke toutes les données dans chaque matrice respective.
    
    for i in range(len(t)-1):
        
        dist[i+1] = d*math.cos(x[i])-h3*math.sin(x[i])    # Distance 
        Ca[i+1] = -m*g*dist[i]                            # Couple de chavirement
 
        xc[i+1] = ((12*hc**2-L**2)/(24*hc))*math.sin(x[i]) - (L**2/(24*hc*math.cos(x[i])))*math.tan(x[i])   # Centre de poussée, coordonnée en x
        xg[i+1] = ((m1*(h1/2-hc)+m2*(h2+h1-hc))/(mtot))*math.sin(x[i])                                      # Centre de gravité, coordonnée en x                                   
 
        distance[i+1] = xc[i] - xg[i]      # Distance entre xc et xg
 
        Cr[i+1] = mtot*g*abs(distance[i])  # Couple de redressement
        sum_c[i+1] = Ca[i] + Cr[i]         # Somme des couples
 
 
        x2[i] = (sum_c[i]-D*x1[i])/I    # Accélération angulaire
        x1[i+1] = x1[i] + x2[i] * dt    # Vitesse angulaire
        x[i+1] = x[i] - x1[i] * dt      # Angle
        x2[i+1] = x2[i]
 
 
### Paramètres du graphique
        
angle_max = math.atan(2*(h1-hc)/L)    # Angle max
angle_s = math.atan(2*hc/L)           # Angle de soulèvement
 
### Variable permettant de passer de radians en degrés
 
rad_to_deg = 57.29578
r = rad_to_deg
 
 
### Graphique des différentes masses [kg] en fonction de la distance [m] à laquelle elles sont déplacées
 
def graphiques_distances_masses():
    
    simulation(m)
    
    # Traçage de chaque fonction pour chaque masse 
    
    simulation(1)
    x_values = [0, 2]
    y_values = [x[0]*r, x[20]*r]
    plt.plot(x_values,y_values,label='m = 1kg')
    plt.legend()
    
    simulation(2)
    y_values2 = [x[0]*r, x[20]*r]
    plt.plot(x_values,y_values2,label='m = 2kg')
    plt.legend()
    
    simulation(3)
    y_values3 = [x[0]*r, x[20]*r]
    plt.plot(x_values,y_values3,label='m = 3kg')
    plt.legend()
    simulation(4)
    y_values4 = [x[0]*r, x[20]*r]
    plt.plot(x_values,y_values4,label='m = 4kg')
    plt.legend()
    
    simulation(5)
    y_values5 = [x[0]*r, x[20]*r]
    plt.plot(x_values,y_values5,label='m = 5kg')
    plt.legend()
    
    # Paramètres généraux du graphique
    
    plt.plot([0,end], [angle_max*r,angle_max*r], '--r', label='submersion') # Ligne en pointillé qui marque le submersion (positif)
    plt.plot([0,end], [-angle_max*r,-angle_max*r], '--r')                   # Ligne en pointillé qui marque le submersion (négatif)
    plt.plot([0,end], [angle_s*r,angle_s*r], '--b', label='soulèvement')    # Ligne en pointillé qui marque le soulèvement (positif)
    plt.plot([0,end], [-angle_s*r,-angle_s*r], '--b')                       # Ligne en pointillé qui marque le soulèvement (négatif)
 
    plt.xlabel('Distance (m)')
    plt.ylabel('Angle (º)')
    plt.xlim(0,2)
    plt.ylim(0,angle_max*r+1)
    plt.show()
 
 
### Graphique de l'angle [degrés] en fonction du temps [s]
    
def graphique_angle_deg():
    
    # Paramètres généraux du graphique
    
    simulation(m)
    plt.plot(t,x*r)
    plt.xlim(0,end)
    plt.xlabel('Temps (s)')
    plt.ylim(0,angle_max*r+1)
    plt.ylabel('Angle (º)')
 
    plt.plot([0,end], [angle_max*r,angle_max*r], '--r', label='submersion') # Ligne en pointillé qui marque le submersion (positif)
    plt.plot([0,end], [-angle_max*r,-angle_max*r], '--r')                   # Ligne en pointillé qui marque le submersion (négatif)
    plt.plot([0,end], [angle_s*r,angle_s*r], '--b', label='soulèvement')    # Ligne en pointillé qui marque le soulèvement (positif)
    plt.plot([0,end], [-angle_s*r,-angle_s*r], '--b')                       # Ligne en pointillé qui marque le soulèvement (négatif)
 
    plt.show()
 
 
### 3 graphique :
#                 1. Angle [rad] en fonction du temps [s]
#                 2. Vitesse angulaire [rad/s] en fonction du temps [s]
#                 3. Accélération angulaire[rad/s**2] en fonction du temps [s]
    
def graphiques():
 
    simulation(m)
    
    # Sous-graphique : 1. Angle [rad] en fonction du temps [s]
    
    plt.figure(1)
 
    plt.subplot(3,1,1)
    plt.plot(t,x, label="Angle")
    plt.legend()
    plt.xlim(0,end)
    
    plt.plot([0,end], [angle_max,angle_max], '--r', label='submersion') # Ligne en pointillé qui marque le submersion (positif)
    plt.plot([0,end], [-angle_max,-angle_max], '--r')                   # Ligne en pointillé qui marque le submersion (négatif)
    plt.plot([0,end], [angle_s,angle_s], '--b', label='soulèvement')    # Ligne en pointillé qui marque le soulèvement (positif)
    plt.plot([0,end], [-angle_s,-angle_s], '--b')                       # Ligne en pointillé qui marque le soulèvement (négatif)
    plt.ylabel('Angle \n (rad)',fontsize = 5)
    
 
    # Sous-graphique : 2. Vitesse angulaire [rad/s] en fonction du temps [s]
    
    plt.subplot(3,1,2)
    plt.plot(t,x1,'g', label="Vitesse angulaire")
    plt.legend()
    plt.xlim(0,end)
    plt.ylabel('Vitesse angulaire \n (rad/s)',fontsize = 5)
    
    # Sous-graphique : 3. Accélération angulaire[rad/s**2] en fonction du temps [s]
    
    plt.subplot(3,1,3)
    plt.plot(t,x2,'#e3902b', label="Accélération angulaire")
    plt.legend()
    plt.xlim(0,end)
    plt.ylabel('Accélération angulaire \n (rad/s^2)',fontsize = 5)
    
    # Paramètres généraux

    plt.xlabel('Temps (s)')
    plt.show()
    
 
### Diagramme de phase : graphique représentant la vitesse angulaire [rad/s] en fonction de l'angle [rad]
 
def diagramme_de_phase():
    
    # Paramètres généraux du graphique
    
    simulation(m)
    plt.title("Diagramme de phase")
    plt.plot(x, x1,linewidth=0.75)
    plt.xlabel("Angle (rad)")
    plt.ylabel("Vitesse angulaire (rad/s)")
    plt.show()
    
 
### Programme principal
 
graphiques_distances_masses()
graphique_angle_deg()
graphiques()
diagramme_de_phase()
