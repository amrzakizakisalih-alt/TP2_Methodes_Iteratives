import numpy as np
from random import random 





# Exemple deformation elastique :
#  - Entree : listeNoeuds (liste coordonnees x,y,z), listeArcs (liste indices), listeMasse (liste double), listeFication (liste entier)
#  - Sortie : tD (matrice), tK (matrice)
def generateMatricePont(listeNoeuds, listeArcs, listeMasse, listeFixation):
    n = len(listeNoeuds)
    K = np.zeros((3*n,3*n))
    D = np.zeros((3*n,3*n))
    F = np.zeros(3*n)
    
    matU = np.zeros((3,3))
    
    
    # Premiere partie : Construction de A
    for aa in listeArcs:
        # Coordonnées des points de l'arc
        p1 = int(aa[0])
        p2 = int(aa[1])
        # Valeur de la raideur
        kk = aa[2]
        
        u = np.array(listeNoeuds[p1]) - np.array(listeNoeuds[p2])
        u = u / np.sqrt(np.dot(u.transpose(),u))
        # Mat (u u^t) :
        for i in range(3):
            for j in range(3):
                matU[i,j] = u[i]*u[j]
        
        indp1 = p1*3
        indp2 = p2*3
        # Completer ICI :
        K[indp1:indp1+3, indp1:indp1+3] += kk * matU  # Contribution pour p1
        K[indp1:indp1+3, indp2:indp2+3] -= kk * matU  # Contribution entre p1 et p2
        K[indp2:indp2+3, indp1:indp1+3] -= kk * matU  # Contribution entre p2 et p1
        K[indp2:indp2+3, indp2:indp2+3] += kk * matU # Contribution pour p2
        
    for j,m in enumerate(listeMasse):
        D[3*j,3*j] = m
        D[3*j+1,3*j+1] = m
        D[3*j+2,3*j+2] = m
        
        F[3*j+2] = -0.1
    
    F = np.dot(D,F)
    
    # Deuxième partie : Prise en compte des fixations    
    listePts = []
    for i in range(n):
        fixe = False
        for e in listeFixation:
            if (i==int(e)):
                fixe = True
                break
        if(fixe == False):
            listePts.append(i)
    
    tn = len(listePts)
    tK = np.zeros((3*tn,3*tn))
    tD = np.zeros((3*tn,3*tn))
    tF = np.zeros(3*tn)                    
        
    for j,e in enumerate(listePts):
        for k,f in enumerate(listePts):
            tK[3*j,3*k] = K[3*e,3*f]
            tK[3*j,3*k+1] = K[3*e,3*f+1]
            tK[3*j,3*k+2] = K[3*e,3*f+2]
            tK[3*j+1,3*k] = K[3*e+1,3*f]
            tK[3*j+1,3*k+1] = K[3*e+1,3*f+1]
            tK[3*j+1,3*k+2] = K[3*e+1,3*f+2]
            tK[3*j+2,3*k] = K[3*e+2,3*f]
            tK[3*j+2,3*k+1] = K[3*e+2,3*f+1]
            tK[3*j+2,3*k+2] = K[3*e+2,3*f+2]
            
            tD[3*j,3*k] = D[3*e,3*f]
            tD[3*j,3*k+1] = D[3*e,3*f+1]
            tD[3*j,3*k+2] = D[3*e,3*f+2]
            tD[3*j+1,3*k] = D[3*e+1,3*f]
            tD[3*j+1,3*k+1] = D[3*e+1,3*f+1]
            tD[3*j+1,3*k+2] = D[3*e+1,3*f+2]
            tD[3*j+2,3*k] = D[3*e+2,3*f]
            tD[3*j+2,3*k+1] = D[3*e+2,3*f+1]
            tD[3*j+2,3*k+2] = D[3*e+2,3*f+2]
        
        tF[3*j] = F[3*e]
        tF[3*j+1] = F[3*e]
        tF[3*j+2] = F[3*e+2]     
       
    return [tD,tK]
    
def generateSourcePont(listeNoeuds, listeArcs, listeMasse, listeFixation,typeSource):
    n = len(listeNoeuds)
    F = np.zeros(3*n)    
    
    for j,m in enumerate(listeMasse):
        if(typeSource=='Pression'):
            F[3*j+2] = -0.002
        if(typeSource=='Cisaillement'):
            F[3*j+2] = -0.002*(listeNoeuds[j][1]-1.0)
            #F[3*j] = -0.1    
    
    # Deuxième partie : Prise en compte des fixations    
    listePts = []
    for i in range(n):
        fixe = False
        for e in listeFixation:
            if (i==int(e)):
                fixe = True
                break
        if(fixe == False):
            listePts.append(i)
    
    tF = np.zeros(len(listePts)*3)
    
    for j,e in enumerate(listePts):
        tF[3*j] = F[3*e]
        tF[3*j+1] = F[3*e]
        tF[3*j+2] = F[3*e+2]     
       
    return tF