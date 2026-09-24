#!/usr/bin/python3
#-*- coding: utf-8 -*-

import numpy as np
import affichage

import algo
import data

import random

"""
    Etude d'un pont (déformation et modes) 
"""

## Attention : Pour executer le programme, il faut utiliser python3 main.py

##----------------------------------
# But du programme :
# 1. Visualisation pont.
# 2. Calcul de la déformation d'un point.
# 3. Calcul des modes propores.
butPrgm = "1"
print("Que souhaitez vous faire ?")
print(" 1. Visualiser la structure du pont.")
print(" 2. Calculer la déformation du pont.")
print(" 3. Calculer les modes.")
print(" 4. Analyser l'effet de la relaxation sur les méthodes itératives.")
butPrgm = input("Entrer le choix : ")

print("Debut")
##----------------------------------

##----------------------------------
nomData = ""
choixData = "1"
print("Choix donnees ?")
print(" 1. Pont simple")
print(" 2. Pont Sydney.")
print(" 3. Pont en Treillis")
choixData = input("Entrer le choix : ")


listeNoeuds = []
listeArcs = []
listeMasse = []
listeFixation = []

if (choixData=="1"):
    listeNoeuds = np.genfromtxt("PontN.data", delimiter="\t")
    listeArcs = np.genfromtxt("PontA.data", delimiter="\t")
    listeMasse = np.genfromtxt("PontM.data", delimiter="\t")
    listeFixation = np.genfromtxt("PontF.data", delimiter="\t")
    
if (choixData=="2"):
    # A vour de créer les fichiers de données.
    listeNoeuds = np.genfromtxt("PontN_2.data", delimiter="\t")
    listeArcs = np.genfromtxt("PontA_2.data", delimiter="\t")
    listeMasse = np.genfromtxt("PontM_2.data", delimiter="\t")
    listeFixation = np.genfromtxt("PontF_2.data", delimiter="\t")
if (choixData=="3"):
   # A vour de créer les fichiers de données.
   listeNoeuds = np.genfromtxt("PontN_3.data", delimiter="\t")
   listeArcs = np.genfromtxt("PontA_3.data", delimiter="\t")
   listeMasse = np.genfromtxt("PontM_3.data", delimiter="\t")
   listeFixation = np.genfromtxt("PontF_3.data", delimiter="\t")    

print("")
##----------------------------------


##----------------------------------
# Visualiser le pont :
if (butPrgm=="1"):
    affichage.afficherPont(listeNoeuds, listeArcs)
    
    [D,K] = data.generateMatricePont(listeNoeuds, listeArcs, listeMasse, listeFixation)
    is_symmetric = np.allclose(K, K.T)
    print("Symétrie de K :", is_symmetric)
    eigvals = np.linalg.eigvals(K)
    print("Valeurs propres de K :", eigvals)
    print("K est semi-définie positive :", np.all(eigvals >= -1e-12))
    cond = np.linalg.cond(K)
    print("Conditionnement de K :", cond)
##----------------------------------    
    
    
##----------------------------------
# Test dessin :
if (butPrgm=="2"):
    [D,K] = data.generateMatricePont(listeNoeuds, listeArcs, listeMasse, listeFixation)
    typeSource = 'Pression' # 'Pression','Cisaillement' 
    print("Quelle type de source souhaitez-vous utiliser ?")
    print("1. Pression uniforme")
    print("2. Cisaillement")
    typeforce = input("Enter le choix : ")
    if(typeforce=='2'):
        typeSource = 'Cisaillement'
    
    F = data.generateSourcePont(listeNoeuds, listeArcs, listeMasse, listeFixation,typeSource)
    
    
    # Nombre de modes calculé :
    meth = 1
    print("Quelle méthode de résolution souhaitez-vous utiliser ?")
    print("1. Jacobi")
    print("2. Gaus-Seidel")
    print("3. Gradient")
    print("4. GC")
    print("5. Comparaison méthodes")
    meth = input("Entrer le choix : ")
    
    
    tol = 1e-5    
    iteMax = 1e3
    relax = 0.6
    sol = F*0.
    residus = []
    # Resolution du probleme 
    if(meth == '1'):
        sol,residus = algo.algoSplit(K,F*0.,F,tol,iteMax,relax,'J')
    if(meth == '2'):
        sol,residus = algo.algoSplit(K,F*0.,F,tol,iteMax,relax,'GS')
    if(meth == '3'):
        sol,residus = algo.algoSplit(K,F*0.,F,tol,iteMax,relax,'G')  
    if(meth == '4'):
        sol,residus = algo.GC(K,F*0.,F,tol,iteMax)
    if (meth=='5'):
        sol1,residus1= algo.algoSplit(K,F*0.,F,tol,iteMax,relax,'J')
        sol2,residus2= algo.algoSplit(K,F*0.,F,tol,iteMax,relax,'GS')
        sol3,residus3= algo.algoSplit(K,F*0.,F,tol,iteMax,relax,'G')  
        sol4,residus4= algo.GC(K,F*0.,F,tol,iteMax)
        amplificationDef = 1.0
        affichage.afficherDeformPont(listeNoeuds,listeArcs,listeFixation,sol1*amplificationDef)
        affichage.afficherDeformPont(listeNoeuds,listeArcs,listeFixation,sol2*amplificationDef) 
        affichage.afficherDeformPont(listeNoeuds,listeArcs,listeFixation,sol3*amplificationDef) 
        affichage.afficherDeformPont(listeNoeuds,listeArcs,listeFixation,sol4*amplificationDef) 
        affichage.afficher([np.arange(1, len(residus1) + 1), np.arange(1, len(residus2) + 1), np.arange(1, len(residus3) + 1), np.arange(1, len(residus4) + 1)],[np.log10(residus1), np.log10(residus2), np.log10(residus3), np.log10(residus4)],['tab:blue', 'tab:red', 'tab:green', 'tab:purple'])
        rang = np.linalg.matrix_rank(K)
    	
        print('Nb. ité.J :',K.shape,len(residus1))
        print("Err. rel. J: ",(residus1[-1]/residus1[0])) 
        print('Nb. ité.GS:',K.shape,len(residus2))
        print("Err. rel. GS: ",(residus2[-1]/residus2[0])) 
        print('Nb. ité. G:',K.shape,len(residus3))
        print("Err. rel. G: ",(residus3[-1]/residus3[0])) 
        print('Nb. ité.GC:',K.shape,len(residus4))
        print("Err. rel. GC: ",(residus4[-1]/residus4[0]))
        
        print("Le rang de la matrice K est",rang)
        is_symmetric = np.allclose(K, K.T)
        print("Symétrie de K :", is_symmetric)
        eigvals = np.linalg.eigvals(K)
        print("K est semi-définie positive :", np.all(eigvals >= -1e-12))
        cond = np.linalg.cond(K)
        print("Conditionnement de K :", cond)
    	
    if not (meth=='5'):
        print('Nb. ité. :',K.shape,len(residus))
# Affichage courbe norme résidus en fonction de l'itération.
        affichage.afficher([np.arange(1,len(residus)+1)],[np.log10(residus)],['tab:blue'])	
    
# Affichage de la déformation du point :
        amplificationDef = 1.0
        affichage.afficherDeformPont(listeNoeuds,listeArcs,listeFixation,sol*amplificationDef)
    
    
##----------------------------------    
    
    
##----------------------------------
# Test dessin :
if (butPrgm=="3"):
    
    [D,K] = data.generateMatricePont(listeNoeuds, listeArcs, listeMasse, listeFixation)
    
    # Nombre de modes calculé :
    nbMode = 1
    print("Combien de modes souhaitez-vous calculer ?")
    nbMode = input("Entrer le choix : ")
    nbMode = int(nbMode)
        
    # Resolution du probleme aux valeurs propres
    (dimX,dimY) = D.shape
    print(dimX,dimY)
    x0 = np.random.rand(dimY)*0.9+1.0
    [LL,VV] = algo.puInv(K,D,x0,nbMode,100,1e-5)
    
     
    # Affinage du calcul des elements propres  
    print("Affinage :")
    for i in range(nbMode):
        print("Mode ",i," / ",nbMode)
        #res = algo.puInv(K - (LL[i]-1e-10)*np.identity(dimX),D,VV[i],1,500,1e-10)
        res = algo.puInv(K - (LL[i]-1e-10)*D,D,VV[i],1,500,1e-10)
        LL[i] = res[0][0]+LL[i]
        VV[i] = np.copy(res[1][0])*1.0

    print("Liste valeurs propres : ")
    print(LL)

    choixMode = "0"
    while(choixMode != "-1"):
        print("Quel mode souhaitez-vous visualiser ? (Entrer un nombre entre 0 et Nmode-1, ou -1 pour arreter le programme)")
        choixMode = input("Entre le choix : ")
        mode = int(choixMode)
        
        if(choixMode == "-1"):
            break
        
        vv = np.random.rand(len(listeNoeuds)*3)*0.0
        listePts = []
        for i in range(len(listeNoeuds)):
            fixe = False
            for e in listeFixation:
                if (i==int(e)):
                    fixe = True
                    break
            if(fixe == False):
                listePts.append(i)
            
        for j,e in enumerate(listePts):
            vv[3*e] = VV[mode][3*j]
            vv[3*e+1] = VV[mode][3*j+1]
            vv[3*e+2] = VV[mode][3*j+2]        
    
        affichage.afficherAnimationPont(listeNoeuds, listeArcs, LL[mode],vv)
        
##----------------------------------
        

elif (butPrgm == "4"):
    [D, K] = data.generateMatricePont(listeNoeuds, listeArcs, listeMasse, listeFixation)
    print("Quelle type de source souhaitez-vous utiliser ?")
    print("1. Pression uniforme")
    print("2. Cisaillement")
    typeforce = input("Enter le choix : ")
    typeSource = 'Pression' if typeforce == '1' else 'Cisaillement'

    F = data.generateSourcePont(listeNoeuds, listeArcs, listeMasse, listeFixation, typeSource)

    tol = 1e-5
    iteMax = 1e3

    # Demande du nombre de divisions pour relax
    print("Combien de divisions pour la variation de relaxation (entre 0 et 2) ?")
    nb_div = int(input("Entrer le nombre de divisions : "))
    
    relax_values = np.linspace(0, 2, nb_div)

    residuals_per_method = {"Jacobi": [], "Gauss-Seidel": [], "Gradient": [], "GC": []}

    # Résolution pour chaque relaxation
    for relax in relax_values:
        sol1, residus1 = algo.algoSplit(K, F * 0., F, tol, iteMax, relax, 'J')
        sol2, residus2 = algo.algoSplit(K, F * 0., F, tol, iteMax, relax, 'GS')
        sol3, residus3 = algo.algoSplit(K, F * 0., F, tol, iteMax, relax, 'G')  
        sol4, residus4 = algo.GC(K, F * 0., F, tol, iteMax)

        # Ajouter tous les résidus pour chaque méthode
        residuals_per_method["Jacobi"].append(residus1[-1] if residus1 else np.nan)
        residuals_per_method["Gauss-Seidel"].append(residus2[-1] if residus2 else np.nan)
        residuals_per_method["Gradient"].append(residus3[-1] if residus3 else np.nan)
        residuals_per_method["GC"].append(residus4[-1] if residus4 else np.nan)

    # Affichage des courbes
    affichage.afficher(
        [relax_values, relax_values, relax_values, relax_values],
        [
            np.log10(residuals_per_method["Jacobi"]),
            np.log10(residuals_per_method["Gauss-Seidel"]),
            np.log10(residuals_per_method["Gradient"]),
            np.log10(residuals_per_method["GC"])
        ],
        ['tab:blue', 'tab:red', 'tab:green', 'tab:purple']
    )

    print("Analyse de l'effet de la relaxation complétée avec succès.")


##----------------------------------



print("FIN \n")
########################################
