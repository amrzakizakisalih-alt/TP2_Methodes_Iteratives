import numpy as np
import copy

import data


#Algorithme du gradient conjugué
#Entree : A (matrice), x0 (vecteur), b (vecteur), iteMax (entier), relax (double), choix (caractere)
#Sortie : x (vecteur), residus (liste double)
def algoSplit(A,x0,b,eps,iteMax,relax,choix):
    (n,p) = A.shape
    x = np.copy(x0)
    r = np.dot(A,x) - b
    w = r*1.0
    residus = [np.sqrt(np.dot(r,r))]
    ite = 1
    while (residus[-1] >= eps*residus[0] and iteMax > ite):
        
        print("Ite : ",ite," / ",iteMax )
        print("Err. rel. : ",residus[-1]/residus[0])
        
        y = np.zeros(n)
        # Resolution de y = M^{-1} r
        if (choix == 'G'):
            # Algorithme du gradient, M = Id / alpha :
            # Completer ICI :
            alpha = np.dot(r, r) / np.dot(np.dot(r, A), r)  # Calcul de alpha
            y = alpha * r  # Mise à jour du vecteur y
            
        if (choix == 'J'):
            # Algorithme de Jacobi, M = D :
            # Completer ICI :
          
            D = np.diag(np.diag(A))
            Dinv=np.linalg.inv(D)
            y=np.dot(Dinv,r)
            
            
            
              
        if (choix == 'GS'):
            # Algorithme de Gauss-Seidel, M = D - E :
            # Completer ICI :
             
            D = np.diag(np.diag(A))
            E = -np.tril(A,-1)
            M=D-E
            Minv = np.linalg.inv(M) # Inversion de D+L
            y=np.dot(Minv,r)

        
        # Relaxation : 
        #Completer ICI :    
        x = (1.0-relax)*x + relax*(x-y)
        
        # Calcul du résidu :
        r = np.dot(A,x) - b
        residus.append(np.sqrt(np.dot(r,r)))
        ite += 1
    print("Ite : ",ite," / ",iteMax )
    print("Err. rel. : ",(residus[-1]/residus[0]))    
    return [x,residus]
    
    
#Algorithme du gradient conjugué
#Entree : A (matrice), x0 (vecteur), b (vecteur), eps (double), iteMax (entier)
#Sortie : x (vecteur), residus (liste double)    
def GC(A,x0,b,eps,iteMax):
    (n,p) = A.shape
    x = np.copy(x0)
    r = np.dot(A,x) - b
    p = -r*1.0
    residus = [np.sqrt(np.dot(r,r))]
    ite = 1
    
    while (residus[-1] >= eps*residus[0] and iteMax > ite):
        
        print("Ite : ",ite," / ",iteMax )
        print("Err. rel. : ",(residus[-1]/residus[0]))
        
        alpha = np.dot(r, r) / np.dot(p, np.dot(A, p))
        x = x + alpha * p
        r_n = r + alpha * np.dot(A, p)
        beta = np.dot(r_n, r_n) / np.dot(r, r)
        p = -r_n + beta * p
        r = r_n
        
        residus.append(np.sqrt(np.dot(r, r)))
        ite += 1
        
        
    print("Ite : ",ite," / ",iteMax )
    print("Err. rel. : ",(residus[-1]/residus[0]))    
    return [x,residus]    
    
    
# Algorithme de la puissance itérée:
# Entree: K (matrice), M (matrice), x0 (vecteur), m (entier), iteMax (entier), err (double)
# Sortie : vp (liste v.p.), V (liste vecteurs propre)    

def puInv(K,M,x0,m,iteMax,err):

    # Factorisation QR de K, faite une fois.
    [Q,R] = qr(K)
    (n,p) = K.shape

    # On determine les modes statiques, c est a dire le noyau de K :
    DDzero = []
    for i in range(n):
        if(np.abs(R[i,i]) < 1e-15):
            DDzero.append(i)

    VVzero = []
    for e in DDzero:
        y = remontee(R[0:e,0:e],R[0:e,e])
        vv = x0*0
        vv[0:e] = -y
        vv[e] = 1
        normvv = np.sqrt(np.dot(np.dot(M,vv),vv))
        vv = vv / normvv
        VVzero.append(vv)

    vp = []
    V = []
    ite = 1
    x00 = np.copy(x0)*1.0
    x1=x0
    for k in range(m):
        print("Etape ",k)
        x0 = x00  + np.random.rand(len(x0))*err

        # On elimine les modes statiques
        for i in VVzero:
            x0 = x0 - np.dot(np.dot(M,x0),i)*i

        # Deflation :
        for i in range(k):
            x0 = x0 - np.dot(np.dot(M,x0),V[i])*V[i]

        # M-Normalisation :
        x0 = x1 / np.sqrt(np.dot(np.dot(M,x0),x0))

        # Calcul de x1 = K^{-1} M x0 :
        ty = np.dot(Q.transpose(),np.dot(M,x0))
        y = remontee(R[0:p,0:p],ty[0:p])

        # M-Normalisation :
        x1 = y / np.sqrt(np.dot(np.dot(M,y),y))

        # On commence les iterations :
        ite = 0
        while(np.abs(np.abs(np.dot(x0,x1))  -1) > err and ite < iteMax):
            # On elimine les modes statiques
            for i in VVzero:
                x1 = x1 - np.dot(np.dot(M,x1),i)*i
            # Deflation :
            for i in range(k):
                # Completer ICI :
                x1 = x1-np.dot(np.dot(M,x1),V[i])*V[i]

            # M-Normalisation :
            # Completer ICI :
            x0 =x1/np.sqrt(np.dot(np.dot(M,x1),x1))
            

            # Calcul de x_{k+1} = K^{-1} M x_k :
            # Completer ICI :
            ty = np.dot(Q.transpose(),np.dot(M,x0))
            y = remontee(R[0:p,0:p],ty[0:p])

            # M-Normalisation :
            x1 = y / np.sqrt(np.dot(np.dot(M,y),y))

            # Iteration suivante :
            ite = ite+1

        print("Nb ite :",ite,np.abs(np.abs(np.dot(x0,x1)) -1))

        # On ajoute le vecteur x0 à la liste de vecteurs propres
        V.append(x0)
        jmax = 0
        emax = -1.0
        for j,e in enumerate(y):
            if np.abs(e) >= emax:
                jmax = j
                emax = np.abs(e)
        # et on ajoute la nouvelle valeur propre.
        vp.append(x0[jmax]/y[jmax])
        print("valeur propre ",vp[-1])
        print("\n")
    return [vp,V]        
    
    

#################################################
# Factorisation QR :
#Entree : A (matrice)
#Sortie : Q (matrice), R (matrice)
def qr(A):
    (n,p) = A.shape
    # Declaration des variables et initialisation
    Q = np.zeros((n,p))
    R = np.zeros((p,p))
    
    
    # Orthonormalisation G.S.
    for i in range(p):
        v = A[:,i]
        for j in range(i):
            scal = np.dot(v,Q[:,j])
            v = v - scal * Q[:,j]
            R[j,i] = scal
        normV = np.sqrt(np.dot(v,v))
        R[i,i] = normV
        if (normV > 1e-15):
            Q[:,i] = v / normV
    return [Q,R]    
    
# Resolution de U x = b
#  - Entree : U (matrice), b (vecteur)
#  - Sortie : x (Vecteur)
def remontee(U,b):
    (n,p) = U.shape
    x = np.zeros(n)
    for l in range(n-1,-1,-1):
        x[l] = b[l]
        for c in range(l+1,n):
            x[l] = x[l] - U[l,c]*x[c]
        if (np.abs(U[l,l]) < 1e-10 ):
            x[l] = 0.0
        else :
            x[l] = x[l] / U[l,l]

    return x 
#################################################       
    
