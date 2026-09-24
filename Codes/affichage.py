import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpImg
from mpl_toolkits.mplot3d import Axes3D

from pylab import *
       
        
def afficherPont(listNoeuds,listArcs):
    #plt.ion()
    #plt.hold('on')
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    
    minX = 0.0
    maxX = 0.0
    minY = 0.0
    maxY = 0.0
    minZ = 0.0
    maxZ = 0.0
    for e in listNoeuds:
        x = e[0]
        y = e[1]
        z = e[2]
        if (minX >= x):
            minX = x
        if (maxX <= x):
            maxX = x
        if (minY >= y):
            minY = y
        if (maxY <= y):
            maxY = y
        if (minZ >= z):
            minZ = z
        if (maxZ <= z):
            maxZ = z

    for j,aa in enumerate(listArcs):
        i1 = int(aa[0])
        i2 = int(aa[1])
        p1 = listNoeuds[i1]
        p2 = listNoeuds[i2]
        ax.plot([p1[0],p2[0]],[p1[1],p2[1]],[p1[2],p2[2]],linewidth=2.0,color="black")
        
    ax.set_xlim3d([minX-0.2, maxX+0.2])
    ax.set_ylim3d([minY-0.2, maxY+0.2])
    ax.set_zlim3d((minZ-0.2, maxZ+0.2))
        
    plt.show()
    
    
def afficherDeformPont(listNoeuds,listArcs,listeFixation,deform):
    listePts = []
    for i in range(len(listNoeuds)):
        fixe = False
        for e in listeFixation:
            if (i==int(e)):
                fixe = True
                break
        if(fixe == False):
            listePts.append(i)
            
            
    minX = np.min(listNoeuds[:,0])
    maxX = np.max(listNoeuds[:,0])
    minY = np.min(listNoeuds[:,1])
    maxY = np.max(listNoeuds[:,1])
    minZ = np.min(listNoeuds[:,2])
    maxZ = np.max(listNoeuds[:,2])
    
    
    vv = np.zeros(len(listNoeuds)*3)    
    for j,e in enumerate(listePts):
        vv[3*e] = deform[3*j]
        vv[3*e+1] = deform[3*j+1]
        vv[3*e+2] = deform[3*j+2]  
        
    fig = plt.figure()    
    plt.clf()
    #ax = fig.gca(projection='3d')
    ax = fig.add_subplot(projection = '3d')
    
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    for j,aa in enumerate(listArcs):
        i1 = int(aa[0])
        i2 = int(aa[1])
        p1 = listNoeuds[i1]
        p2 = listNoeuds[i2]
        
        
        dp1 = np.array([vv[i1*3], vv[i1*3+1], vv[i1*3+2]])
        dp2 = np.array([vv[i2*3], vv[i2*3+1], vv[i2*3+2]])
        
        n1 = np.sqrt(np.dot(dp1,dp1)) / (0.1)
        n2 = np.sqrt(np.dot(dp2,dp2)) / (0.1)
        n1 = min(n1,1.0)
        n2 = min(n2,1.0)
        
        dx = (p2[0]+vv[i2*3] - (p1[0]+vv[i1*3]))
        dy = (p2[1]+vv[i2*3+1] - (p1[1]+vv[i1*3+1]))
        dz = (p2[2]+vv[i2*3+2] - (p1[2]+vv[i1*3+2]))
        
        ax.plot([p1[0]+vv[i1*3],p1[0]+vv[i1*3]+dx/3],[p1[1]+vv[i1*3+1],p1[1]+vv[i1*3+1]+dy/3],[p1[2]+vv[i1*3+2],p1[2]+vv[i1*3+2]+dz/3],linewidth=2.0,color=(n1,0,0))
        ax.plot([p1[0]+vv[i1*3]+dx/3,p1[0]+vv[i1*3]+2*dx/3],[p1[1]+vv[i1*3+1]+dy/3,p1[1]+vv[i1*3+1]+2*dy/3],[p1[2]+vv[i1*3+2]+dz/3,p1[2]+vv[i1*3+2]+2*dz/3],linewidth=2.0,color=((n1+n2)*0.5,0,0))
        ax.plot([p1[0]+vv[i1*3]+2*dx/3,p1[0]+vv[i1*3]+dx],[p1[1]+vv[i1*3+1]+2*dy/3,p1[1]+vv[i1*3+1]+dy],[p1[2]+vv[i1*3+2]+2*dz/3,p1[2]+vv[i1*3+2]+dz],linewidth=2.0,color=(n2,0,0))
            
    ax.quiver(listNoeuds[:,0],listNoeuds[:,1],listNoeuds[:,2],vv[::3],vv[1::3],vv[2::3],color='k')
    ax.set_xlim3d([minX-0.2, maxX+0.2])
    ax.set_ylim3d([minY-0.2, maxY+0.2])
    ax.set_zlim3d((minZ-0.2, maxZ+0.2))    
    
    plt.show()    
          
    
def afficherAnimationPont(listNoeuds,listArcs,ll,vv):
    plt.ion()
    #plt.hold('off')
    fig = plt.figure()
    #ax = fig.gca(projection='3d')
    mm = max(vv)
    
    minX = 0.0
    maxX = 0.0
    minY = 0.0
    maxY = 0.0
    minZ = 0.0
    maxZ = 0.0
    for e in listNoeuds:
        x = e[0]
        y = e[1]
        z = e[2]
        if (minX >= x):
            minX = x
        if (maxX <= x):
            maxX = x
        if (minY >= y):
            minY = y
        if (maxY <= y):
            maxY = y
        if (minZ >= z):
            minZ = z
        if (maxZ <= z):
            maxZ = z
    
    for i in range(87):
        tt = np.sin((2*i*np.pi)/50)
        vt = vv*tt*0.2/mm
        plt.clf()
        #ax = fig.gca(projection='3d')
        ax = fig.add_subplot(projection = '3d')
        
        ax.grid(False)
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        for j,aa in enumerate(listArcs):
            i1 = int(aa[0])
            i2 = int(aa[1])
            p1 = listNoeuds[i1]
            p2 = listNoeuds[i2]
            
            
            dp1 = np.array([vt[i1*3], vt[i1*3+1], vt[i1*3+2]])
            dp2 = np.array([vt[i2*3], vt[i2*3+1], vt[i2*3+2]])
            
            n1 = np.sqrt(np.dot(dp1,dp1)) / (0.2)
            n2 = np.sqrt(np.dot(dp2,dp2)) / (0.2)
            n1 = min(n1,1.0)
            n2 = min(n2,1.0)
            
            dx = (p2[0]+vt[i2*3] - (p1[0]+vt[i1*3]))
            dy = (p2[1]+vt[i2*3+1] - (p1[1]+vt[i1*3+1]))
            dz = (p2[2]+vt[i2*3+2] - (p1[2]+vt[i1*3+2]))
            ax.plot([p1[0]+vt[i1*3],p1[0]+vt[i1*3]+dx/3],[p1[1]+vt[i1*3+1],p1[1]+vt[i1*3+1]+dy/3],[p1[2]+vt[i1*3+2],p1[2]+vt[i1*3+2]+dz/3],linewidth=2.0,color=(n1,0,0))
            ax.plot([p1[0]+vt[i1*3]+dx/3,p1[0]+vt[i1*3]+2*dx/3],[p1[1]+vt[i1*3+1]+dy/3,p1[1]+vt[i1*3+1]+2*dy/3],[p1[2]+vt[i1*3+2]+dz/3,p1[2]+vt[i1*3+2]+2*dz/3],linewidth=2.0,color=((n1+n2)*0.5,0,0))
            ax.plot([p1[0]+vt[i1*3]+2*dx/3,p1[0]+vt[i1*3]+dx],[p1[1]+vt[i1*3+1]+2*dy/3,p1[1]+vt[i1*3+1]+dy],[p1[2]+vt[i1*3+2]+2*dz/3,p1[2]+vt[i1*3+2]+dz],linewidth=2.0,color=(n2,0,0))
            
 
        ax.set_xlim3d([minX-0.2, maxX+0.2])
        ax.set_ylim3d([minY-0.2, maxY+0.2])
        ax.set_zlim3d((minZ-0.2, maxZ+0.2)) 
        
        plt.show()
        plt.pause(1e-5)  
            
    #plt.hold('on')
    



def afficher(list_x,list_y,attribut):
    for j,absic in enumerate(list_x):
        plt.plot(absic,list_y[j],attribut[j],linewidth=2.0,)
    plt.show()
