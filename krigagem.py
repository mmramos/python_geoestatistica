import numpy as np

def kri(y,x,X,Y,V,grid,C0,C1,a,scale):
    
    AY2 = [y]+Y
    AX2 = [x]+X

    D = np.array([[0]*len(AY2)]*(len(AX2)),float)
    
    # D = Matriz de distancias 
    
    for j in xrange (len(AY2)):
        for i in xrange (len(AX2)):
            D[j,i] = scale*(np.sqrt((AY2[j] - AY2[i])**2 + (AX2[j] - AX2[i])**2))
            
            
    MC = np.array([[1]*(len(AX2))]*(len(AY2)),float)
    MC[len(AY2)-1,len(AX2)-1] = 0
    
    # MC = matriz de covariancias
    
    for j in xrange (len(AY2)-1):
        for i in xrange (len(AX2)-1):
            if D[j,i] == 0:
                MC[j,i] = C0 + C1
            else:
                MC[j,i] = C1*np.exp(-3*(D[j+1,i+1])/a)
                
    # NC = matriz de covariancias com o ponto
    
    NC = [1]*len(AX2)
    
    for i in xrange (1,len(AX2)):
        NC[i-1] = C1*np.exp(-3*(D[0,i])/a)
        
    # CI = inversa da matriz de covariancias
    
    CI = np.linalg.inv(MC)
    
    # W =  matriz de ponderacoes + parametro de lagrange
    
    W = np.dot(CI,NC)
    
    V0 = []
    for i in xrange (len(X)):
        V0.append(V[i]*W[i])
    
    VX0 = sum(V0)
    
    return(VX0)

#------------------------Esferica-----------------------#

def kriges(y,x,X,Y,V,grid,C0,C1,a,scale):
    
    AY2 = [y]+Y
    AX2 = [x]+X
    
    D = np.array([[0]*len(AY2)]*(len(AX2)),float)
    
    # D = Matriz de distancias 
    
    for j in range (len(AY2)):
        for i in range (len(AX2)):
            D[j,i] = scale*(np.sqrt((AY2[j] - AY2[i])**2 + (AX2[j] - AX2[i])**2))
            
            
    MC = np.array([[-1]*(len(AX2))]*(len(AY2)),float)
    MC[len(AY2)-1,len(AX2)-1] = 0
    
    # MC = matriz de covariancias
    
    for j in range (len(AY2)-1):
        for i in range (len(AX2)-1):
            if D[j,i] > a:
                MC[j,i] = C0 + C1
            else:
                MC[j,i] = C0 + C1*(1.5*(D[j+1,i+1] / a) + 0.5*(D[j+1,i+1] / a))
                
    # NC = matriz de covariancias com o ponto
    
    NC = [1]*len(AX2)
    
    for i in range (1,len(AX2)):
        if D[j,i] > a:
            NC[i-1] = C0 + C1
        else:
            NC[i-1] = C0 + C1*(1.5*(D[0,i] / a) + 0.5*(D[0,i] / a))
        
    # CI = inversa da matriz de covariancias
    
    CI = np.linalg.inv(MC)
    
    # W =  matriz de ponderacoes + parametro de lagrange
    
    W = np.dot(CI,NC)
    
    V0 = []
    for i in range (len(X)):
        V0.append(V[i]*W[i])
    
    VX0 = sum(V0)
    
    return(VX0)

#------------------------Exponencial-----------------------#

def krigex(y,x,X,Y,V,grid,C0,C1,a,scale):
    
    AY2 = [y]+Y
    AX2 = [x]+X
    
    D = np.array([[0]*len(AY2)]*(len(AX2)),float)
    
    # D = Matriz de distancias 
    
    for j in range (len(AY2)):
        for i in range (len(AX2)):
            D[j,i] = scale*(np.sqrt((AY2[j] - AY2[i])**2 + (AX2[j] - AX2[i])**2))
            
            
    MC = np.array([[-1]*(len(AX2))]*(len(AY2)),float)
    MC[len(AY2)-1,len(AX2)-1] = 0
    
    # MC = matriz de covariancias
    
    for j in range (len(AY2)-1):
        for i in range (len(AX2)-1):
            MC[j,i] = C0 + C1*(1 - np.exp( -D[j+1,i+1] / a))
                
    # NC = matriz de covariancias com o ponto
    
    NC = [1]*len(AX2)
    
    for i in range (1,len(AX2)):
        NC[i-1] = C0 + C1*(1 - np.exp( -D[0,i] / a))
        
    # CI = inversa da matriz de covariancias
    
    CI = np.linalg.inv(MC)
    
    # W =  matriz de ponderacoes + parametro de lagrange
    
    W = np.dot(CI,NC)
    
    V0 = []
    for i in range (len(X)):
        V0.append(V[i]*W[i])
    
    VX0 = sum(V0)
    
    return(VX0)

#-----------------------Gaussiana-------------------------#

def krigaus(y,x,X,Y,V,grid,C0,C1,a,scale):
    
    AY2 = [y]+Y
    AX2 = [x]+X
    
    D = np.array([[0]*len(AY2)]*(len(AX2)),float)
    
    # D = Matriz de distancias 
    
    for j in range (len(AY2)):
        for i in range (len(AX2)):
            D[j,i] = scale*(np.sqrt((AY2[j] - AY2[i])**2 + (AX2[j] - AX2[i])**2))
            
            
    MC = np.array([[-1]*(len(AX2))]*(len(AY2)),float)
    MC[len(AY2)-1,len(AX2)-1] = 0
    
    # MC = matriz de covariancias
    
    for j in range (len(AY2)-1):
        for i in range (len(AX2)-1):
            if D[j,i] > a:
                MC[j,i] = C0 + C1
            else:
                MC[j,i] = C0 + C1*(1 - np.exp( -(D[j+1,i+1]**2) / (a**2)))
                
    # NC = matriz de covariancias com o ponto
    
    NC = [1]*len(AX2)
    
    for i in range (1,len(AX2)):
        if D[j,i] > a:
            NC[i-1] = C0 + C1
        else:
            NC[i-1] = C0 + C1*(1 - np.exp( -(D[0,i]**2) / (a**2)))
        
    # CI = inversa da matriz de covariancias
    
    CI = np.linalg.inv(MC)
    
    # W =  matriz de ponderacoes + parametro de lagrange
    
    W = np.dot(CI,NC)
    
    V0 = []
    for i in range (len(X)):
        V0.append(V[i]*W[i])
    
    VX0 = sum(V0)
    
    return(VX0)
