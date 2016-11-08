import numpy as np

def var(X,Y,V,num0,scale):

    #-----------------------Matriz Triangular de Distancias---------------------------#

    D = np.array([[0]*len(Y)]*(len(X)),float) # D = matriz de distancias 
    
    for j in range (len(Y)):
        for i in range (len(X)):
            D[j,i] = scale*(np.sqrt((Y[j] - Y[i])**2 + (X[j] - X[i])**2))

    print np.max(D),' max D'


    DU = np.triu(D) # funcao matriz retangular/Matriz de distancias sem repeticao
    print np.max(DU),' max DU'
    
    #----------------------------------Variograma-------------------------------------#

    #-------------------------------------------------------------------------------
    num = num0+1

    rang = np.linspace(np.min(DU),np.max(DU),num) # escolhe o range de valores
    print 'range determinado pelos valores maximos e minimos:',len(rang),'elementos'
    print
    print rang

    range1 = [] # range da base
    range2 = [] # range da ponta

    #-------------------------------------------------------------------------------

    for i in range(num-1):
        range1.append(rang[i])
        range2.append(rang[i+1])

    #-------------------------------------------------------------------------------
    
    SH = np.shape(D)
    print SH,'forma da matriz de distancias'

    d1 = [] # Codigo
    d2 = [] # Range superior
    d3 = [] # Range Inferior
    d4 = [] # Valor da distancia
    d5 = [] # F(V[j],V[i])

    for k in range (num-1):
        for j in range (SH[0]):
            for i in range (SH[1]):
                if DU[j,i] > range1[k] and DU[j,i] <= range2[k]: # maior ao inves de maior e igual, retira os zeros
                    d1.append(k)
                    d2.append(range1[k])
                    d3.append(range2[k])
                    d4.append(DU[j,i])
                    d5.append((V[j] - V[i])**2)
            

    cod1 = np.linspace(min(d1),max(d1),max(d1)+1) # ou max(d1)+1, depende do valor de 'num'
    print
    print 'as ',num0,' categorias escolhidas em codigo '
    print
    print cod1
    print
    print '!!!ATENCAO!!!: estes numeros devem ser inteiros. Exemplo: (0.  1.  2. ...  n.)'
    print 'se nao forem inteiros, altere o valor de num0'

    vall = [0.0]*len(cod1)
    
    #------------------------------------------------------------------#
    SSD = [1]*(len(d1))

    ES = np.array([[0.0]*len(d1)]*len(cod1),float)
    EO = np.array([[0.0]*len(d1)]*len(cod1),float)


    for j in range (len(cod1)):
        for i in range (len(d1)):
            if d1[i] == j:
                ES[j,i] = d5[i]
                EO[j,i] = SSD[i]
            
    #------------------------------------------------------------------#

    SS  = ES.sum(axis=1) # soma dos valores das propriedades
    SSO = EO.sum(axis=1) # numero de ocorrencias da pripriedade no intervalo
    
    VI = []
    for i in range (len(SS)):
        VI.append((1/(2.0000*SSO[i]))*SS[i])
    
    rang_mid = [] # range de valor medio
    for i in range(len(range1)):
        rang_mid.append( (range1[i]+range2[i])/2.0 )  
    
    print
    print SS
    print SSO

    result = np.array([[0.0]*2]*len(VI),float)

    for i in range (len(VI)):
        result[i,0] = rang_mid[i]
        result[i,1] = VI[i]
    
    return(result)
    
    return(SS)
