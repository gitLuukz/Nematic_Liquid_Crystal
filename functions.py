import numpy as np

import plotfunctions as pf



def grid(width,length,N,theta):
    """
    Returns the starting configuration of the system

    Parameters
    ----------
    width : TYPE
        Width of the grid
    length : TYPE
        Length of the grid. Needs to be even to have a correct grid with periodic boundary conditions
    N : TYPE
        Iterations that will be computed in the system 

    Returns
    -------
    R : TYPE
        Matrix with x,y coordinates for first step and is empty for other steps
    dx : TYPE
        first direction of x component of the molecus
    dy : TYPE
        first direction of y component of the molecus
    theta : TYPE
        matrix containing rotation of all molecules for N time steps size: (#molecules x N)

    """
    x1 = np.linspace(0,width-1,width)+3/4
    x2 = np.linspace(0,width-1,width)+1/4
    
    # length should be an even number to make sure that the grid will be correct when looked at with periodic boundaries
    
    
    R = np.zeros((2,width*length,N)) 
    for i in range(length):
        indices = (np.linspace(0,width-1,width)+ i*width).astype(int)
        y = np.zeros(width) + i + 2/4
        if i%2 == 0:
            R[0,indices,0] = x1
            R[1,indices,0] = y
        elif i%2 == 1: 
            R[0,indices,0] = x2
            R[1,indices,0] = y
    
    #thetagrid = np.zeros((width*length,N))
    dx,dy = angle(theta[:,0])
    return R, dx,dy

def angle(theta):
    """
    

    Parameters
    ----------
    Returns parameters that can be used to plot the configuration.
    
    theta : TYPE
        An array containing the angles at step n of N

    Returns
    -------
    dx : TYPE
        the x component of arrows
    dy : TYPE
        the y component of arrows

    """
    dx = np.sin(theta)
    dy = np.cos(theta)
    return dx,dy

def update(Uiold,Uinew,Ui,UT,theta,length,width,epsilon,T,mu,N,M,delta,k,acceptance_r,Cv,NT,kT,T2,T1):
    UT[0] = -6*M
    avUT = np.zeros((N,NT))
    scaledenergy = np.zeros((N,NT))
    scaledT = np.zeros((NT))
    for tr in range(NT):
        theta = np.zeros((length,width,N))
        Uiold = np.zeros((length,width,N+1))
        Uinew = np.zeros((length,width,N))
        Ui = np.zeros((length,width,N))
        
        for n in range(N):    
            for i in range(length):
                for j in range(width):
                    thetaold = theta[i,j,n-1]
                    Uinew[i,j,n],theta[i,j,n] = interactionenergy(epsilon,mu,thetaold,delta,i,j,n,length,width,theta)
                    dUi = Uinew[i,j,n] - Uiold[i,j,n]
                    if dUi < float(0):
                        Ui[i,j,n] = Uiold[i,j,n]+dUi
                        Uiold[i,j,n+1] = Uinew[i,j,n]
                        acceptance_r[n,tr] += 1
                    
                    else: 
                        if np.exp(-(dUi)/(kT[tr])) > np.random.uniform(0,1):
                            Ui[i,j,n] = Uiold[i,j,n]+dUi
                            Uiold[i,j,n+1] = Uinew[i,j,n]
                            acceptance_r[n,tr] += 1
                        
                        else: 
                            Ui[i,j,n] = Uiold[i,j,n]
                            Uiold[i,j,n+1] = Uiold[i,j,n]
                            theta[i,j,n] = theta[i,j,n-1]
            
            
            T2[n,tr] = apolar(theta[:,:,n-1],M)
            T1[n,tr] = polar(theta[:,:,n-1],M) 
            #Calc properties  
            if float(n)>=(N-N/2):
                
                Cv[n,tr] = (np.mean(Ui[:,:,n]**2) - (np.mean(Ui[:,:,n]))**2)/(kT[tr]**2) # square of first M shouldn ot be there # there was also a k**2 that did not belong there
                #cv(Ui[:,:,n],M)#
                # R,dx,dy = grid(width, length, N, theta)
                
                #fignr = n+1
                #pf.orientation(fignr, R[:,:,0], dx, dy, width, length)    
            
                avUT[n,tr] = UT[n+1,tr]/M
                scaledenergy[n,tr] = avUT[n,tr]#/epsilon
                scaledT[tr] = kT[tr]#/epsilon
            UT[n+1,tr] = 1/2*internalenergy(Ui[:,:,n])
    acceptance_r = acceptance_r/M
                
    return theta, UT,avUT,scaledenergy, acceptance_r,Cv,T2,T1,Ui
    
    
def internalenergy(Ui):
    UT = np.sum(Ui)
    return UT

def cv(Ui,M):
    icv = 2/(3*M) - (np.mean(Ui**2) - (np.mean(Ui))**2)/(np.mean(Ui))**2
    cv = 1/icv
    return cv
    
def interactionenergy(epsilon,mu,thetaold,delta,i,j,n,length,width,theta):
    thetaparticle = thetachange(thetaold,delta)
    #print('thetaparticle = ',thetaparticle)
    [thetaijneighbour,l,w] = thetainteract(i,j,n,length,width,theta) 
    # print(thetaijneighbour)
    thetaij = thetaparticle - thetaijneighbour
    #print('thetaij = ', thetaij)
    Ui = pairpotential(thetaij,mu,epsilon)
    #Ui = np.sum(Uij)
    return Ui,thetaparticle
    
def pairpotential(thetaij,mu,epsilon):
    #mu = 0 #yes if apolar
    Ui = np.sum(-epsilon*(np.cos(2*thetaij)+mu*np.cos(thetaij)))
    return Ui    
    
def thetachange(thetaold,delta):
    thetanew = thetaold + (np.random.uniform(0,2)-1)*delta#*180/np.pi    
    return thetanew
    
def thetainteract(i,j,n,length,width,theta):
    w = 2
    l = 3
    thetaij = np.zeros((l,w))
#    if i-1< float(0):
#        if j+1>=width:
#            thetaij[l-3,w-1] = theta[i-1+length,j+1-width,n] 
#        else:
#            thetaij[l-3,w-2] = theta[i-1+length,j,n]
#    if i+1>=length:
#        if j+1>=width:
#            thetaij[l-1,w-1] = theta[i+1-length,j+1-width,n]
#        else:
#            thetaij[l-1,w-2] = theta[i+1-length,j,n]
#    else:
#        if j-1<float(0):
#            thetaij[l-1,w-2] = theta[i,j-1+width,n]
#        if j+1>=width:
#            thetaij[l-2,w-1] = theta[i,j+1-width,n]
#        else:
#            thetaij[l-3,w-2] = theta[i-1,j,n]
#            thetaij[l-3,w-1] = theta[i-1,j+1,n]
#            thetaij[l-2,w-2] = theta[i,j-1,n]
#            thetaij[l-2,w-1] = theta[i,j+1,n]
#            thetaij[l-1,w-2] = theta[i+1,j,n]
#            thetaij[l-1,w-1] = theta[i+1,j+1,n]
#            
    
    thetatest = np.zeros((3,3))
#    print(i,j,(j+1)*(j<length-1)+0)
    a = (i+1)*(i<length-1)+0
    b = (j+1)*(j<width-1)+0
    thetatest[0,:] = np.array([theta[i-1,j-1,(n)-(j<1 or i<1)], theta[i-1,j,(n)-(i<1)], theta[i-1,b,(n-1)+(j+1>width)]]) 
    thetatest[1,:] = np.array([theta[i,j-1,(n)-(j<1)], theta[i,j,n], theta[i,b,(n-1)+(j+1>width)]])            
    thetatest[2,:] = np.array([theta[a,j-1,(n-1)+(i+1>length and j>0)], theta[a,j,(n-1)+(i+1>length)], theta[a,b,(n-1)+(i+1>length or j+1>width)]])         
    
    if (i % 2) == 0:
        thetaij[0,0] = thetatest[0,0]
        thetaij[0,1] = thetatest[0,1]
        thetaij[1,0] = thetatest[1,0]
        thetaij[1,1] = thetatest[1,2]
        thetaij[2,0] = thetatest[2,0]
        thetaij[2,1] = thetatest[2,1]

    if (i % 2) == 1:
        thetaij[0,0] = thetatest[0,1]
        thetaij[0,1] = thetatest[0,2]
        thetaij[1,0] = thetatest[1,0]
        thetaij[1,1] = thetatest[1,2]
        thetaij[2,0] = thetatest[2,1]
        thetaij[2,1] = thetatest[2,2]
    #print('neighbour = ', thetaij)
    return thetaij,l,w
    

def apolar(theta,M):
    tanalfa = np.sum(np.sin(2*theta[:,:]))/np.sum(np.cos(2*theta[:,:]))
    alfa = 1/2*np.arctan(tanalfa)    
    T2 = 1/M*np.sum(np.cos(2*(theta[:,:]-alfa)))
    return T2    
    
def polar(theta,M): 
    tanalfap = np.sum(np.sin(theta[:,:]))/np.sum(np.cos(theta[:,:]))
    alfap = np.arctan(tanalfap)    
    T1 = 1/M*np.sum(np.cos(theta[:,:]-alfap))
    return T1
    
    
def error(a):
    N = len(a)
    x = np.zeros(N-1)
    for t in range(N-1):
        n = np.linspace(0,N-t-1,N-t,dtype = int)
        # print(n)
        top  = ( (N-t)* np.sum(a[n]*a[n+t]) - np.sum(a[n])*np.sum(a[n+t]) ) 
        bot = ( np.sqrt( (N-t)*np.sum(a[n]**2) - (np.sum(a[n]))**2) * np.sqrt((N-t)*np.sum(a[n+t]**2) - np.sum(a[n+t])**2) )
        # print(top,bot)
        # print(n,n+t)
        x[t] = top/bot
    return x

def function(t,tau):
    return np.exp(-t/tau)