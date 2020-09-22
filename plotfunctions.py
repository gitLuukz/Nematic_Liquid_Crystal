import matplotlib.pyplot as plt

def orientation(fignr,R,dx,dy,width,length):
    """

    Parameters
    ----------
    fignr : 
        Number of the figure.
    R : 
        2xn matrix containing locations of moleculse
    dx : 
        x component of the direction of the arrows
    dy : 
        y component of the direction of the arrows
    width : 
        Grid width
    length : 
        Grid length

    Returns
    -------
    Plot of the configuration at a time step

    """
    plt.figure(fignr,figsize = (5,5))
    plt.quiver( R[0,:],R[1,:]-1/2, dx ,dy)
    plt.grid(which = 'both')
    plt.axis([-1,width+1,-1,length+1])
    
    
def specificheat(Trange,Cv,fignr,sigma):
    plt.figure(fignr,figsize = (5,5))
    plt.errorbar(Trange,Cv,sigma, marker = '.', color = 'blue', ls = '')
    plt.xlim([0,4])
    plt.title('Specific heat')
    plt.xlabel('T*')
    plt.ylabel('$C_v*$')
    plt.savefig('Cv.png')
    
def apolar(Trange,T2,fignr,sigma):
    plt.figure(fignr,figsize = (5,5))
    plt.errorbar(Trange,T2,sigma, marker = '.', color = 'blue', ls = '')
    plt.xlim([0,4])
    plt.ylim([-0.25,1.2])
    plt.title('Apolar')
    plt.xlabel('T*')
    plt.ylabel('$T2*$')
    plt.savefig('T2.png')
    
    
def polar(Trange,T1,fignr,sigma):
    plt.figure(fignr,figsize = (5,5))
    plt.errorbar(Trange,T1,sigma, marker = '.', color = 'blue', ls = '')
    plt.xlim([0,8])
    plt.ylim([-0.25,1.2])
    plt.title('Polar')
    plt.xlabel('T*')
    plt.ylabel('$T1*$')
    plt.savefig('T1.png')
    
def both(Trange,T1,sigma1,T2,sigma2,fignr):
    plt.figure(fignr,figsize= (5,5))
    plt.errorbar(Trange,T1,sigma1, marker = '.', color = 'blue', ls = '',label = 'polar')
    plt.errorbar(Trange,T2,sigma2, marker = '+', color = 'black', ls = '',label = 'Apolar')
    plt.xlim([0,8])
    plt.ylim([-0.25,1.2])
    plt.legend()
    plt.xlabel('T*')
    plt.title('Long-range order parameter')
    plt.savefig('polar_apolar.png')
    
def Energy(kT, U, sigmaU, fignr):
    plt.figure(fignr,figsize = (5,5))
    plt.errorbar(kT, U, sigmaU, marker = '.', color = 'blue', ls = '')
    plt.xlabel('T*')
    plt.ylabel('U*')
    plt.xlim([0,8])
    plt.ylim([-1200,0])
    plt.title('Energy')
    plt.savefig('Scaled_energy.png')
    
def energy_time(U,fignr):
    plt.figure(fignr,figsize = (5,5))
    plt.plot(U)
    plt.xlabel('Iteration')
    plt.ylabel('Energy')
    plt.title('Energy over iterations')
    plt.savefig('energy_iterations.png')