from celluloid import Camera
import numpy as np
import matplotlib.pyplot as pp

N = 2000 # Number of people
q = 0.0 # Percentage of vaccinated people
D = 95 # Size of Domain

T_recover = 90         # time to recover from illness
T_immune = 90         # time after illness in which one is immune
Spreading_distance = 1.4 # distance at which infected individuals infect others
T = 600                # total running time
    
def dist(A,B):
    # A = array of points
    # B = single point
    return np.sqrt((A[:,:1]-B[0])**2 + (A[:,1:2]-B[1])**2)


lambda_f = 1

def random_walk(A):
    delta_x = np.random.uniform(-lambda_f,lambda_f,A.shape[0])[:,np.newaxis]
    delta_y = np.random.uniform(-lambda_f,lambda_f,A.shape[0])[:,np.newaxis]
    A[:,:1] = A[:,:1]+delta_x
    A[:,1:2] = A[:,1:2]+delta_y
    # Reflection at the boundary:
    A[A[:,0]>D,:1] = -2*A[A[:,0]>D,:1]+3*D
    A[A[:,1]>D,1:2] = 3*D-2*A[A[:,1]>D,1:2]
    A[A[:,0]<0,:1] = -A[A[:,0]<0,:1]
    A[A[:,1]<0,1:2] = -A[A[:,1]<0,1:2]

    return A

# Susceptibles:
S = np.random.rand(N,2)*D

I = np.zeros(N, dtype='int')
I[0]=1
# 0 = Susceptible
# 1 = Infected
# 2 = Immune

Times = np.zeros(N)
Times[0]=T_recover

Times_immune = np.zeros(N)


f, (P1,P2)=pp.subplots(ncols=2, figsize=(12,5))
camera = Camera(f)

L = []
m=0
n=1
ctr=0
while len(I[I==1])!=0 and ctr<1000:
# for r in range(T):
    print(n)
    L.append(len(I[I==1]))
    fill = [np.nan for k in range(T-len(L)+1)]
    P1.plot(range(len(L+fill)),L+fill, color='r')
    P1.axis([0, T, 0, N])
    P2.plot(S[:,0],S[:,1],'.', c='b')
    P2.plot(S[I==1,0],S[I==1,1],'.', c='r')
    P2.plot(S[I==2,0],S[I==2,1],'.', c='y')
    camera.snap()
    # pp.plot(Vx,Vy,'.', c='y')
    while n!=m:
        n = len(I[I==1])
        infected = S[np.where(I==1)]   
        J = np.copy(I)         
        for i in range(len(infected)):
            Infection_indices = np.where((dist(S,infected[i])<Spreading_distance)*(I!=2)[:,np.newaxis])[0]
            I[Infection_indices] = 1 
        m = len(I[I==1])
        Times[np.where(J!=I)] = T_recover
    n=n-1

    S=random_walk(S)
    TT = np.copy(Times)>0
    Times[Times>0] = Times[Times>0]-1
    TT2 = Times>0
    
    J = np.copy(I)
    I[TT!=TT2] = 2
    
    Times_immune[J!=I]=T_immune
    TT = np.copy(Times_immune)>0
    Times_immune[Times_immune>0] = Times_immune[Times_immune>0]-1
    TT2 = Times_immune>0
    I[TT!=TT2] = 0
    
    
    
    
    

anim = camera.animate() # repeat=False
P2.legend(['Healthy','Infected', 'Immune'], loc='upper right', 
          bbox_to_anchor=(0.8, 1.1), ncol=3)
P1.legend(['Total number infected'], loc='upper right', 
          bbox_to_anchor=(0.7, 1.1), ncol=3)

# anim.save('flatten_curve3.html')

