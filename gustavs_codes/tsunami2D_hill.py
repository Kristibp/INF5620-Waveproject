from wave2D_du0 import *

#
#     Profiling: python -m cProfile tsunami2D_hill.py
#

#import cProfile
#cProfile.run('foo()')

class Simulations:
    def __init__(self, option=None, version='scalar'):
        self.option = option
        self.version = version
        
    def simulate(self, option):   
        if option == 'gaussian_bottom':
            self.solver = self.simulate_gaussian_bottom()
        elif option == 'cosine_hat':
            self.solver = self.simulate_cosine_hat()
        elif option == 'box':
            self.solver = self.simulate_box()
            
    def simulate_gaussian_bottom(self):
        """
        Simulate an earthquake-generated tsunami over a
        elliptic gaussian-shaped bottom.
        """
        
        case = 2
        
        if case == 1:
            
            #     Length and with of simulated region
            Lx = 800.0     # Units: (m)
            Ly = 430.0
            Nx = 80         # Grid points in x direction
            Ny = 40         # Grid points in y direction
            T = 10.0         # Total simulated time (s)
            dt = 0.1        # Time step (s)
            
            g = 9.81        # Acceleration of gravity (m/s^2)
            
            I0 = 0.0
            Ia = 140.0
            Im = 0.0
            Is = 40.2
            B0 = -500.0
            Ba = 450.0
            Bmx = 0.3*Lx
            Bmy = 0.5*Ly
            Bs = 100.0
            bxy = 0.8
            
            zmin = -100.0
            zmax = 150.0
        
        elif case == 2:
                        
            #     Length and with of simulated region
            Lx = 800.0     # Units: (m)
            Ly = 430.0
            Nx = 80         # Grid points in x direction
            Ny = 40         # Grid points in y direction
            T = 10.0         # Total simulated time (s)
            dt = 0.1        # Time step (s)
            
            g = 9.81        # Acceleration of gravity (m/s^2)
            
            I0 = 0.0
            Ia = 140.0
            Im = 0.0
            Is = 40.2
            B0 = -500.0
            Ba = 250.0
            Bmx = 0.3*Lx
            Bmy = 0.5*Ly
            Bs = 50.0
            bxy = 2.0
            
            zmin = -500.0#-50.0
            zmax = 150.0


        #     Initial water surface with tsunami
        def I(x, y):
            return I0 + Ia*exp(-((x-Im)/Is)**2)
    
        #     Shape of sea bottom
        def bottom(x, y):
            return B0 + Ba*exp(-((x-Bmx)/Bs)**2 \
                                    -((y-Bmy)/(bxy*Bs))**2)
    
        def q(x, y):
            return -g*bottom(x, y)
    
        def plot_u(u, x, xv, y, yv, t, n):
    
            #surf(xv, yv, u, show=True, zlim=[-2.4, 2.4])
            #surf(xv, yv, u, zlim=[-50.0, 300.0])
            #axis([0.0,400.0,0.0,430.0,-500.0,300.0])
            #show()
            surf(xv, yv, u, zlim=[zmin, zmax])
            #hold('on')
            
            #print 'xvshape = ',xv.shape[0],',yvshape = ', yv.shape[1]
            b_a = zeros((xv.shape[0],yv.shape[1]))
            b_a[:,:] = bottom(xv, yv)
            #print 'b_a.shape=',b_a.shape 
            #surf(xv, yv, b_a)
            mesh(xv, yv, b_a, zlim=[zmin, zmax])
            #axis([0.0,400.0,0.0,430.0,-500.0,300.0])
            #show()
            #time.sleep(1.0)
            filename = 'tmp_gaussian2%04d.png' % n
            savefig(filename)
        
        #     Construct problem object
        problem = Problem(I=I, V=None, f=None, q=q, b=0.0, Lx=Lx, 
                          Ly=Ly, T=T)

        #     Construct solver object
        solver = Solver(problem=problem, Nx=Nx, Ny=Ny,
                        dt=dt, user_action=plot_u, 
                        version=self.version)
        #     Solve the PDE
        solver.solve()   

    def simulate_cosine_hat(self):
        """
        Simulate an earthquake-generated tsunami over 
        a sea bottom with a cosine hat hill.
        """

        Lx = 1000.0
        Ly = 430.0
        Nx = 50
        Ny = 40
        T = 14.0
        dt = 0.15
        
        g = 9.81        # Acceleration of gravity (m/s^2)
            
        I0 = 0.0
        Ia = 140.0
        Im = 0.0
        Is = 40.2
        
        B0 = -300.0
        Ba = 200.0
        Bmx = 0.3*Lx
        Bmy = 0.5*Ly
        Bs = 200.0
            
        zmin = -350.0
        zmax = 150.0
        
        #import time               # Measure CPU time
        t0 = time.clock()
        
        #     Initial water surface with tsunami
        def I(x, y):
            return I0 + Ia*exp(-((x-Im)/Is)**2)
        
        #     Shape of sea bottom
        def bottom(x, y):
            if type(x) == float64:
                if sqrt((x-Bmx)**2 + (y-Bmy)**2) <= Bs:
                    b_value = B0 + Ba*cos(0.5*pi*(x-Bmx)/Bs) \
                        *cos(0.5*pi*(y-Bmy)/Bs)
                else:
                    b_value = B0
            else:
                b_value = zeros((x.shape[0], y.shape[1]))
                for i in range(0, x.shape[0]):
                    for j in range(0, y.shape[1]):
                        xx = x[i,0]
                        yy = y[0,j]
                        if sqrt((xx-Bmx)**2 + (yy-Bmy)**2) <= Bs:
                            b_value[i,j] = B0 + Ba*cos(0.5*pi*(xx-Bmx)/Bs) \
                                *cos(0.5*pi*(yy-Bmy)/Bs)
                        else:
                            b_value[i,j] = B0
            return b_value
    
        def q(x, y):
            return -g*bottom(x, y)

        def plot_u(u, x, xv, y, yv, t, n):
    
            #surf(xv, yv, u, show=True, zlim=[-2.4, 2.4])
            #surf(xv, yv, u, zlim=[-50.0, 300.0])
            #axis([0.0,400.0,0.0,430.0,-500.0,300.0])
            #show()
            surf(xv, yv, u, zlim=[zmin, zmax])
            #hold('on')
            #1 + 1 
            #print 'xvshape = ',xv.shape[0],',yvshape = ', yv.shape[1]
            #b_a = zeros((xv.shape[0],yv.shape[1]))
            #print 'xv=',xv,',yv=',yv
            #b_a[:,:] = bottom(xv, yv)
            #print 'b_a=',b_a
            #print 'shape(b_a)=',shape(b_a)
            #print 'b_a.shape=',b_a.shape 
            #surf(xv, yv, b_a)
            #mesh(xv, yv, b_a)
            #mesh(xv, yv, b_a, zlim=[zmin, zmax])
            #axis([0.0,400.0,0.0,430.0,-500.0,300.0])
            #show()
            #time.sleep(5.0)
            #filename = 'tmp_cos_hat%04d.png' % n
            #savefig(filename)
        
        #     Construct problem object
        problem = Problem(I=I, V=None, f=None, q=q, b=0.0, Lx=Lx, 
                          Ly=Ly, T=T)

        #     Construct solver object
        solver = Solver(problem=problem, Nx=Nx, Ny=Ny,
                        dt=dt, user_action=plot_u, 
                        version=self.version)
        #     Solve the PDE
        solver.solve()   
    
        t1 = time.clock()
        print 'used time: ', t1-t0 


    def simulate_box(self):
        """
        Simulate an earthquake-generated tsunami over
        a sea bottom with a box-shaped hill.
        """
        
        Lx = 1000.0
        Ly = 430.0
        Nx = 50
        Ny = 40
        T = 10.0
        dt = 0.15
        
        g = 9.81        # Acceleration of gravity (m/s^2)
            
        I0 = 0.0
        Ia = 140.0
        Im = 0.0
        Is = 40.2
        
        B0 = -300.0
        Ba = 200.0
        Bmx = 0.3*Lx
        Bmy = 0.5*Ly
        Bs = 70.0
        b = 1.3
            
        zmin = -320.0
        zmax = 150.0
        
        t0 = time.clock()

        #     Initial water surface with tsunami
        def I(x, y):
            return I0 + Ia*exp(-((x-Im)/Is)**2)
        
        #     shape of sea bottom
        def bottom(x, y):
            if type(x) == float64:
                if (x >= Bmx-Bs) and (x <= Bmx+Bs) \
                        and (y >= Bmy-b*Bs) and (y <= Bmy+b*Bs):
                    b_value = B0 + Ba
                else:
                    b_value = B0
            else:
                b_value = zeros((x.shape[0], y.shape[1]))
                for i in range(0, x.shape[0]):
                    for j in range(0, y.shape[1]):
                        xx = x[i,0]
                        yy = y[0,j]
                        if (xx >= Bmx-Bs) and (xx <= Bmx+Bs) \
                                and (yy >= Bmy-b*Bs) and (yy <= Bmy+b*Bs):
                            b_value[i,j] = B0 + Ba
                        else:
                            b_value[i,j] = B0
            return b_value
        
        def q(x, y):
            return -g*bottom(x, y)


        def plot_u(u, x, xv, y, yv, t, n):
            
            #surf(xv, yv, u, show=True, zlim=[-2.4, 2.4])
            #surf(xv, yv, u, zlim=[-50.0, 300.0])
            #axis([0.0,400.0,0.0,430.0,-500.0,300.0])
            #show()
            surf(xv, yv, u, zlim=[zmin, zmax])
            #show()
            #hold('on')
            #1 + 1 
            #print 'xvshape = ',xv.shape[0],',yvshape = ', yv.shape[1]
            
            b_a = zeros((xv.shape[0],yv.shape[1]))
            #print 'xv=',xv,',yv=',yv
            b_a[:,:] = bottom(xv, yv)
            #print 'b_a=',b_a
            #print 'shape(b_a)=',shape(b_a)
            #print 'b_a.shape=',b_a.shape 
            #surf(xv, yv, b_a)
            #mesh(xv, yv, b_a)
            mesh(xv, yv, b_a, zlim=[zmin, zmax])
            
            #axis([0.0,400.0,0.0,430.0,-500.0,300.0])
            #show()
            #time.sleep(5.0)
            filename = 'tmp_box%04d.png' % n
            savefig(filename)
        

        #     Construct problem object
        problem = Problem(I=I, V=None, f=None, q=q, b=0.0, Lx=Lx, 
                          Ly=Ly, T=T)
        
        #     Construct solver object
        solver = Solver(problem=problem, Nx=Nx, Ny=Ny,
                        dt=dt, user_action=plot_u, 
                        version=self.version)
        
        #     Solve the PDE
        solver.solve()   
        
        t1 = time.clock()
        print 'used time: ', t1-t0 



def main(): 
    sim = Simulations(version='vectorized')
    sim.simulate('box')

if __name__ == '__main__':
    main()
