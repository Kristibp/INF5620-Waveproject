import scitools.std as st
from wave2D_du0 import *

#
#     Profiling: python -m cProfile tsunami2D_hill.py
#

#import cProfile
#cProfile.run('foo()')

class Simulations:
    def __init__(self, option=None, version='scalar', \
                     q_average='arithmetic'):
        self.option = option
        self.version = version
        self.q_average = q_average
        
    def simulate(self, option):   
        if option == 'gaussian':
            self.solver = self.simulate_gaussian()
        elif option == 'cosine_hat':
            self.solver = self.simulate_cosine_hat()
        elif option == 'box':
            self.solver = self.simulate_box()
        elif option == 'gaussian_1d':
            self.solver = self.simulate_gaussian_1d()
        elif option == 'gaussian_1d_slit':
            self.solver = self.simulate_gaussian_1d_slit()
        elif option == 'box_1d':
            self.solver = self.simulate_box_1d()
        elif option == 'box_1d_slit':
            self.solver = self.simulate_box_1d_slit()
        elif option == 'smooth_step':
            self.solver = self.simulate_smooth_step()
        elif option == 'exp_pol_1d':
            self.solver = self.simulate_exp_pol_1d()

    
    def simulate_gaussian(self):
        """
        Simulate an earthquake-generated tsunami over a
        elliptic gaussian-shaped bottom.
        """
        
        case = 1
        
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
            Ba = 400.0
            Bmx = 0.5*Lx
            Bmy = 0.8*Ly
            Bs = 200.0
            bxy = 1.8
            
            zmin = -500.0
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
            Ba =400.0
            Bmx = 0.3*Lx
            Bmy = 0.2*Ly
            Bs = 50.0
            bxy = 1.5
            
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
            st.surf(xv, yv, u, zlim=[zmin, zmax])
            #hold('on')
            
            #print 'xvshape = ',xv.shape[0],',yvshape = ', yv.shape[1]
            b_a = zeros((xv.shape[0],yv.shape[1]))
            b_a[:,:] = bottom(xv, yv)
            #print 'b_a.shape=',b_a.shape 
            #surf(xv, yv, b_a)
            st.mesh(xv, yv, b_a, zlim=[zmin, zmax])
            #axis([0.0,400.0,0.0,430.0,-500.0,300.0])
            #show()
            #time.sleep(1.0)
            filename = 'tmp_gaussian2%04d.png' % n
            st.savefig(filename)
        
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
        Nx = 80
        Ny = 40
        T = 14.0
        dt = 0.1
        
        g = 9.81        # Acceleration of gravity (m/s^2)
            
        I0 = 0.0
        Ia = 140.0
        Im = 0.0
        Is = 40.2
        
        B0 = -300.0
        Ba = 250.0
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
            st.surf(xv, yv, u, zlim=[zmin, zmax])
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
        T = 14.0
        dt = 0.15
        
        g = 9.81        # Acceleration of gravity (m/s^2)
        
        I0 = 0.0
        Ia = 140.0
        Im = 0.0
        Is = 40.2
        
        B0 = -300.0
        Ba = 400.0
        Bmx = 0.3*Lx
        Bmy = 0.5*Ly
        Bs = 35.0
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
            st.surf(xv, yv, u, zlim=[zmin, zmax])
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
            st.mesh(xv, yv, b_a, zlim=[zmin, zmax])
            
            #title('Bs = 40')
            #axis([0.0,400.0,0.0,430.0,-500.0,300.0])
            #show()
            #time.sleep(5.0)
            #filename = 'tmp_box%04d.png' % n
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


            
    def simulate_box_1d(self):
        """
        Simulate an one-dimensional earthquake-generated 
        tsunami over a sea bottom with a box-shaped hill.
        """
        
        Lx = 800.0
        Ly = 10.0
        Nx = 200
        Ny = 3
        T = 14.0
        dt = 0.04
        
        g = 9.81        # Acceleration of gravity (m/s^2)
        
        I0 = 0.0
        Ia = 140.0
        Im = 0.0
        Is = 40.2
        
        B0 = -300.0
        Ba = 275.0
        Bmx = 0.75*Lx
        Bs = 0.25*Lx
            
        zmin = -320.0
        zmax = 150.0
        
        t0 = time.clock()

        #     Initial water surface with tsunami
        def I(x, y):
            return I0 + Ia*exp(-((x-Im)/Is)**2)
        
        #     shape of sea bottom
        def bottom(x, y):
            if type(x) == float64:
                if (x >= Bmx-Bs) and (x <= Bmx+Bs):
                    b_value = B0 + Ba
                else:
                    b_value = B0
            else:
                b_value = zeros((x.shape[0], y.shape[1]))
                for i in range(0, x.shape[0]):
                    for j in range(0, y.shape[1]):
                        xx = x[i,0]
                        yy = y[0,j]
                        if (xx >= Bmx-Bs) and (xx <= Bmx+Bs):
                            b_value[i,j] = B0 + Ba
                        else:
                            b_value[i,j] = B0
            return b_value
        
        def q(x, y):
            return -g*bottom(x, y)


        def plot_u(u, x, xv, y, yv, t, n):
          
            b_a = zeros((xv.shape[0],yv.shape[1]))
            b_a[:,:] = bottom(xv, yv)
            st.plot(xv[:,0], u[:,0], '-', xv[:,0], b_a[:,0], '--', ylim=[zmin, zmax])
            #show()
            #hold('on')
            
            #axis([0.0,400.0,0.0,430.0,-500.0,300.0])
            #show()
            #time.sleep(5.0)
            filename = 'tmp_1d_box%04d.png' % n
            st.savefig(filename)
        

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


    def simulate_smooth_step(self):
        """
        Simulate an one-dimensional earthquake-generated 
        tsunami over a sea bottom with a box-shaped hill.
        """
        
        Lx = 800.0
        Ly = 10.0
        Nx = 200
        Ny = 3
        T = 14.0
        dt = 0.04
        
        g = 9.81        # Acceleration of gravity (m/s^2)
        
        I0 = 0.0
        Ia = 140.0
        Im = 0.0
        Is = 40.2
        
        B0 = -300.0
        Ba = 275.0
        Bmx = 0.5*Lx
        Bs = 0.005*Lx
            
        zmin = -320.0
        zmax = 150.0
        
        t0 = time.clock()

        #     Initial water surface with tsunami
        def I(x, y):
            return I0 + Ia*exp(-((x-Im)/Is)**2)
      
        #     Shape of sea bottom
        def bottom(x, y):
            return B0 + Ba/(1.0+exp(-(x-Bmx)/Bs))
        
        def q(x, y):
            return -g*bottom(x, y)


        def plot_u(u, x, xv, y, yv, t, n):
          
            b_a = zeros((xv.shape[0],yv.shape[1]))
            b_a[:,:] = bottom(xv, yv)
            st.plot(xv[:,0], u[:,0], '-', xv[:,0], b_a[:,0], '--', ylim=[zmin, zmax])
            #show()
            #hold('on')
            
            #axis([0.0,400.0,0.0,430.0,-500.0,300.0])
            #show()
            #time.sleep(5.0)
            filename = 'tmp_1d_smooth_step%04d.png' % n
            st.savefig(filename)
        

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


                
    def simulate_gaussian_1d(self):
        """
        Simulate an one-dimensional earthquake-generated 
        tsunami over a sea bottom with a gaussian shaped hill.
        """
        
        Lx = 800.0
        Ly = 10.0
        Nx = 200
        Ny = 3
        T = 14.0
        dt = 0.04
        
        g = 9.81        # Acceleration of gravity (m/s^2)
        
        I0 = 0.0
        Ia = 140.0
        Im = 0.0
        Is = 40.2
        
        B0 = -300.0
        Ba = 275.0
        Bmx = 0.5*Lx
        Bs = 35.0
            
        zmin = -320.0
        zmax = 150.0
        
        t0 = time.clock()

        #     Initial water surface with tsunami
        def I(x, y):
            return I0 + Ia*exp(-((x-Im)/Is)**2)

    
        #     Shape of sea bottom
        def bottom(x, y):
            return B0 + Ba*exp(-((x-Bmx)/Bs)**2)
        
        def q(x, y):
            return -g*bottom(x, y)


        def plot_u(u, x, xv, y, yv, t, n):
          
            b_a = zeros((xv.shape[0],yv.shape[1]))
            b_a[:,:] = bottom(xv, yv)
            st.plot(xv[:,0], u[:,0], '-', xv[:,0], b_a[:,0], '--', ylim=[zmin, zmax])
            #show()
            #hold('on')
            
            #axis([0.0,400.0,0.0,430.0,-500.0,300.0])
            #show()
            #time.sleep(5.0)
            filename = 'tmp_1d_gaussian%04d.png' % n
            st.savefig(filename)
        

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



    def simulate_exp_pol_1d(self):
        """
        Simulate an one-dimensional earthquake-generated 
        tsunami over a sea bottom with a box-shaped hill.
        """
        
        Lx = 800.0
        Ly = 10.0
        Nx = 100
        Ny = 3
        T = 17.0
        dt = 0.05
        
        g = 9.81        # Acceleration of gravity (m/s^2)
        
        I0 = 0.0
        Ia = 140.0
        Im = 0.0
        Is = 40.2
        
        B0 = -300.0
        Ba = 200.0
        Bmx = 0.5*Lx
        Bs = 0.1*Lx
        
        zmin = -320.0
        zmax = 150.0
        
        t0 = time.clock()
        
        #     Initial water surface with tsunami
        def I(x, y):
            return I0 + Ia*exp(-((x-Im)/Is)**2)
        
        #     Shape of sea bottom
        def bottom(x, y):
            xp = (x-Bmx)/Bs
            return B0 + Ba*(1 + 1.7*xp**2)*exp(-xp**2)
        
        def bottom2(x, y):
            return B0 + 1.15*Ba/( (1.0+exp(-(x-290.0)/(0.025*Lx))) \
                                *(1.0+exp((x-510.0)/(0.025*Lx))) )

        #     shape of sea bottom
        # def bottom3(x, y):
        #     if type(x) == float64:
        #         if ((x >= Bmx-Bs) and (x <= Bmx-Bs2)) \
        #                 or ((x >= Bmx+Bs2))):
        #             b_value = B0 + Ba
        #         else:
        #             b_value = B0
        #     else:
        #         b_value = zeros((x.shape[0], y.shape[1]))
        #         for i in range(0, x.shape[0]):
        #             for j in range(0, y.shape[1]):
        #                 xx = x[i,0]
        #                 yy = y[0,j]
        #                 if (xx >= Bmx-Bs) and (xx <= Bmx+Bs):
        #                     b_value[i,j] = B0 + Ba
        #                 else:
        #                     b_value[i,j] = B0
        #     return b_value
        

        def q(x, y):
            return -g*bottom3(x, y)
        
        def plot_u(u, x, xv, y, yv, t, n):
          
            b_a = zeros((xv.shape[0],yv.shape[1]))
            b_a[:,:] = bottom(xv, yv)
            b_a2 = zeros((xv.shape[0],yv.shape[1]))
            b_a2[:,:] = bottom2(xv, yv)
            #b_a3 = zeros((xv.shape[0],yv.shape[1]))
            #b_a3[:,:] = bottom3(xv, yv)
            st.plot(xv[:,0], u[:,0], '-', xv[:,0], b_a[:,0], \
                     xv[:,0], b_a2[:,0], '--', \
                     ylim=[zmin, zmax])
            #show()
            #hold('on')
            
            #axis([0.0,400.0,0.0,430.0,-500.0,300.0])
            #show()
            #time.sleep(5.0)
            filename = 'tmp_1d_exp_pol%04d.png' % n
            st.savefig(filename)
            

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


    def simulate_gaussian_1d_slit(self):
        """
        Simulate an one-dimensional earthquake-generated 
        tsunami over a sea bottom with a gaussian shaped hill with a gap.
        """
        
        Lx = 800.0
        Ly = 10.0
        Nx = 200
        Ny = 3
        T = 14.0
        dt = 0.04
        
        g = 9.81        # Acceleration of gravity (m/s^2)
        
        I0 = 0.0
        Ia = 140.0
        Im = 0.0
        Is = 40.2
        
        B0 = -300.0
        Ba = 275.0
        Bmx = 0.5*Lx
        Bs = 35.0
            
        zmin = -320.0
        zmax = 150.0
        
        t0 = time.clock()

        #epsilon = size of the gap
        #eps = 0
        eps = 10

        #     Initial water surface with tsunami
        def I(x, y):
            return I0 + Ia*exp(-((x-Im)/Is)**2)

        def ordinaryBottom(x,y):
            return B0 + Ba*exp(-((x-Bmx)/Bs)**2)

        def hole(x,y):
            return B0

        #     Shape of sea bottom
        def bottom(x, y):
            if type(x) == float64:
                if (x >= Bmx-10) and (x <= Bmx+10):
                    b_value = hole(x,y)
                else:
                    b_value = ordinaryBottom(x,y)
            else:
                b_value = zeros((x.shape[0], y.shape[1]))
                for i in range(0, x.shape[0]):
                    for j in range(0, y.shape[1]):
                        xx = x[i,0]
                        yy = y[0,j]
                        if (xx >= Bmx-eps) and (xx <= Bmx+eps):
                            b_value[i,j] = hole(xx,yy)
                        else:
                            b_value[i,j] = ordinaryBottom(xx,yy)
            return b_value
        
        def q(x, y):
            return -g*bottom(x, y)


        def plot_u(u, x, xv, y, yv, t, n):
          
            b_a = zeros((xv.shape[0],yv.shape[1]))
            b_a[:,:] = bottom(xv, yv)
            st.plot(xv[:,0], u[:,0], '-', xv[:,0], b_a[:,0], '--', ylim=[zmin, zmax])
            #show()
            #hold('on')
            
            #axis([0.0,400.0,0.0,430.0,-500.0,300.0])
            #show()
            #time.sleep(5.0)
            filename = 'tmp_1d_gaussian_slit%04d.png' % n
            st.savefig(filename)
        

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

    def simulate_box_1d_slit(self):
        """
        Simulate an one-dimensional earthquake-generated 
        tsunami over a sea bottom with a box-shaped hill with a gap.
        """
        
        Lx = 800.0
        Ly = 10.0
        Nx = 200
        Ny = 3
        T = 14.0
        dt = 0.04
        
        g = 9.81        # Acceleration of gravity (m/s^2)
        
        I0 = 0.0
        Ia = 140.0
        Im = 0.0
        Is = 40.2
        
        B0 = -300.0
        Ba = 275.0
        Bmx = 0.5*Lx
        Bs = 0.10*Lx
            
        zmin = -320.0
        zmax = 150.0
        
        t0 = time.clock()

        #epsilon = size of the gap
        #eps = 0
        eps = 10

        #     Initial water surface with tsunami
        def I(x, y):
            return I0 + Ia*exp(-((x-Im)/Is)**2)
        
        #     shape of sea bottom

        def ordinaryBottom(x,y):
            if (x >= Bmx-Bs) and (x <= Bmx+Bs):
                b_value = B0 + Ba
            else:
                b_value = B0
            return b_value
                
        def hole(x,y):
            return B0

        #     Shape of sea bottom
        def bottom(x, y):
            if type(x) == float64:
                if (x >= Bmx-eps) and (x <= Bmx+eps):
                    b_value = hole(x,y)
                else:
                    b_value = ordinaryBottom(x,y)
            else:
                b_value = zeros((x.shape[0], y.shape[1]))
                for i in range(0, x.shape[0]):
                    for j in range(0, y.shape[1]):
                        xx = x[i,0]
                        yy = y[0,j]
                        if (xx >= Bmx-eps) and (xx <= Bmx+eps)
                            b_value[i,j] = hole(xx,yy)
                        else:
                            b_value[i,j] = ordinaryBottom(xx,yy)
            return b_value
        
        def q(x, y):
            return -g*bottom(x, y)


        def plot_u(u, x, xv, y, yv, t, n):
          
            b_a = zeros((xv.shape[0],yv.shape[1]))
            b_a[:,:] = bottom(xv, yv)
            st.plot(xv[:,0], u[:,0], '-', xv[:,0], b_a[:,0], '--', ylim=[zmin, zmax])
            #show()
            #hold('on')
            
            #axis([0.0,400.0,0.0,430.0,-500.0,300.0])
            #show()
            #time.sleep(5.0)
            filename = 'tmp_1d_box_slit%04d.png' % n
            st.savefig(filename)
        

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
    sim = Simulations(version='vectorized',\
                          q_average='arithmetic')
    #sim.simulate('gaussian_1d')
    #sim.simulate('box_1d')
    #sim.simulate('smooth_step')
    #sim.simulate('exp_pol_1d')
    #sim.simulate('gaussian_1d_slit')
    sim.simulate('box_1d_slit')

if __name__ == '__main__':
    main()
