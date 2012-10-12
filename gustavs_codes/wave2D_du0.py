#!/usr/bin/env python

"""
A finite difference scheme is used to solve the 2D wave
equation 

u_tt + b*u_t = (q*u_x)_x + (q*u_y)_y + f(x, y, t) 

on the domain (0, Lx)*(0, Ly). The boundary condition is 
of Neumann type, with du/dn = 0, where n is the normal vector 
on the boundary. The initial conditions are u(x, y, 0) = 
I(x, y) and u_t(x, y, 0) = V(x, y).

Nx and Ny are the total number of mesh cells in the x and y
directions. The mesh points are numbered as (0,0), (1,0), (2,0),
..., (Nx,0), (0,1), (1,1), ..., (Nx, Ny).

dt is the time step. If dt<=0, an optimal time step is used.
T is the stop time for the simulation.

I, V, f are functions: I(x,y), V(x,y), f(x,y,t). V and f
can be specified as None or 0, resulting in V=0 and f=0.

user_action: function of (u, x, y, t, n) called at each time
level (x and y are one-dimensional coordinate vectors).
This function allows the calling code to plot the solution,
compute errors, etc.

"""

import time
from scitools.std import *
from math import log
from numpy import *


class Problem:
    def __init__(self, I, V, f, q, b, Lx, Ly, T):
        """
        I = I(x, y) = initial condition u(x, y, t=0)
        V = V(x, y) = initial condition u_t(x, y, t=0)
        f = f(x, y, t) 
        q = q(x, y) 
        b = constant
        Lx = width of domain in the coordinate x
        Ly = width of domain in the coordinate y
        T = stopping time for simulation
        """
        self.I = I
        self.V = V
        self.f = f
        self.q = q
        self.b = b
        self.Lx = Lx
        self.Ly = Ly
        self.T = T
        
        
class Solver:
    
    def __init__(self, problem, Nx, Ny, dt=0.01, \
                     user_action=None, version='scalar', \
                     dt_safety_factor=0.9):
        """
        Nx = number of grid points in space coordinate x
        Ny = number of grid points in space coordinate y
        dt = time step
        user_action = action to do after every time step
        version = computational option (serial, vectorized,
                  fortran, etc.)
        """
        self.problem = problem
        self.Nx = Nx
        self.Ny = Ny
        self.dt = dt
        self.user_action = user_action
        self.version = version
        self.dt_safety_factor = dt_safety_factor
        
    def solve(self):
        self.dt, self.time_used, self.x, self.y, \
            self.xv, self.yv, self.u_last, \
            self.t_last = self.solve_2d_wave()
        

    def solve_2d_wave(self):
        """
        Solve the 2D wave PDE using a finite difference 
        method.
        """
        
        I = self.problem.I
        V = self.problem.V
        f = self.problem.f
        q = self.problem.q
        b = self.problem.b
        Lx = self.problem.Lx
        Ly = self.problem.Ly
        T = self.problem.T
        Nx = self.Nx
        Ny = self.Ny
        dt = self.dt
        user_action = self.user_action
        version = self.version
        dt_safety_factor = self.dt_safety_factor
        
        if version == 'cython':
            try:
                #import pyximport; pyximport.install()
                import wave2D_u0_loop_cy as compiled_loops
            except ImportError, e:
                print 'No module wave2D_u0_loop_cy. Run make_wave2D.sh!'
            print e
            sys.exit(1)
        elif version == 'f77':
            try:
                import wave2D_u0_loop_f77 as compiled_loops
            except ImportError:
                print 'No module wave2D_u0_loop_f77. Run make_wave2D.sh!'
            sys.exit(1)
        elif version == 'c_f2py':
            try:
                import wave2D_u0_loop_c_f2py as compiled_loops
            except ImportError:
                print 'No module wave2D_u0_loop_c_f2py. Run make_wave2D.sh!'
            sys.exit(1)
        elif version == 'c_cy':
            try:
                import wave2D_u0_loop_c_cy as compiled_loops
            except ImportError, e:
                print 'No module wave2D_u0_loop_c_cy. Run make_wave2D.sh!'
            print e
            sys.exit(1)
            
        if version == 'scalar':
            advance = self.advance_scalar
        elif version == 'vectorized':
            advance = self.advance_vectorized
        elif version == 'c_cy':
            advance = compiled_loops.advance_cwrap
        else:
            advance = compiled_loops.advance
                
        import time               # Measure CPU time
        t0 = time.clock()
        
        x = linspace(0, Lx, Nx+1) # Mesh points in x dir
        y = linspace(0, Ly, Ny+1) # Mesh points in y dir
        dx = x[1] - x[0]
        dy = y[1] - y[0]
        
        xv = x[:,newaxis]         # For vectorized function 
        yv = y[newaxis,:]         # evaluations
    
        N = int(round(T/float(dt)))
        t = linspace(0, N*dt, N+1)          # Mesh points in time
        #     Help variables
        Cx2 = 0.5*(dt/dx)**2;  Cy2 = 0.5*(dt/dy)**2  
        dt2 = dt**2
        tb = 0.5*dt*b
        
        #     Allow f and V to be None or 0
        if f is None or f == 0:
            f = lambda x, y, t: 0.0
            #f = (lambda x, y, t: 0) if version == 'scalar' else \
            #    lambda x, y, t: zeros((xv.shape[0], yv.shape[1]))
                # or simpler: x*y*0
        if V is None or V == 0:
            V = lambda x, y: 0.0
            #V = (lambda x, y: 0) if version == 'scalar' else \
            #    lambda x, y: zeros((xv.shape[0], yv.shape[1]))
        
        order = 'Fortran' if version == 'f77' else 'C'
        u   = zeros((Nx+1,Ny+1), order=order)   # Solution array
        u_1 = zeros((Nx+1,Ny+1), order=order)   # Solution at t-dt
        u_2 = zeros((Nx+1,Ny+1), order=order)   # Solution at t-2*dt
        f_a = zeros((Nx+1,Ny+1), order=order)   # For compiled loops
        V_a = zeros((Nx+1,Ny+1), order=order)
        q_a = zeros((Nx+1,Ny+1), order=order)
        
        V_a[:,:] = V(xv, yv)
        q_a[:,:] = q(xv, yv)

        max_q = float(q_a.max())

        stability_limit = (1/sqrt(max_q))*(1/sqrt(1/dx**2 + 1/dy**2))
        print 'stability limit dt = ', stability_limit
        if dt <= 0:               
            dt = dt_safety_factor*stability_limit
        elif dt > stability_limit:
            print 'error: dt=%g exceeds the stability limit %g' % \
                (dt, stability_limit)
        
        #     Set initial condition
        if version == 'scalar':
            for i in range(0, Nx+1):
                for j in range(0, Ny+1):
                    u_1[i,j] = I(x[i], y[j])
        else: # use vectorized version
            u_1[:,:] = I(xv, yv)
        
        if user_action is not None:
            user_action(u_1, x, xv, y, yv, t, 0)
            
        print '---'
        #     General evaluation of u[i,j] for the first step
        def u_step(i, j, ip1, im1, jp1, jm1, xip1, xim1, yjp1, yjm1):
            uij = u_1[i,j] + (1.0 - tb)*dt*V(x[i], y[j]) \
                + 0.5*Cx2*( (q(x[i], y[j]) + q(xip1, y[j])) \
                *(u_1[ip1,j] - u_1[i,j]) \
                - (q(xim1, y[j]) + q(x[i], y[j])) \
                *(u_1[i,j] - u_1[im1,j]) ) \
                + 0.5*Cy2*( (q(x[i], y[j]) + q(x[i], yjp1)) \
                *(u_1[i,jp1] - u_1[i,j]) \
                - (q(x[i], yjm1) + q(x[i], y[j])) \
                *(u_1[i,j] - u_1[i,jm1]) ) + 0.5*f(x[i], y[j], t[0])*dt2
        
            return uij
    
        n = 0
        if version == 'scalar':
            
            #     First step for interior points
            for i in range(1, Nx):
                for j in range(1, Ny):
                    ip1 = i+1
                    im1 = i-1
                    jp1 = j+1
                    jm1 = j-1
                    xip1 = x[i+1]
                    xim1 = x[i-1]
                    yjp1 = y[j+1]
                    yjm1 = y[j-1]
                    u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                                        xip1, xim1, yjp1, yjm1)
            
        else:  # use vectorized version
            f_a[:,:] = f(xv, yv, t[n])  # precompute, size as u
            
            #     First step for interior points
            u[1:-1,1:-1] = u_1[1:-1,1:-1] + (1.0 - tb)*dt*V_a[1:-1,1:-1] \
                + 0.5*Cx2*( (q_a[1:-1,1:-1] + q_a[2:,1:-1]) \
                *(u_1[2:,1:-1] - u_1[1:-1,1:-1]) \
                - (q_a[:-2,1:-1] + q_a[1:-1,1:-1]) \
                *(u_1[1:-1,1:-1] - u_1[:-2,1:-1]) ) \
                + 0.5*Cy2*( (q_a[1:-1,1:-1] + q_a[1:-1,2:]) \
                *(u_1[1:-1,2:] - u_1[1:-1,1:-1]) \
                - (q_a[1:-1,:-2] + q_a[1:-1,1:-1] ) \
                *(u_1[1:-1,1:-1] - u_1[1:-1,:-2]) ) + 0.5*dt2*f_a[1:-1,1:-1]  
        
    
        #
        #     The boundaries are always treated with scalar 
        #     summation.
        #    
        #     First step for boundary points        
        j = 0
        for i in range(1, Nx):
            ip1 = i+1
            im1 = i-1
            jp1 = j+1
            jm1 = jp1
            xip1 = x[i+1]
            xim1 = x[i-1]
            yjp1 = y[j+1]
            yjm1 = y[j] - dy
            u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                                xip1, xim1, yjp1, yjm1)
        
        j = Ny
        for i in range(1, Nx):
            ip1 = i+1
            im1 = i-1
            jm1 = j-1
            jp1 = jm1
            xip1 = x[i+1]
            xim1 = x[i-1]
            yjp1 = y[j] + dy
            yjm1 = y[j-1] 
            u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                                xip1, xim1, yjp1, yjm1)
        
        i = 0
        for j in range(1, Ny):
            ip1 = i+1
            im1 = ip1
            jp1 = j+1
            jm1 = j-1
            xip1 = x[i+1]
            xim1 = x[i] - dx
            yjp1 = y[j+1]
            yjm1 = y[j-1] 
            u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                                xip1, xim1, yjp1, yjm1)
            
        i = Nx
        for  j in range(1, Ny):
            im1 = i-1
            ip1 = im1
            jp1 = j+1
            jm1 = j-1
            xip1 = x[i] + dx
            xim1 = x[i-1] 
            yjp1 = y[j+1]
            yjm1 = y[j-1]
            u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                                xip1, xim1, yjp1, yjm1)

        #     First step for corner points
        i = 0 
        j = 0
        ip1 = i+1
        im1 = ip1
        jp1 = j+1
        jm1 = jp1
        xip1 = x[i+1] 
        xim1 = x[i] - dx 
        yjp1 = y[j+1]
        yjm1 = y[j] - dy
        u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                            xip1, xim1, yjp1, yjm1)

        i = Nx
        j = 0
        im1 = i-1
        ip1 = im1
        jp1 = j+1
        jm1 = jp1
        xip1 = x[i] + dx 
        xim1 = x[i-1] 
        yjp1 = y[j+1]
        yjm1 = y[j] - dy
        u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                            xip1, xim1, yjp1, yjm1)
        
        i = Nx
        j = Ny
        im1 = i-1
        ip1 = im1
        jm1 = j-1
        jp1 = jm1
        xip1 = x[i] + dx 
        xim1 = x[i-1] 
        yjp1 = y[j] + dy
        yjm1 = y[j-1] 
        u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                            xip1, xim1, yjp1, yjm1)

        i = 0
        j = Ny
        ip1 = i+1
        im1 = ip1
        jm1 = j-1
        jp1 = jm1
        xip1 = x[i+1] 
        xim1 = x[i] - dx 
        yjp1 = y[j] + dy
        yjm1 = y[j-1] 
        u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                            xip1, xim1, yjp1, yjm1)

        
        if user_action is not None:
            user_action(u, x, xv, y, yv, t, 1)
            #print '>> u=',u
        #print '###'
        
        u_2[:,:] = u_1; u_1[:,:] = u
    
        for n in range(1, N):
            if version == 'scalar':
                # use f(x,y,t) function
                u = advance(u, u_1, u_2, f, q, x, y, dx, dy, \
                                t, n, Cx2, Cy2, tb, dt2)
            else:
                f_a[:,:] = f(xv, yv, t[n])  # precompute, size as u
                u = advance(u, u_1, u_2, f_a, q_a, f, q, x, y, \
                                dx, dy, t, n, Cx2, Cy2, tb, dt2)

            if version == 'f77':
                for a in 'u', 'u_1', 'u_2', 'f_a':
                    if not isfortran(eval(a)):
                        print '%s: not Fortran storage!' % a

            if user_action is not None:
                if user_action(u, x, xv, y, yv, t, n+1):
                    break
            #print '==='
            u_2[:,:], u_1[:,:] = u_1, u
            
        t1 = time.clock()
        #     dt might be computed in this function, 
        #     so return the value
        return dt, t1 - t0, x, y, xv, yv, u, t[N]
    

    def advance_vectorized(self, u, u_1, u_2, f_a, q_a, f, q, x, y, \
                               dx, dy, t, n, Cx2, Cy2, tb, dt2):
        """
        Advance one time step in the finite difference scheme.
        This is a version with vector summation.
        """
        Nx = u.shape[0]-1;  Ny = u.shape[1]-1
        
        u[1:-1,1:-1] = ( 2.0*u_1[1:-1,1:-1] - (1.0-tb)*u_2[1:-1,1:-1] \
                            + Cx2*( (q_a[1:-1,1:-1] + q_a[2:,1:-1]) \
                                       *(u_1[2:,1:-1] - u_1[1:-1,1:-1]) \
                                       - (q_a[:-2,1:-1] + q_a[1:-1,1:-1]) \
                                        *(u_1[1:-1,1:-1] - u_1[:-2,1:-1]) ) \
                            + Cy2*( (q_a[1:-1,1:-1] + q_a[1:-1,2:]) \
                                        *(u_1[1:-1,2:] - u_1[1:-1,1:-1]) \
                                        - (q_a[1:-1,:-2] + q_a[1:-1,1:-1]) \
                                        *(u_1[1:-1,1:-1] - u_1[1:-1,:-2]) ) \
                            + f_a[1:-1,1:-1]*dt2 )*(1.0/(1.0+tb))
    
        u = self.advance_boundary(u, u_1, u_2, f, q, x, y, \
                                      dx, dy, t, n, Cx2, Cy2, tb, \
                                      dt2, Nx, Ny)
        
        return u
    
    def advance_scalar(self, u, u_1, u_2, f, q, x, y, dx, dy, \
                           t, n, Cx2, Cy2, tb, dt2):
        """
        Advance one time step in the finite difference scheme.
        This is a version with scalar summation.
        """

        def u_step(i, j, ip1, im1, jp1, jm1, xip1, xim1, yjp1, yjm1):
            uij = (2.0*u_1[i,j] - (1.0-tb)*u_2[i,j] \
                   + Cx2*((q(x[i], y[j]) + q(xip1, y[j])) \
                              *(u_1[ip1,j] - u_1[i,j]) \
                              - (q(xim1, y[j]) + q(x[i], y[j])) \
                              *(u_1[i,j] - u_1[im1,j])) \
                   + Cy2*((q(x[i], y[j]) + q(x[i], yjp1)) \
                              *(u_1[i,jp1] - u_1[i,j]) \
                              - (q(x[i], yjm1) + q(x[i], y[j])) \
                              *(u_1[i,j] - u_1[i,jm1])) \
                   + f(x[i], y[j], t[n])*dt2)/(1.0 + tb)
            return uij

        Nx = u.shape[0]-1;  Ny = u.shape[1]-1
        #     Loop over interior points
        for i in range(1, Nx):
            for j in range(1, Ny):
                ip1 = i+1
                im1 = i-1
                jp1 = j+1
                jm1 = j-1
                xip1 = x[i+1]
                xim1 = x[i-1]
                yjp1 = y[j+1]
                yjm1 = y[j-1]
                u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                                    xip1, xim1, yjp1, yjm1)
        
        u = self.advance_boundary(u, u_1, u_2, f, q, x, y, \
                             dx, dy, t, n, Cx2, Cy2, tb, \
                             dt2, Nx, Ny)
        
        return u


    
    def advance_boundary(self, u, u_1, u_2, f, q, x, y, \
                             dx, dy, t, n, Cx2, Cy2, tb, \
                             dt2, Nx, Ny):

        def u_step(i, j, ip1, im1, jp1, jm1, xip1, xim1, yjp1, yjm1):
            uij = (2.0*u_1[i,j] - (1.0-tb)*u_2[i,j] \
                   + Cx2*((q(x[i], y[j]) + q(xip1, y[j])) \
                              *(u_1[ip1,j] - u_1[i,j]) \
                              - (q(xim1, y[j]) + q(x[i], y[j])) \
                              *(u_1[i,j] - u_1[im1,j])) \
                   + Cy2*((q(x[i], y[j]) + q(x[i], yjp1)) \
                              *(u_1[i,jp1] - u_1[i,j]) \
                              - (q(x[i], yjm1) + q(x[i], y[j])) \
                              *(u_1[i,j] - u_1[i,jm1])) \
                   + f(x[i], y[j], t[n])*dt2)/(1.0 + tb)
            return uij

        
        #     Boundary points
        j = 0
        for i in range(1, Nx):
            ip1 = i+1
            im1 = i-1
            jp1 = j+1
            jm1 = jp1
            xip1 = x[i+1]
            xim1 = x[i-1]
            yjp1 = y[j+1]
            yjm1 = y[j] - dy
            u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                                xip1, xim1, yjp1, yjm1)

        j = Ny
        for i in range(1, Nx):
            ip1 = i+1
            im1 = i-1
            jm1 = j-1
            jp1 = jm1
            xip1 = x[i+1]
            xim1 = x[i-1]
            yjp1 = y[j] + dy
            yjm1 = y[j-1]
            u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                                xip1, xim1, yjp1, yjm1)

        i = 0
        for j in range(1, Ny):
            ip1 = i+1
            im1 = ip1
            jp1 = j+1
            jm1 = j-1
            xip1 = x[i+1]
            xim1 = x[i] - dx
            yjp1 = y[j+1]
            yjm1 = y[j-1]
            u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                                xip1, xim1, yjp1, yjm1)

    
        i = Nx
        for j in range(1, Ny):
            im1 = i-1
            ip1 = im1
            jp1 = j+1
            jm1 = j-1
            xip1 = x[i] + dx
            xim1 = x[i-1] 
            yjp1 = y[j+1]
            yjm1 = y[j-1]
            u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                                xip1, xim1, yjp1, yjm1)

        #     Corner points
        i = 0
        j = 0
        ip1 = i+1
        im1 = ip1
        jp1 = j+1
        jm1 = jp1
        xip1 = x[i+1] 
        xim1 = x[i] - dx 
        yjp1 = y[j+1]
        yjm1 = y[j] - dy
        u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                            xip1, xim1, yjp1, yjm1)
    
        i = Nx
        j = 0
        im1 = i-1
        ip1 = im1
        jp1 = j+1
        jm1 = jp1
        xip1 = x[i] + dx 
        xim1 = x[i-1]  
        yjp1 = y[j+1]
        yjm1 = y[j] - dy
        u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                            xip1, xim1, yjp1, yjm1)

        i = Nx
        j = Ny
        im1 = i-1
        ip1 = im1
        jm1 = j-1
        jp1 = jm1
        xip1 = x[i] + dx 
        xim1 = x[i-1]  
        yjp1 = y[j] + dy
        yjm1 = y[j-1] 
        u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                            xip1, xim1, yjp1, yjm1)

        i = 0
        j = Ny
        ip1 = i+1
        im1 = ip1
        jm1 = j-1
        jp1 = jm1
        xip1 = x[i+1]  
        xim1 = x[i] - dx  
        yjp1 = y[j] + dy
        yjm1 = y[j-1] 
        u[i,j] = u_step(i, j, ip1, im1, jp1, jm1, \
                            xip1, xim1, yjp1, yjm1)
        
        return u


class Tests:
    def __init__(self, option=None, version='scalar'):
        self.option = option
        self.version = version

    def do_test(self, option):
        if option == 'constant':
            self.solver = self.test_constant()
        elif option == '1d_plug_x':
            self.solver = self.test_1d_plug_x()
        elif option == '1d_plug_y':
            self.solver = self.test_1d_plug_y()
        elif option == '1d_gaussian_x':
            self.solver = self.test_1d_gaussian_x()
        elif option == '1d_gaussian_y':
            self.solver = self.test_1d_gaussian_y()
        elif option == 'standing_wave':
            self.solver = self.test_standing_wave()
        elif option == 'convergence_sw':
            self.solver = self.test_convergence_sw()
        elif option == 'convergence_sw2':
            self.solver = self.test_convergence_sw2()

    def test_constant(self):
        """
        Test with input that should give a constant solution 
        u(x, y, t) = u0. 
        """
        
        Lx = 10.0
        Ly = 10.0
        u0 = 2.0
        b = 3.3
        Nx = 5
        Ny = 5
        T = 0.4
        dt = 0.1
        
        def q(x, y):
            return x**2 + y**2

        def I(x, y):
            return u0
    
        def V(x, y):
            return 0.0 
            # print 'type(x)=',type(x)
            # if type(x) == float:
            #     vv = 0.0
            # elif type(x) == numpy.ndarray:
            #     print 'shape(x)[0]=',shape(x)[0]
            #     vv = zeros((shape(x)[0], shape(y)[1]))
            # else:
            #     print 'WARNING! Function is not defined'
            #     print 'for the input format.'
            #     sys.exit()
            # return vv
        
        def f(x, y, t):
            return 0.0
        
        def action(u, x, xv, y, yv, t, n):
            #print t
            print u
            
        #     Construct problem object
        problem = Problem(I=I, V=V, f=f, q=q, b=b, Lx=Lx, 
                          Ly=Ly, T=T)
        
        #     Construct solver object
        solver = Solver(problem=problem, Nx=Nx, Ny=Ny,
                        dt=dt, user_action=action, 
                        version=self.version)
        #     Solve the PDE
        solver.solve()
        

    def test_1d_plug_x(self):
        """
        Test the solver with a one-dimensional plug wave.
        """
        
        Lx = 1.0
        Ly = 0.5
        b = 0.0
        c = 0.5
        Nx = 50
        Ny = 3
        T = 2.0
        dt = 0.04
        print 'C=',c*dt/(Lx/Nx)
        
        def q(x, y):
            return c**2
        
        def f(x, y, t):
            return 0.0
            
        def V(x, y):
            return 0.0
        
        def I(x, y):
            if type(x) == float64:
                if (abs(x-0.5*Lx) > 0.1):
                    ivalue = 0.0
                else:
                    ivalue = 1.0
            else:
                ivalue = zeros((x.shape[0], y.shape[1]))
                for i in range(0, x.shape[0]):
                    if (abs(x[i,0]-0.5*Lx) > 0.1):
                        ivalue[i,:] = 0.0
                    else:
                        ivalue[i,:] = 1.0
            return ivalue
        
        def print_u(u, x, xv, y, yv, t, n):
            #print t
            print u

        def plot_u(u, x, xv, y, yv, t, n):
            if t[n] == 0:
                time.sleep(1)
            #mesh(x, y, u, title='t=%g' % t[n])
            plot(x, u[:,1], title='t=%g' % t[n], ylim=[-2.4, 2.4])
            filename = 'tmp_%04d.png' % n
            savefig(filename)
        
            time.sleep(0.4)
  
        #     Construct problem object
        problem = Problem(I=I, V=V, f=f, q=q, b=b, Lx=Lx, 
                          Ly=Ly, T=T)
        
        #     Construct solver object
        solver = Solver(problem=problem, Nx=Nx, Ny=Ny,
                        dt=dt, user_action=plot_u, 
                        version=self.version)
        #     Solve the PDE
        solver.solve()
    

    def test_1d_plug_y(self):
        """
        Test the solver with a one-dimensional plug wave.
        """
        
        Lx = 0.5
        Ly = 1.0
        b = 0.0
        c = 0.5
        Nx = 3
        Ny = 50
        T = 2.0
        dt = 0.04
        print 'C=',c*dt/(Ly/Ny)
        
        def q(x, y):
            return c**2
            
        def f(x, y, t):
            return 0.0
        
        def V(x, y):
            return 0.0

        def I(x, y):
            if type(y) == float64:
                if (abs(y-0.5*Ly) > 0.1):
                    ivalue = 0.0
                else:
                    ivalue = 1.0
            else:
                ivalue = zeros((x.shape[0], y.shape[1]))
                for j in range(0, y.shape[1]):
                    if (abs(y[0,j]-0.5*Ly) > 0.1):
                        ivalue[:,j] = 0.0
                    else:
                        ivalue[:,j] = 1.0
            return ivalue
        
        
        def print_u(u, x, xv, y, yv, t, n):
            #print t
            print u

        def plot_u(u, x, xv, y, yv, t, n):
            if t[n] == 0:
                time.sleep(1)
            #mesh(x, y, u, title='t=%g' % t[n])
            plot(y, u[1,:], title='t=%g' % t[n], ylim=[-2.4, 2.4])
            filename = 'tmp_%04d.png' % n
            savefig(filename)
            
            time.sleep(0.4)

        #     Construct problem object
        problem = Problem(I=I, V=V, f=f, q=q, b=b, Lx=Lx, 
                          Ly=Ly, T=T)
        
        #     Construct solver object
        solver = Solver(problem=problem, Nx=Nx, Ny=Ny,
                        dt=dt, user_action=plot_u, 
                        version=self.version)
        #     Solve the PDE
        solver.solve()


    def test_1d_gaussian_x(self):
        """
        Test the solver with a one-dimensional gaussian wave.
        """
        
        Lx = 1.0
        Ly = 0.5
        b = 0.0
        c = 0.5
        a = 500.0
        Nx = 100
        Ny = 3
        T = 2.4
        dt = 0.02
        print 'c=',c,',dt=',dt,',Nx=',Nx,',Lx=',Lx
        print 'C=',c*dt*Nx/Lx
        
        def q(x, y):
            return c**2
        
        def f(x, y, t):
            return 0.0
        
        def V(x, y):
            return 0.0
    
        def I(x, y):
            return exp(-a*(x-0.5*Lx)**2)
    
        def u_exact(x, t, a, Lx):
            return 0.5*exp(-a*(x-c*t-0.5*Lx)**2) \
                + 0.5*exp(-a*(x+c*t-0.5*Lx)**2)
            #return 0.5*exp(-a*(abs(x-0.5*Lx)-c*t)**2) \
            #    + 0.5*exp(-a*(abs(x-0.5*Lx)+c*t)**2)
    
        def print_u(u, x, xv, y, yv, t, n):
            #print t
            #print 'n = ', n
            print u[:,1]-u_exact(x, t[n], a, Lx)  
                
        def plot_u(u, x, xv, y, yv, t, n):
            if t[n] == 0:
                time.sleep(1)
            plot(x, u[:,1],'--',x, u_exact(x, t[n], a, Lx),'o')
            filename = 'tmp_%04d.png' % n
            savefig(filename)
            time.sleep(0.1)
        
        #     Construct problem object
        problem = Problem(I=I, V=V, f=f, q=q, b=b, Lx=Lx, 
                          Ly=Ly, T=T)
        
        #     Construct solver object
        solver = Solver(problem=problem, Nx=Nx, Ny=Ny,
                        dt=dt, user_action=plot_u, 
                        version=self.version)
        #     Solve the PDE
        solver.solve()
        
        
    def test_1d_gaussian_y(self):
        """
        Test the solver with a one-dimensional gaussian wave.
        """
        
        Lx = 0.5
        Ly = 1.0
        b = 0.0
        c = 0.5
        a = 500.0
        Nx = 3
        Ny = 100
        T = 2.4
        dt = 0.02
        print 'c=',c,',dt=',dt,',Ny=',Ny,',Ly=',Ly
        print 'C=',c*dt*Ny/Ly
        
        def q(x, y):
            return c**2
        
        def f(x, y, t):
            return 0.0
        
        def V(x, y):
            return 0.0
    
        def I(x, y):
            return exp(-a*(y-0.5*Ly)**2)
    
        def u_exact(y, t, a, Ly):
            return 0.5*exp(-a*(y-c*t-0.5*Ly)**2) \
                + 0.5*exp(-a*(y+c*t-0.5*Ly)**2)
    
        def print_u(u, x, xv, y, yv, t, n):
            #print t
            #print 'n = ', n
            print u[1,:]-u_exact(y, t[n], a, Ly)  #len(u) #len(u_exact(x, t[n], a, Lx))
        
        def plot_u(u, x, xv, y, yv, t, n):
            if t[n] == 0:
                time.sleep(1)
            #mesh(x, y, u, title='t=%g' % t[n])
            #plot(x, u[:,1], 'o', '-', title='t=%g' % t[n], ylim=[-2.4, 2.4])
            #plot(x, u_exact(x, t[n], a, Lx), 'o', ylim=[-2.4, 2.4])
            plot(y, u[1,:],'--',y, u_exact(y, t[n], a, Ly),'o')
            #plt.plot(x, u_exact(x, t[n], a, Lx))
            #plt.show()
            filename = 'tmp_%04d.png' % n
            savefig(filename)
            time.sleep(0.1)

        #     Construct problem object
        problem = Problem(I=I, V=V, f=f, q=q, b=b, Lx=Lx, 
                          Ly=Ly, T=T)
        
        #     Construct solver object
        solver = Solver(problem=problem, Nx=Nx, Ny=Ny,
                        dt=dt, user_action=plot_u, 
                        version=self.version)
        #     Solve the PDE
        solver.solve()
        
        
    def test_standing_wave(self):
        """
        Test the solver with a manufactured solution, which
        here is the standing wave
        
        u(x,y,t) = exp(-b*t)*cos(wt)*cos(mx*x*pi/Lx)
        *cos(my*y*pi/Ly).
        
        """
            
        b = 2.3#2.0
        w = 10.7
        c = 1.3
        mx = 1
        my = 1
        Lx = 2.3
        Ly = 2.3
        Nx = 40
        Ny = 40
        T = 0.02#10.0
        dt = 0.002
        print 'dt=', dt
        print 'stability = ', c**2*dt**2/((Lx/Nx)**2) \
            + c**2*dt**2/((Ly/Ny)**2)
        
        def u_exact(x, y, t):
            return exp(-b*t)*cos(w*t)*cos(mx*x*pi/Lx) \
                *cos(my*y*pi/Ly)
        
        def q(x, y):
            return c**2
        
        def I(x, y):
            return cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)
            
        def V(x, y):
            return -b*I(x, y)
        
        def f(x, y, t):
            return (-w**2 + b*w*tan(w*t) + c**2*((mx*pi/Lx)**2 \
                    + (my*pi/Ly)**2))*u_exact(x, y, t)
        
        def print_u(u, x, xv, y, yv, t, n):
            u_e = zeros((len(x), len(y)))
            for j in range (0, len(y)):
                for i in range (0, len(x)):
                    u_e[i,j] = u_exact(x[i], y[j], t[n])
            difference = abs(u[:,:]-u_e[:,:]).max()
            #difference = abs(u[1:len(x)-1,1:len(y)-1]
            #                 -u_e[1:len(x)-1,1:len(y)-1]).max()
            print difference
                        
        def plot_u(u, x, xv, y, yv, t, n):
            #u_e = zeros((len(x), len(y)))
            #for j in range (0, len(y)):
            #    for i in range (0, len(x)):
            #        u_e[i,j] = u_exact(x[i], y[j], t[n])
            surf(xv, yv, u, show=True, zlim=[-2.4, 2.4])
            #surf(xv, yv, u_e, show=True, zlim=[-2.4, 2.4])
            #time.sleep(0.05)
        
        #     Construct problem object
        problem = Problem(I=I, V=V, f=f, q=q, b=b, Lx=Lx, 
                          Ly=Ly, T=T)
        
        #     Construct solver object
        solver = Solver(problem=problem, Nx=Nx, Ny=Ny,
                        dt=dt, user_action=print_u, 
                        version=self.version)
        #     Solve the PDE
        solver.solve()


    def test_convergence_sw(self):
        """
        Test the solver with a manufactured solution, which
        here is the standing wave
        
        u(x,y,t) = exp(-b*t)*cos(wt)*cos(mx*x*pi/Lx)
        *cos(my*y*pi/Ly).
        
        """  
        
        b = 2.3#2.0
        w = 1.7#1.7
        c = 1.3
        mx = 1
        my = 1
        Lx = 2.3
        Ly = 2.3
        Ft = 1.0
        Fx = 30.0
        Fy = 30.0 
        T = 0.01
        error_h = 1000.0
        
        #print 'stability = ', c**2*dt**2/((Lx/Nx)**2) \
        #    + c**2*dt**2/((Ly/Ny)**2)
        
        def q(x, y):
            return c**2
        
        def I(x, y):
            return cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)
            
        def V(x, y):
            return -b*I(x, y)
            
        def f(x, y, t):
            return (-w**2 + b*w*tan(w*t)  \
                         + c**2*((mx*pi/Lx)**2 + (my*pi/Ly)**2)) \
                         *u_exact(x, y, t)
            
        #     Construct problem object
        problem = Problem(I=I, V=V, f=f, q=q, b=b, Lx=Lx, 
                          Ly=Ly, T=T)
        
        def u_exact(x, y, t):
            """
            The exact solution of u(x, y, t).
            """
            return exp(-b*t)*cos(w*t)*cos(mx*x*pi/Lx) \
                *cos(my*y*pi/Ly)
        
        def calc_diff(u, x, y, t):
            """
            Get the largest error between the numerical
            and exact solutions.
            """
            u_e = zeros((len(x), len(y)))
            for j in range (0, len(y)):
                for i in range (0, len(x)):
                    u_e[i,j] = u_exact(x[i], y[j], t)
            error_h = abs(u[:,:]-u_e[:,:]).max()
            #print u[:,:]-u_e[:,:]
            return error_h

        def print_u():
            print 

        errors = []
        h_values = 0.01, 0.001, 0.0001, 0.00001
        #h_values = 0.1, 0.05, 0.01
        
        for h in h_values:
                
            dt = Ft*h
            #     dx = Lx/Nx = Fx*h, dy= Ly/Ny = Fy*h
            Nx = int(round(Lx/(Fx*h)))
            Ny = int(round(Ly/(Fy*h))) 
            #Nx = 20
            #Ny = 20
            
            print 'dt=',dt,',Nx=',Nx,',Ny=',Ny
            #     Construct solver object
            solver = Solver(problem=problem, Nx=Nx, Ny=Ny,
                            dt=dt, user_action=None, 
                            version=self.version)
            #     Solve the PDE
            solver.solve()
            
            error_h = calc_diff(solver.u_last, solver.x, 
                                solver.y, solver.t_last)
            
            print '   error_h = ', error_h
            errors.append(error_h)
            
        nh = len(h_values)
        
        #     Calculate convergence rates
        r = [log(errors[i-1]/errors[i])/ \
                 log(h_values[i-1]/h_values[i]) \
                 for i in range(1, nh, 1)]
        
        for i in range(1, nh):
            print h_values[i-1], r[i-1]


    def test_convergence_sw2(self):
        """
        Test the solver with a manufactured solution, which
        here is the standing wave
        
        u(x,y,t) = exp(-b*t)*cos(wt)*cos(mx*x*pi/Lx)
        *cos(my*y*pi/Ly).
        
        """  
        
        b = 2.3#2.0
        w = 1.7#1.7
        c = 1.3
        mx = 1
        my = 1
        Lx = 2.3
        Ly = 2.3
        Ft = 1.0
        Fx = 30.0
        Fy = 30.0 
        T = 0.01
        error_h = 1000.0
        
        #print 'stability = ', c**2*dt**2/((Lx/Nx)**2) \
        #    + c**2*dt**2/((Ly/Ny)**2)
        
        def q(x, y):
            return exp(-x-y)
        
        def I(x, y):
            return cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)
            
        def V(x, y):
            return -b*I(x, y)
            
        def f(x, y, t):
            cx = mx*pi/Lx
            cy = my*pi/Ly
            return (-w**2 + b*w*tan(w*t)  \
                         + exp(-x-y)*(cx**2 + cy**2 \
                         - cx*tan(cx*x) - cy*tan(cy*y))) \
                         *u_exact(x, y, t)
            
        #     Construct problem object
        problem = Problem(I=I, V=V, f=f, q=q, b=b, Lx=Lx, 
                          Ly=Ly, T=T)
        
        def u_exact(x, y, t):
            """
            The exact solution of u(x, y, t).
            """
            return exp(-b*t)*cos(w*t)*cos(mx*x*pi/Lx) \
                *cos(my*y*pi/Ly)
        
        def calc_diff(u, x, y, t):
            """
            Get the largest error between the numerical
            and exact solutions.
            """
            u_e = zeros((len(x), len(y)))
            for j in range (0, len(y)):
                for i in range (0, len(x)):
                    u_e[i,j] = u_exact(x[i], y[j], t)
            error_h = abs(u[:,:]-u_e[:,:]).max()
            #print u[:,:]-u_e[:,:]
            return error_h

        def print_u():
            print 

        errors = []
        h_values = 0.01, 0.001, 0.0001#, 0.00001
        #h_values = 0.1, 0.05, 0.01
        
        for h in h_values:
                
            dt = Ft*h
            #     dx = Lx/Nx = Fx*h, dy= Ly/Ny = Fy*h
            Nx = int(round(Lx/(Fx*h)))
            Ny = int(round(Ly/(Fy*h))) 
            #Nx = 20
            #Ny = 20
            
            print 'dt=',dt,',Nx=',Nx,',Ny=',Ny
            #     Construct solver object
            solver = Solver(problem=problem, Nx=Nx, Ny=Ny,
                            dt=dt, user_action=None, 
                            version=self.version)
            #     Solve the PDE
            solver.solve()
            
            error_h = calc_diff(solver.u_last, solver.x, 
                                solver.y, solver.t_last)
            
            print '   error_h = ', error_h
            errors.append(error_h)
            
        nh = len(h_values)
        
        #     Calculate convergence rates
        r = [log(errors[i-1]/errors[i])/ \
                 log(h_values[i-1]/h_values[i]) \
                 for i in range(1, nh, 1)]
        
        for i in range(1, nh):
            print h_values[i-1], r[i-1]


          
def main():    
    t = Tests(version='vectorized')
    #t.do_test('1d_gaussian_y') 
    t.do_test('convergence_sw2')
    #t.do_test('standing_wave')
    
if __name__ == '__main__':
    main()
