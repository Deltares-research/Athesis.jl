# First order explicit upwind scheme for the non-linear Shallow-Water Equations
import Pkg
ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"
Pkg.add("Plots")
using Plots

function thomas_alg(a,b,c,d)
    #
    # Thomas algorithm
    # Name of the function: thomas_alg
    # Name of the file:     thomas_alg.m
    # Arguments:
    #   a = vector of lower diagonal elements
    #   b = vector of diagonal elements
    #   c = vector of upper diagonal elements
    #   d = known right hand side vector
    # Results / output:
    #   x = solution of linear system
    #

    c[1]=c[1]/b[1];
    d[1]=d[1]/b[1];
    #
    N = length(d);      # number of elements in the vector
    #
    x = zeros(Float64,N);
    #
    # Part I: Forward elimination
    for i = 2:N
        K    = 1/(b[i]-c[i-1]*a[i]);
        c[i] = c[i]*K;
        d[i] = (d[i]-a[i]*d[i-1])*K;
    end
    #
    # Part II: Backward elimination
    x[N] = d[N];
    #
    for i = N-1:-1:1
        x[i] = d[i]-c[i]*x[i+1];
    end
    return x
end

function IntTrajectory(xstart, x, dx, dt, q, BCL, BCR)
    # Name of the fuction is IntTrajectory
    # Name of the file (mandatory) is IntTrajectory.m
    # (Can also be a vector function for all points in the domain together, using as vectors xstart and q)
    # Input / arguments:
    #   xstart = starting point of integration trajectory (scalar)
    #   x      = vector of spatial coordinates
    #   dx     = mesh spacing (for the equidistant grid)
    #   dt     = final time of the ODE integrator = physical time
    #   q      = vector of state variables
    #   BCL    = boundary condition on the left
    #   BCR    = boundary condition on the right
    #
    CFL  = 0.1;          # local Courant-Friedrichs-Levy number
    xL   = x[1];         # location of left boundary
    xR   = x[end];       # location of the right boundary
    tau  = 0;            # initialize local time tau
    x0   = xstart;       # starting point of the characteristic
    NMAX = 10000;        # max. number of integration steps
    #
    q0 = 0.0
    #
    for n = 1:NMAX
        if (x0 <= xL)
            q0 = BCL;
            break
        end
        if (x0 >= xR)
            q0 = BCR;
            break
        end
        j    = ceil(Int64,(x0-xL)/dx);                    # current nearest index j associated to point x0
        p    = (x0-x[j])/dx;                          # remainder (real number 0 <= p < 1)
        q0   = (1-p)*q[j] + p*q[j+1];               # linear interpolation of q at x0
        qmax = max( abs(q[j]),abs(q[j+1]) );        # max. local velocity
        dtau = CFL*dx/qmax;                         # local time step
        if (tau + dtau > dt)
           dtau = dt - tau;                         # match final time dt EXACTLY
        end
        if (tau >= dt)
            break
        end
        #x1   = x0 - dtau*q0;                                     # update x0 with first order explicit Euler method
        x1   = x0 - dtau*q0 + 0.5*dtau^2*(q[j+1]-q[j])/dx*q0;    # update x0 with second order explicit Euler method
        x0   = x1;
        tau  = tau + dtau;                          # update the local time
    end
    return q0
end

function timeloop(NMAX, times, time0, dx, grid0, grid, state, bc, physics, checking)
    # Loop in time

    # Unpack the time constants
    dt0 = time0["dt0"]
    tend = time0["tend"]

    # Unpack the grid constants
    imax = grid0["imax"]
    IMAX = grid0["IMAX"]

    # Unpack the grid arrays
    x = grid["x"]
    xb = grid["xb"]

    # Unpack the state
    h = state["h"]
    u = state["u"]
    hu = state["hu"]
    hi = state["hi"]
    eta = state["eta"]
    etaj = state["etaj"]

    # Unpack the boundary conditions
    BCR = bc["BCR"]
    BCL = bc["BCL"]

    # Unpack the physics
    g = physics["grav"]

    # Unpack checking variables
    masstot = checking["masstot"]
    momtot  = checking["momtot"]

    b = Array{Float64,1}(undef, imax)
    d = Array{Float64,1}(undef, imax)
    ld = Array{Float64,1}(undef, imax)
    ud = Array{Float64,1}(undef, imax)

    eta1 = Array{Float64,1}(undef, imax)
    u1 = Array{Float64,1}(undef, IMAX)
    Fu = Array{Float64,1}(undef, IMAX)

    for n = 1:NMAX
    #    amax = max(abs(q))          # Maximum wave propagation speed
    #    dt = C*dx/amax              # Time step
        dt = min(dt0, tend-times[n])
        hu[1] = max(0, h[1]+eta[1])
        #H(1) = 0
        for i = 2:IMAX-1
            hu[i] = max(0, h[i]+eta[i-1],h[i]+eta[i])
        end
        hu[IMAX] = max(0, h[imax]+eta[imax])
        #H(IMAX) = 0
        #Htot(n) = sum(H)
        if (times[n] >= tend)
        #    Hdiff = Htot(n) - Htot(1)
            break
        end

        for i = 1:IMAX
            if (hu[i] <= 0.0)
                u[i] = 0.0                                  # empty (dry)interfaces do not have a velocity
            end
        end

        # Loop in space over all velocity points
        for i = 1:IMAX
            Fu[i] = IntTrajectory(x[i], x, dx, dt, u, BCL, BCR)
        end

        for i = 1:imax
            d[i]  = eta[i] - dt/dx*( hu[i+1]*Fu[i+1] - hu[i]*Fu[i] )
            ld[i] = -g*dt^2/dx^2*hu[i]
            ud[i] = -g*dt^2/dx^2*hu[i+1]
        end
        ld[1]    = 0
        ud[imax] = 0
        b = ones(Float64, imax) - ld - ud
        eta1 = thomas_alg(ld', b', ud', d')'
        u1[1] = BCL
        for i = 2:IMAX-1
            if (hu[i] > 0.01)
                u1[i]   = Fu[i] - g*dt/dx*(eta1[i]-eta1[i-1])
            else
                u1[i] = 0
            end
        end
        u1[IMAX] = BCR
        #
        times[n+1] = times[n] + dt
        for i = 1:imax
            eta[i] = max(eta1[i], -hi[i])
        end   # overwrite old solution with new one
        u  = u1
        #masstot[n+1]   = sum(max(eta1, -hi) + hi)
        masstot[n+1]   = sum(eta + hi)
        etaj[1] = eta[1]
        for i = 2:IMAX-1
            if (i > 1 && i < IMAX)
                etaj[i] = 0.5*(eta[i-1]+eta[i])
            end
        end
        etaj[IMAX]  = eta[imax]
        for i = 1:IMAX
            hu[i] = max(0.0,h[i] + etaj[i])
        end
        momtot[n+1] = sum(hu.*u)
        # Plot
        #subplot(4,1,1)
        plot(xb,[eta, -hi])
        filnam = string("wl\_t=",string(times[n+1]),".png")
        png(filnam)
        #title(sprintf('Current time t = #f',times[n]))
        #title(sprintf('Water level'))
        #axis([xL,xR,-2,2])
        #title(sprintf('Current time t = #f',times[n]))
        #subplot(4,1,2)
        plot(x,u)
        filnam = string("u\_t=",string(times[n+1]),".png")
        png(filnam)
        #axis([xL,xR,-1,1])
        #subplot(4,1,3)
        plot(times,masstot)
        filnam = string("mass\_t=",string(times[n+1]),".png")
        png(filnam)

        #axis([0,tend,40,50])
        #subplot(4,1,4)
        plot(times,momtot)
        filnam = string("momentum\_t=",string(times[n+1]),".png")
        png(filnam)
        #axis([0,tend,-50,50])
        show = true
        #sleep(5)
        #gui()
    end

end

function swe1d_nonlin()

    # Definition of the computational domain
    xL = 0.0                          # left border of spatial domain
    xR = 100.0                          # right border of spatial domain
    t0 = 0.0                          # initial time
    tend = 10.0                      # final time
    # Definition of the advective speed
    a = 1                           # advective speed
      # Discretization parameters
    IMAX = 100                      # number of velocity points in x-direction
    imax = IMAX-1                    # number of pressure points in x-direction
    dx = (xR-xL)/IMAX               # mesh size / mesh spacing
    dt0 = 0.1                       # time step
    #C = 0.9                         # Courant number
    # Definition of the gravity constant
    g = 9.81
    # Define the meshes
    #x  = range(xL, xR, length=IMAX)                # Distribute IMAX equidistant points between xL and xR for the velocity points
    x = Array{Float64,1}(undef, IMAX)
    #x = [xL:dx:xR]
    #xb = range(xL+dx/2, xR-dx/2, length=imax)      # Distribute IMAX equidistant points between xL and xR for the pressure points
    xb = Array{Float64,1}(undef, imax)
    eta = Array{Float64,1}(undef, imax)
    for i = 1:imax
        x[i]  =  xL + dx*(i-1)
        xb[i] =  xL + dx*(i-1/2)
        if x[i] < (xR-xL)/2.0
            eta[i] = 5.0
        else
            eta[i] = 1.0
        end
    end
    # Define the boundary conditions
    # For the velocities
    BCL = 0.0
    BCR = 0.0
    # Define the inital condition
    # First the free surface
    # Now initialize the velocities and the topography (at velocity points)
    u = Array{Float64,1}(undef, IMAX)
    h = Array{Float64,1}(undef, IMAX)
    etaj = Array{Float64,1}(undef, IMAX)
    hu = Array{Float64,1}(undef, IMAX)
    for i = 1:IMAX
        u[i] = 0                # Initial velocity
    #    h(i) = 0                # Bottom profile
        h[i] = 0.0
        if (i > 1 && i < IMAX)
            etaj[i] = 0.5*(eta[i-1]+eta[i])
            hu[i]   = max(0,h[i] + etaj[i])
        end
    end
    hu[1]    = max(0,h[1] + etaj[2])
    hu[IMAX] = max(0,h[IMAX] + etaj[imax])
    hi = Array{Float64,1}(undef, imax)
    for i = 1:imax
        hi[i]      = 0.5*(h[i]+h[i+1])
    end
    eta = max(eta,-hi)

    NMAX = 10000
    masstot = Array{Float64,1}(undef, NMAX)
    momtot  = Array{Float64,1}(undef, NMAX)
    masstot[1] = sum(eta + hi)
    momtot[1] = sum(hu.*u)
    #subplot(2,1,1)
    #plot(xb,eta,'o-')
    #subplot(2,1,2)
    #plot(x,u,'o-')
    # Initialize the current time with the initial time t0
    times = Array{Float64,1}(undef, NMAX)
    times[1] = t0
    i = 1
    # Define boundary conditions
    # Boundary condition on the left: q(xL,t) = q(xL,0)

    time0 = Dict{String,Float64}("tend"=>tend,"dt0"=>dt0)
    grid0 = Dict{String,Int64}("imax"=>imax,"IMAX"=>IMAX)
    grid  = Dict{String,Array{Float64,1}}("x"=>x,"xb"=>xb)
    state = Dict{String,Array{Float64,1}}("h"=>h,"u"=>u, "hi"=>hi, "hu"=>hu, "eta"=>eta, "etaj"=>etaj)
    bc    = Dict{String,Float64}("BCL"=>BCL,"BCR"=>BCR)
    physics = Dict{String,Float64}("grav"=>g)
    checking = Dict{String,Array{Float64,1}}("masstot"=>masstot,"momtot"=>momtot)

    # Loop in time
    timeloop(NMAX, times, time0, dx, grid0, grid, state, bc, physics, checking)
end

swe1d_nonlin()
