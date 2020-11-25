import flopy
import os
import numpy as np

model_dir = os.path.join(os.getcwd(),"data")
model_name = 'hooghoudt'

# tdis
nper = 1
tdis_rc = []
for i in range(nper):
    tdis_rc.append((1., 1, 1))

# solver data
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-9, 1e-3, 0.97

# model spatial discretization
nlay, nrow, ncol = 1, 258, 258

# dimenions
Lx = 80.
Ly = 80.
Lz = 1.

# cell spacing
delr = Lx/(ncol-1)
delc = Ly/(nrow-1)
area = delr * delc

# top/bot of the aquifer
top = 1.0
tops = np.linspace(top, top-Lz, nlay+1)

# hydraulic conductivity
hk = 10.

# boundary heads
h1 = 1.
h2 = 1.

# initial
hi = 0.

# build chd stress period data
left_bnd = [[(ilay,irow,0), h1] for ilay in range(nlay) for irow in range(nrow)]
right_bnd = [[(ilay,irow,ncol-1), h2] for ilay in range(nlay) for irow in range(nrow)]
bnd_data = left_bnd + right_bnd
chd_spd = {0: bnd_data}

# average recharge rate
avg_rch = 5.0e-04

# build recharge spd
rch_spd = {}
for n in range(nper):
    rch_spd[n] = avg_rch


sim = flopy.mf6.MFSimulation(sim_name="simple", version='mf6',
                                 exe_name='mf6',
                                 sim_ws=model_dir, memory_print_option='all')
# create tdis package

tdis = flopy.mf6.ModflowTdis(sim, time_units='DAYS',
                             nper=nper, perioddata=tdis_rc)

# create iterative model solution and register the gwf model with it
ims = flopy.mf6.ModflowIms(sim,
                           print_option='SUMMARY',
                           outer_hclose=hclose,
                           outer_maximum=nouter,
                           under_relaxation='DBD',
                           inner_maximum=ninner,
                           inner_hclose=hclose, rcloserecord=rclose,
                           linear_acceleration='BICGSTAB',
                           relaxation_factor=relax)

# create gwf model
gwf = flopy.mf6.ModflowGwf(sim, modelname=model_name, save_flows=True)

dis = flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=nrow, ncol=ncol,
                              delr=delr, delc=delc,
                              top=tops[0], botm=tops[1:])

# initial conditions
ic = flopy.mf6.ModflowGwfic(gwf, strt=hi)

# node property flow
npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=True,
                              icelltype=0,
                              k=hk)

# chd file
chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

# recharge file
rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=rch_spd)

# output control
oc = flopy.mf6.ModflowGwfoc(gwf,
                            head_filerecord='{}.hds'.format(model_name),
                            headprintrecord=[
                                ('COLUMNS', 10, 'WIDTH', 15,
                                 'DIGITS', 6, 'GENERAL')],
                            saverecord=[('HEAD', 'ALL')],
                            printrecord=[('BUDGET', 'ALL')])

sim.write_simulation()