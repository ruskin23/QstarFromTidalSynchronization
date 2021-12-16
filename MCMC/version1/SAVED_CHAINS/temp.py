import matplotlib.pyplot as plt
import numpy

s=['13','86','96']

porb=[]
pspin=[]
mass=[]
eccentricity=[]
r=[]

porb_s=[]
pspin_s=[]
mass_s=[]
eccentricity_s=[]
r_s=[]

with open('SpinlogQCatalog_el0.4.txt','r') as f:
    for lines in f:
        x=lines.split()
        if x[0] in s:
            porb=float(x[6])
            pspin=float(x[12])
            omega_orb=2*numpy.pi/porb
            omega_spin=2*numpy.pi/pspin
            tidal_frequency=2*abs(omega_orb-omega_spin)
            print(f'System = {x[0]} I={tidal_frequency/(2*omega_spin)}')

# with open('SpinlogQCatalog_el0.4.txt','r') as f:
#     next(f)
#     for lines in f:
#         x=lines.split()
#         if x[0] in s:
#             porb_s.append(float(x[6]))
#             pspin_s.append(float(x[12]))
#             eccentricity_s.append(float(x[8]))
#             mass_s.append( float(x[14])*float(x[15]))
#             r_s.append(float(x[12])/float(x[6]))
#         else:
#             porb.append(float(x[6]))
#             pspin.append(float(x[12]))
#             eccentricity.append(float(x[8]))
#             mass.append(float(x[14])*float(x[15]))
#             r.append(float(x[12])/float(x[6]))

# plt.scatter(porb_s,eccentricity_s,c='r')
# plt.scatter(porb,eccentricity,c='g')
# # plt.axhline(y=1,linestyle='--')
# plt.show()









