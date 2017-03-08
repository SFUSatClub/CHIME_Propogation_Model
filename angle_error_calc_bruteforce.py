

from sympy import *
from sympy import lambdify, Matrix
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import cm
import pylab
import math
import csv





   
#constants:

lamda = 700E6 			       # wavelength of transmitted signal 


def relative_symbolic_error():

	''' 
	note: we are calculating the attitude with respect to antenna.
	we are using symmetry in antenna pattern to justify usage of only one angle to describe the attitude
 
	Distance  = distance between cubesat and CHIME
	Phi = angle with respect to antenna's pointing 
	
	dDistance  = fractional error in Distance
	dPhi = fractional error in Phi

	'''

	Distance, Phi = symbols('Distance Phi')
	dDistance, dPhi = symbols('dDistance dPhi')	
	        
	#P = 10 ** (0.2 * exp(-(Phi/0.445) ** 2)) * (lamda/(4 * math.pi * Distance)) # this power function is evaluated at the antenna, propagation power drops by 1/r^2?
	P = 10 ** (0.2 * exp(-(Phi/0.445) ** 2)) * (lamda/(4 * math.pi * Distance ** 2))

	rel_err = sqrt((diff(P, Distance) * dDistance) ** 2 + (diff(P, Phi) * dPhi) ** 2)     # assumes no corrolation between Phi, Distance  
	
        return rel_err
    

def expected_signal_power(Distance_min, Distance_max, Phi_min, Phi_max, sample):

	'''
	Calculating signal power for specified ranges of Distance (distance to CHIME) and Phi
	'''
	fname = raw_input('Please enter the file name:')
	with open('%s.csv' % fname, 'w') as csvfile:
    		fieldnames = ['Power', 'Phi', 'Distance', 'dPower', 'dPhi', 'dDistance']
        	writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        	writer.writeheader()
	   	Distance_list     = np.linspace(Distance_min, Distance_max, sample).tolist()
	   	Phi_list    = np.linspace(Phi_min, Phi_max, sample).tolist()

		for i_Distance in Distance_list:
			for i_Phi in Phi_list:
				P = (0.2 * exp(-(i_Phi/0.445) ** 2)) * (lamda/(4 * math.pi * i_Distance ** 2)) # this power function is evaluated at the antenna
				writer.writerow({'Power': P, 'Phi': i_Phi, 'Distance': i_Distance})



def expected_power_error(Distance_min=700000, Distance_max=800000, Distance_err_min=1000, Distance_err_max=10000, Phi_min=0, Phi_max=1, Phi_err_min=0.01, Phi_err_max=0.4, sample=100, err_sample=1000, power_error_tolerance=0.01, error_required='dPhi'):

	'''
	This calculates everything: (P, dP, dPhi or dDistance)
	error required for either dPhi or dDistance while dP = power_error_tolerance (e.g. dP = 0.01 or 1% error in Power) through BRUTE FODistanceE.
	parameters: Distance_min/max range, Distance_error min/max range, Phi_min/max range, Phi_err_min/max range, sample: spacing between numbers (Phi and Distance)
	err_sample: spacing between numbers (Phi_err and Distance_err)

	Note: "error_required" means you are either looking for dPhi (error in Phi) or dDistance (error in Distance). Pick either "dPhi" or "dDistance"
	'''

	Distance, Phi = symbols('Distance Phi')
	dDistance, dPhi = symbols('dDistance dPhi')
	P = (0.2 * exp(-((Phi)/0.445) ** 2)) * (lamda/(4 * math.pi * Distance ** 2))
	rel_err = sqrt((diff(P, Distance) * dDistance) ** 2 + (diff(P, Phi) * dPhi) ** 2)

	fname = raw_input('Please enter the file name:')
	with open('%s.csv' % fname, 'w') as csvfile:
    		fieldnames = ['Power', 'Phi', 'Distance', 'dPower', 'dPhi', 'dDistance']
        	writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        	writer.writeheader()

		Distance_list     = np.linspace(Distance_min, Distance_max, sample).tolist()
		Phi_list    = np.linspace(Phi_min, Phi_max, sample).tolist()

		Distance_err_list  = np.linspace(Distance_err_min, Distance_err_max, err_sample)
		Phi_err_list = np.linspace(Phi_err_min, Phi_err_max, err_sample)


		for i_Distance in Distance_list:

			for i_Phi in Phi_list:

				eval_P  = P.subs([('Distance', i_Distance), ('Phi', i_Phi)])

				if error_required == 'dPhi':

					for j_Distance_err in Distance_err_list:
						prm_rel_err = rel_err.subs([('Distance', i_Distance), ('Phi', i_Phi), ('dDistance', j_Distance_err)])
						solved_err = solve(prm_rel_err - power_error_tolerance, error_required)
						print P, power_error_tolerance, i_Distance, i_Phi, solved_err, j_Distance_err 
						writer.writerow({'Distance': i_Distance, 'Phi': i_Phi, 'Power' : eval_P, 'dPower': power_error_tolerance, 'dPhi': solved_err, 'dDistance': j_Distance_err})	

				elif error_required == 'dDistance':

					for j_Phi_err in Phi_err_list:
						prm_rel_err = rel_err.subs([('Distance', i_Distance), ('Phi', i_Phi), ('dPhi', j_Phi_err)])
						solved_err = solve(prm_rel_err - power_error_tolerance, error_required)
						print P, power_error_tolerance, i_Distance, i_Phi, j_Phi_err, solved_err
						writer.writerow({'Distance': i_Distance, 'Phi': i_Phi, 'Power' : P, 'dPower': power_error_tolerance, 'dPhi': j_Phi_err, 'dDistance': solved_err})	
	
				else:
					print "invalid input: please select either 'dPhi' or 'dDistance'"	



	
			
			

def plot(surf_fun, d_Phi, d_Distance, sample):

    fig = plt.figure(figsize=(16,16))
    ax  = fig.gca(projection='3d')
    
    dDistance, dPhi = np.linspace(d_Distance[0], d_Distance[1], sample), np.linspace(d_Phi[0], d_Phi[1], sample)
    dDistance, dPhi = np.meshgrid(dDistance, dPhi)

    f = lambdify(('dDistance', 'dPhi'), surf_fun, 'numpy') 
    surf = f(dDistance, dPhi)
    
    ax.plot_surface(dDistance, dPhi, surf, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.set_zlim(0, 0.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    ax.set_title('error surface for P')
    ax.set_xlabel('dDistance')
    ax.set_ylabel('dPhi')
    ax.set_zlabel('dP')



#rel_err = relative_symbolic_error()
#prm_rel_err = rel_err.subs([('Phi', 0.75), ('Distance', 800000)])
#surf_fun = solve(prm_rel_err - 0.01, 'dPhi')[1]

#surf_fun = solve(rel_err - 0.01, 'dPhi')[1]
#plot(surf_fun, d_Distance=(9000, 10000), d_Phi=(0, 0.17), sample=100)
#pylab.show()



exp_err_dPhi = expected_power_error()

