
import sympy as sp
from sympy import *
from sympy import lambdify, Matrix
from sympy.plotting import plot3d
#from mpmath import *
#from sympy import Symbol, solveset, S, erf, log, sqrt



import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import cm
import pylab
import math
import csv

#from pyPdf import PdfFileWriter, PdfFileReader





   
#constants:

lamda = 0.3747 			       # wavelength of transmitted signal 


def relative_symbolic_error():

	''' 
	note: we are calculating the attitude with respect to antenna.
	we are using symmetry in antenna pattern to justify usage of only one angle to describe the attitude
 
	Distance  = distance between cubesat and CHIME
	Theta = angle with respect to antenna (zenith angle)
	
	dDistance  = fractional error in Distance
	dTheta = fractional error in Theta

	'''

	Distance, Theta = symbols('Distance Theta')
	dDistance, dTheta = symbols('dDistance dTheta')	
	        
	#P = 10 ** (0.2 * exp(-(Theta/0.445) ** 2)) * (lamda/(4 * math.pi * Distance)) # this power function is evaluated at the antenna, propagation power drops by 1/r^2?
	#P = (10 ** (0.2 * exp(-((Theta)/0.445) ** 2)))* (lamda/(4 * math.pi * Distance)) ** 2
	P = (1.644 * 0.71) * (15.0/(math.pi)) * ((1.0/Distance) ** 2) * ((sp.cos((math.pi/2) * sp.cos(Theta)))/(sp.sin(Theta))) ** 4 
	rel_err = sqrt((diff(P, Distance) * dDistance) ** 2 + (diff(P, Theta) * dTheta) ** 2) / P  # assumes no corrolation between Theta, Distance  
	print rel_err
        return rel_err
    



def expected_error_values(Distance_min=500000, Distance_max=800000, Distance_err_min=2000, Distance_err_max=2000, Theta_min=0.0, Theta_max=0.2, Theta_err_min=0.01, Theta_err_max=0.4, sample=10, err_sample=100, power_error_tolerance=0.01, error_required='dTheta'):

	'''
	This calculates everything: (P, dP, dTheta or dDistance for given ranges of Theta and Distance)
	The calculation solves for dTheta or dDistance while dP = power_error_tolerance (e.g. dP = 0.01 or 1% fractional error in Power) through BRUTE FORCE.
	parameters: Distance_min/max range, Distance_error min/max range, Theta_min/max range, Theta_err_min/max range, sample: spacing between numbers (Theta and Distance)
	err_sample: spacing between numbers (Theta_err and Distance_err)

	Note: "error_required" means you are either looking for dTheta (error in Theta) or dDistance (error in Distance). Pick either "dTheta" or "dDistance"
	'''


	Distance, Theta = symbols('Distance Theta')
	dDistance, dTheta = symbols('dDistance dTheta')
	P = (10 ** (0.2 * exp(-((Theta)/0.445) ** 2)))* (lamda/(4 * math.pi * Distance)) ** 2 
	rel_err = sqrt((diff(P, Distance) * dDistance) ** 2 + (diff(P, Theta) * dTheta) ** 2) / P

	fname = raw_input('Please enter the file name you want to dump the data into with no extension:')
	with open('%s.csv' % fname, 'w') as csvfile:
    		fieldnames = ['Power', 'Theta', 'Distance', 'dPower', 'dTheta', 'dDistance']
        	writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        	writer.writeheader()

		Distance_list     = np.linspace(Distance_min, Distance_max, sample).tolist()
		Theta_list    = np.linspace(Theta_min, Theta_max, sample).tolist()

		Distance_err_list  = np.linspace(Distance_err_min, Distance_err_max, err_sample)
		Theta_err_list = np.linspace(Theta_err_min, Theta_err_max, err_sample)


		for i_Distance in Distance_list:

			for i_Theta in Theta_list:

				eval_P  = P.subs([('Distance', i_Distance), ('Theta', i_Theta)])
				

				if error_required == 'dTheta':
					
					#print 'Solving for %s when dP = %s over the domains: Distance_min = %s, Distance_max = %s, Theta_min = %s, Theta_max = %s, 	   						Distance_error_min = %s, Distance_error_max = %s .' % (error_required, str(power_error_tolerance), 						str(Distance_min), str(Distance_max), str(Theta_min), str(Theta_max), str(Distance_err_min), str(Distance_err_max))
					for j_Distance_err in Distance_err_list:
						prm_rel_err = rel_err.subs([('Distance', i_Distance), ('Theta', i_Theta), ('dDistance', j_Distance_err)])
						solved_err = solve(prm_rel_err - power_error_tolerance, error_required)
						print power_error_tolerance, i_Distance, i_Theta, solved_err, j_Distance_err 
						writer.writerow({'Distance': i_Distance, 'Theta': i_Theta, 'dPower': power_error_tolerance, 									'dTheta': solved_err, 'dDistance': j_Distance_err})	

				elif error_required == 'dDistance':

					#print 'Solving for %s when dP = %s over the domains: Distance_min = %s, Distance_max = %s, Theta_min = %s, Theta_max = %s, 	   						Theta_error_min = %s, Theta_error_max = %s .' % (error_required, str(power_error_tolerance), 						str(Distance_min), str(Distance_max), str(Theta_min), str(Theta_max), str(Theta_err_min), str(Theta_err_max))

					for j_Theta_err in Theta_err_list:
						prm_rel_err = rel_err.subs([('Distance', i_Distance), ('Theta', i_Theta), ('dTheta', j_Theta_err)])
						solved_err = solve(prm_rel_err - power_error_tolerance, error_required)
						print power_error_tolerance, i_Distance, i_Theta, j_Theta_err, solved_err
						writer.writerow({'Distance': i_Distance, 'Theta': i_Theta, 'dPower': power_error_tolerance, 									'dTheta': j_Theta_err,'dDistance': solved_err})	
	
				else:
					print "invalid input: please select either 'dTheta' or 'dDistance'"	



	


def expected_error_plot(fn='dP/P', dTheta_range=(0.015,0.03), dDist_range=(1000,2000), Theta_range=(1.44, 1.57), dist=800000, power_error_tolerance=0.01, sample=100):

	'''
	Depending on user input for fn which can be: dP/P or dTheta or dDistance
	This plots relative errors dP/P as a funciton (dTheta, dDistance, Theta, Distance) or
	plots dTheta for dP/P = power_error_tolerance as function of dDistance or 
	dDistance for dP/P = power_error_tolerance as function dTheta
	'''

	
	Theta = np.linspace(Theta_range[0], Theta_range[1], 4).tolist()



	power_error_tolerance = 0	
	dDistance , dTheta = symbols('dDistance dTheta')
	rel_err = relative_symbolic_error()


	for i_Theta in Theta:

		Theta_degree = i_Theta * 180/math.pi
		Theta_degree = str(Theta_degree)[:4]
	
		if fn == 'dP/P':


    			dDistance, dTheta = np.linspace(dDist_range[0], dDist_range[1], sample), np.linspace(dTheta_range[0], dTheta_range[1], sample)
			dDistance, dTheta = np.meshgrid(dDistance, dTheta)


			prm_rel_err = rel_err.subs([('Distance', dist), ('Theta', i_Theta)])
    			prm_rel_err_lambdified = lambdify(('dDistance', 'dTheta'), prm_rel_err, 'numpy')
			surf = prm_rel_err_lambdified(dDistance, dTheta)

			plt.subplot(2, 2, Theta.index(i_Theta)+1)
			plt.pcolormesh(dTheta, dDistance, surf, cmap='RdBu', vmin=0.005, vmax=0.011)
			plt.axis([dTheta.min(), dTheta.max(), dDistance.min(), dDistance.max()])
    			plt.ylabel('$\delta$ Distance (m)', fontsize=18)
    			plt.xlabel(r'$\delta$ $\theta$ (rads)', fontsize=18)
			plt.title(r'$\frac{dA}{A}$ at |$\theta$| = %s degrees with respect to the Nadir.' % Theta_degree, fontsize=20, y=1.08)
			plt.tick_params(axis='both', which='major', labelsize=16)
			plt.tick_params(axis='both', which='minor', labelsize=16)
			plt.colorbar()
		elif fn == 'dTheta':
			dDistance = np.linspace(dDist_range[0], dDist_range[1], sample)
			prm_rel_err = rel_err.subs([('Distance', dist), ('Theta', i_Theta)])
			solved_err = solve(prm_rel_err - power_error_tolerance, fn)
			solved_err_lambdified = lambdify(('dDistance'), solved_err, 'numpy')
			f = solved_err_lambdified(dDistance)		
			fn_plot = plt.plot(dDistance, f)
			plt.subplot(2, 2, Theta.index(i_Theta)+1)
			plt.setp(lines, color='r', linewidth=2.0)
			plt.xlabel('$\delta$ Distance (m)', fontsize=16)
			plt.ylabel(r'$\delta$ $\theta$ (rads)', fontsize=16)
			plt.title(r'$delta$ $\theta$ at |$\theta$| = %s as a function $\delta$ Distance while $\frac{\delta A}{A}$ = %s' % i_Theta, power_error_tolerance, y=1.08)

		elif fn == 'dDistance':
			dTheta = np.linspace(dTheta_range[0], dTheta_range[1], sample)
			prm_rel_err = rel_err.subs([('Distance', dist), ('Theta', i_Theta)])
			solved_err = solve(prm_rel_err - power_error_tolerance, fn)
			solved_err_lambdified = lambdify(('dTheta'), solved_err, 'numpy')
			f = solved_err_lambdified(dTheta)
			fn_plot = plt.plot(dTheta, f)
			plt.subplot(2, 2, Theta.index(i_Theta)+1)
			plt.setp(lines, color='r', linewidth=2.0)
			plt.ylabel('$\delta$ Distance (m)', fontsize=16)
			plt.xlabel(r'$\delta$ $\theta$ (rads)', fontsize=16)
			plt.title(r'$delta$ $\theta$ at |$\theta$| = %s as a function $\delta$ Distance while $\frac{\delta A}{A}$ = %s' % i_Theta,power_error_tolerance, y=1.08)
	plt.subplots_adjust(wspace=0.5, hspace=0.5)	
	pylab.show()
	


def append_pdf(input,output):
    [output.addPage(input.getPage(page_num)) for page_num in range(input.numPages)]

# Creating an object where pdf pages are appended to
#output = PdfFileWriter()

# Appending two pdf-pages from two different files
#append_pdf(PdfFileReader(open("SamplePage1.pdf","rb")),output)
#append_pdf(PdfFileReader(open("SamplePage2.pdf","rb")),output)

# Writing all the collected pages to a file
#output.write(open("CombinedPages.pdf","wb"))



exp_err_plot = expected_error_plot()



