

from sympy import *
from sympy import lambdify, Matrix
from sympy.plotting import plot3d

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
	Phi = angle with respect to antenna's pointing 
	
	dDistance  = fractional error in Distance
	dPhi = fractional error in Phi

	'''

	Distance, Phi = symbols('Distance Phi')
	dDistance, dPhi = symbols('dDistance dPhi')	
	        
	#P = 10 ** (0.2 * exp(-(Phi/0.445) ** 2)) * (lamda/(4 * math.pi * Distance)) # this power function is evaluated at the antenna, propagation power drops by 1/r^2?
	P = (10 ** (0.2 * exp(-((Phi)/0.445) ** 2)))* (lamda/(4 * math.pi * Distance)) ** 2 
	rel_err = sqrt((diff(P, Distance) * dDistance) ** 2 + (diff(P, Phi) * dPhi) ** 2) / P  # assumes no corrolation between Phi, Distance  
	
        return rel_err
    



def expected_error_values(Distance_min=500000, Distance_max=800000, Distance_err_min=2000, Distance_err_max=2000, Phi_min=0.0, Phi_max=0.2, Phi_err_min=0.01, Phi_err_max=0.4, sample=10, err_sample=100, power_error_tolerance=0.01, error_required='dPhi'):

	'''
	This calculates everything: (P, dP, dPhi or dDistance for given ranges of Phi and Distance)
	The calculation solves for dPhi or dDistance while dP = power_error_tolerance (e.g. dP = 0.01 or 1% fractional error in Power) through BRUTE FORCE.
	parameters: Distance_min/max range, Distance_error min/max range, Phi_min/max range, Phi_err_min/max range, sample: spacing between numbers (Phi and Distance)
	err_sample: spacing between numbers (Phi_err and Distance_err)

	Note: "error_required" means you are either looking for dPhi (error in Phi) or dDistance (error in Distance). Pick either "dPhi" or "dDistance"
	'''


	Distance, Phi = symbols('Distance Phi')
	dDistance, dPhi = symbols('dDistance dPhi')
	P = (10 ** (0.2 * exp(-((Phi)/0.445) ** 2)))* (lamda/(4 * math.pi * Distance)) ** 2 
	rel_err = sqrt((diff(P, Distance) * dDistance) ** 2 + (diff(P, Phi) * dPhi) ** 2) / P

	fname = raw_input('Please enter the file name you want to dump the data into with no extension:')
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
					
					#print 'Solving for %s when dP = %s over the domains: Distance_min = %s, Distance_max = %s, Phi_min = %s, Phi_max = %s, 	   						Distance_error_min = %s, Distance_error_max = %s .' % (error_required, str(power_error_tolerance), 						str(Distance_min), str(Distance_max), str(Phi_min), str(Phi_max), str(Distance_err_min), str(Distance_err_max))
					for j_Distance_err in Distance_err_list:
						prm_rel_err = rel_err.subs([('Distance', i_Distance), ('Phi', i_Phi), ('dDistance', j_Distance_err)])
						solved_err = solve(prm_rel_err - power_error_tolerance, error_required)
						print power_error_tolerance, i_Distance, i_Phi, solved_err, j_Distance_err 
						writer.writerow({'Distance': i_Distance, 'Phi': i_Phi, 'dPower': power_error_tolerance, 									'dPhi': solved_err, 'dDistance': j_Distance_err})	

				elif error_required == 'dDistance':

					#print 'Solving for %s when dP = %s over the domains: Distance_min = %s, Distance_max = %s, Phi_min = %s, Phi_max = %s, 	   						Phi_error_min = %s, Phi_error_max = %s .' % (error_required, str(power_error_tolerance), 						str(Distance_min), str(Distance_max), str(Phi_min), str(Phi_max), str(Phi_err_min), str(Phi_err_max))

					for j_Phi_err in Phi_err_list:
						prm_rel_err = rel_err.subs([('Distance', i_Distance), ('Phi', i_Phi), ('dPhi', j_Phi_err)])
						solved_err = solve(prm_rel_err - power_error_tolerance, error_required)
						print power_error_tolerance, i_Distance, i_Phi, j_Phi_err, solved_err
						writer.writerow({'Distance': i_Distance, 'Phi': i_Phi, 'dPower': power_error_tolerance, 									'dPhi': j_Phi_err,'dDistance': solved_err})	
	
				else:
					print "invalid input: please select either 'dPhi' or 'dDistance'"	



	


def expected_error_plot(fn='dP/P', dPhi_range=(0.015,0.03), dDist_range=(1000,2000), phi_range=(0,0.088), dist=800000, power_error_tolerance=0.01, sample=100):

	'''
	Depending on user input for fn which can be: dP/P or dPhi or dDistance
	This plots relative errors dP/P as a funciton (dPhi, dDistance, Phi, Distance) or
	plots dPhi for dP/P = power_error_tolerance as function of dDistance or 
	dDistance for dP/P = power_error_tolerance as function dPhi
	'''

	
	phi = np.linspace(phi_range[0], phi_range[1], 4).tolist()



	power_error_tolerance = 0	
	dDistance , dPhi = symbols('dDistance dPhi')
	rel_err = relative_symbolic_error()


	for i_phi in phi:

		phi_degree = i_phi * 180/math.pi
		phi_degree = str(phi_degree)[:4]
	
		if fn == 'dP/P':


    			dDistance, dPhi = np.linspace(dDist_range[0], dDist_range[1], sample), np.linspace(dPhi_range[0], dPhi_range[1], sample)
			dDistance, dPhi = np.meshgrid(dDistance, dPhi)


			prm_rel_err = rel_err.subs([('Distance', dist), ('Phi', i_phi)])
    			prm_rel_err_lambdified = lambdify(('dDistance', 'dPhi'), prm_rel_err, 'numpy')
			surf = prm_rel_err_lambdified(dDistance, dPhi)

			plt.subplot(2, 2, phi.index(i_phi)+1)
			plt.pcolormesh(dPhi, dDistance, surf, cmap='RdBu', vmin=0.005, vmax=0.01)
			plt.axis([dPhi.min(), dPhi.max(), dDistance.min(), dDistance.max()])
    			plt.ylabel('$\delta$ Distance (m)', fontsize=18)
    			plt.xlabel('$\delta$ $\phi$ (rads)', fontsize=18)
			plt.title(r'$\frac{dP}{P}$ at |$\phi$| = %s degrees with respect to the Nadir.' % phi_degree, fontsize=20, y=1.08)
			plt.tick_params(axis='both', which='major', labelsize=16)
			plt.tick_params(axis='both', which='minor', labelsize=16)
			plt.colorbar()
		elif fn == 'dPhi':
			dDistance = np.linspace(dDist_range[0], dDist_range[1], sample)
			prm_rel_err = rel_err.subs([('Distance', dist), ('Phi', i_phi)])
			solved_err = solve(prm_rel_err - power_error_tolerance, fn)
			solved_err_lambdified = lambdify(('dDistance'), solved_err, 'numpy')
			f = solved_err_lambdified(dDistance)		
			fn_plot = plt.plot(dDistance, f)
			plt.subplot(2, 2, phi.index(i_phi)+1)
			plt.setp(lines, color='r', linewidth=2.0)
			plt.xlabel('$\delta$ Distance (m)', fontsize=16)
			plt.ylabel('$\delta$ $\phi$ (rads)', fontsize=16)
			plt.title(r'$delta$ $\phi$ at |$\phi$| = %s as a function $\delta$ Distance while $\frac{\delta P}{P}$ = %s' % i_phi, power_error_tolerance, y=1.08)

		elif fn == 'dDistance':
			dPhi = np.linspace(dPhi_range[0], dPhi_range[1], sample)
			prm_rel_err = rel_err.subs([('Distance', dist), ('Phi', i_phi)])
			solved_err = solve(prm_rel_err - power_error_tolerance, fn)
			solved_err_lambdified = lambdify(('dPhi'), solved_err, 'numpy')
			f = solved_err_lambdified(dPhi)
			fn_plot = plt.plot(dPhi, f)
			plt.subplot(2, 2, phi.index(i_phi)+1)
			plt.setp(lines, color='r', linewidth=2.0)
			plt.ylabel('$\delta$ Distance (m)', fontsize=16)
			plt.xlabel('$\delta$ $\phi$ (rads)', fontsize=16)
			plt.title(r'$delta$ $\phi$ at |$\phi$| = %s as a function $\delta$ Distance while $\frac{\delta P}{P}$ = %s' % i_phi,power_error_tolerance, y=1.08)
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



