#!/usr/bin/env python3

import numpy as np
from scipy.io import netcdf
import sys, os

def main(file):
	print()
	print("Usage: "+file+" quasisymmetry_out.*.nc")

	#if len(sys.argv) != 2:
	#	print("Error! You must specify 1 argument: the quasisymmetry_out.*.nc file")
	#	exit(1)

	def toString(ncVar):
		temp = [c.decode('UTF-8') for c in ncVar]
		return (''.join(temp)).strip()

	filename = file#sys.argv[1]
	print("Reading filename "+filename)
	f = netcdf.netcdf_file(filename,mode='r',mmap=False)
	general_option = toString(f.variables['general_option'][()])
	if general_option != "single":
		print("Error! This script is designed for plotting single runs, but the quasisymmetry_out file you provided is a scan.")
		f.close()
		exit(1)

	nfp = f.variables['nfp'][()]
	B0 = f.variables['B0'][()]
	r = f.variables['r'][()]
	eta_bar = f.variables['eta_bar'][()]
	mpol = f.variables['mpol'][()]
	ntor = f.variables['ntor'][()]
	RBC = f.variables['RBC'][()]
	RBS = f.variables['RBS'][()]
	ZBC = f.variables['ZBC'][()]
	ZBS = f.variables['ZBS'][()]
	R0c = f.variables['R0c'][()]
	R0s = f.variables['R0s'][()]
	Z0c = f.variables['Z0c'][()]
	Z0s = f.variables['Z0s'][()]
	print("RBC.shape:",RBC.shape)

	phi = f.variables['phi'][()]
	X1c = f.variables['X1c'][()]
	Y1c = f.variables['Y1c'][()]
	Y1s = f.variables['Y1s'][()]
	sigma = f.variables['sigma'][()]
	curvature = f.variables['curvature'][()]
	torsion = f.variables['torsion'][()]
	elongation = f.variables['elongation'][()]
	elongation_in_Rz_plane = f.variables['elongation_in_Rz_plane'][()]
	modBinv_sqrt_half_grad_B_colon_grad_B = f.variables['modBinv_sqrt_half_grad_B_colon_grad_B'][()]

	#order_r_option = f.variables["order_r_option"][()]
	#order_r_option = ''.join(str(f.variables["order_r_option"][()]))
	order_r_option = toString(f.variables["order_r_option"][()])
	print("order_r_option:",order_r_option)
	order_r_squared = (order_r_option != 'r1' and order_r_option != 'r1_compute_B2')
	print("order_r_squared:",order_r_squared)
	order_r_cubed = (order_r_option != 'r1' and order_r_option != 'r1_compute_B2' and order_r_option != 'r2')
	print("order_r_cubed:",order_r_cubed)
	order_r1_compute_B2 = order_r_option == 'r1_compute_B2'

	if order_r1_compute_B2:
		B20 = f.variables['B20'][()]
		B2s_array = f.variables['B2s_array'][()]
		B2c_array = f.variables['B2c_array'][()]
		B02 = f.variables['B02'][()]

	if order_r_squared:
		B2s = f.variables['B2s'][()]
		B2c = f.variables['B2c'][()]
		B20_mean = f.variables['B20_mean'][()]
		B20 = f.variables['B20'][()]
		X20 = f.variables['X20'][()]
		X2s = f.variables['X2s'][()]
		X2c = f.variables['X2c'][()]
		Y20 = f.variables['Y20'][()]
		Y2s = f.variables['Y2s'][()]
		Y2c = f.variables['Y2c'][()]
		Z20 = f.variables['Z20'][()]
		Z2s = f.variables['Z2s'][()]
		Z2c = f.variables['Z2c'][()]
		r_singularity_vs_zeta = f.variables['r_singularity_vs_zeta'][()]
		try:
			r_singularity_basic_vs_zeta = f.variables['r_singularity_basic_vs_zeta'][()]
		except:
			r_singularity_basic_vs_zeta = r_singularity_vs_zeta # Old output files might not have this field
		for j in range(len(phi)):
			if r_singularity_vs_zeta[j] > 1.0e10:
				r_singularity_vs_zeta[j] = np.nan
			if r_singularity_basic_vs_zeta[j] > 1.0e10:
				r_singularity_basic_vs_zeta[j] = np.nan

	if order_r_cubed:
		X3s1 = f.variables['X3s1'][()]
		X3c1 = f.variables['X3c1'][()]
		X3s3 = f.variables['X3s3'][()]
		X3c3 = f.variables['X3c3'][()]
		Y3s1 = f.variables['Y3s1'][()]
		Y3c1 = f.variables['Y3c1'][()]
		Y3s3 = f.variables['Y3s3'][()]
		Y3c3 = f.variables['Y3c3'][()]
		Z3s1 = f.variables['Z3s1'][()]
		Z3c1 = f.variables['Z3c1'][()]
		Z3s3 = f.variables['Z3s3'][()]
		Z3c3 = f.variables['Z3c3'][()]
		try:
			B3s1 = f.variables['B3s1'][()]
			B3c1 = f.variables['B3c1'][()]
			B3s3 = f.variables['B3s3'][()]
			B3c3 = f.variables['B3c3'][()]
		except:
			B3s1 = phi * 0
			B3s3 = phi * 0
			B3c1 = phi * 0
			B3c3 = phi * 0

	f.close()

	my_xlim = [0,phi[-1]]

	N_theta = 150
	N_phi = 8
	Nfig  = N_phi

	theta1D = np.linspace(0,2*np.pi,N_theta)
	phi1D = np.linspace(0,2*np.pi/nfp,N_phi,endpoint=False)

	phi2D,theta2D = np.meshgrid(phi1D,theta1D)
	#print "theta2D:"
	#print theta2D

	R = np.zeros((N_theta,N_phi))
	z = np.zeros((N_theta,N_phi))
	for m in range(mpol+1):
		for jn in range(ntor*2+1):
			n = jn-ntor
			angle = m * theta2D - nfp * n * phi2D
			sinangle = np.sin(angle)
			cosangle = np.cos(angle)
			R += RBC[m,jn] * cosangle + RBS[m,jn] * sinangle
			z += ZBC[m,jn] * cosangle + ZBS[m,jn] * sinangle

	R0 = np.zeros(N_phi)
	z0 = np.zeros(N_phi)
	for n in range(len(R0c)):
		angle = nfp * n * phi1D
		sinangle = np.sin(angle)
		cosangle = np.cos(angle)
		R0 += R0c[n] * cosangle + R0s[n] * sinangle
		z0 += Z0c[n] * cosangle + Z0s[n] * sinangle

	#exit(0)

	import matplotlib.pyplot as plt

	fig = plt.figure(figsize=(16,7))
	fig.patch.set_facecolor('white')

	if order_r_cubed:
		numRows = 5
		numCols = 7
	elif order_r_squared:
		numRows = 4
		numCols = 5
	elif order_r1_compute_B2:
		numRows = 3
		numCols = 4
	else:
		numRows = 2
		numCols = 4
	plotNum = 1

	plt.subplot(numRows,numCols,plotNum)
	plotNum += 1
	plt.plot(phi,curvature)
	plt.title('curvature')
	plt.xlabel('$\phi$')
	plt.xlim(my_xlim)

	plt.subplot(numRows,numCols,plotNum)
	plotNum += 1
	plt.plot(phi,torsion)
	plt.title('torsion')
	plt.xlabel('$\phi$')
	plt.xlim(my_xlim)

	plt.subplot(numRows,numCols,plotNum)
	plotNum += 1
	plt.plot(phi,X1c)
	plt.title('X1c')
	plt.xlabel('$\phi$')
	plt.xlim(my_xlim)

	plt.subplot(numRows,numCols,plotNum)
	plotNum += 1
	plt.plot(phi,Y1s)
	plt.title('Y1s')
	plt.xlabel('$\phi$')
	plt.xlim(my_xlim)

	plt.subplot(numRows,numCols,plotNum)
	plotNum += 1
	plt.plot(phi,Y1c)
	plt.title('Y1c')
	plt.xlabel('$\phi$')
	plt.xlim(my_xlim)

	plt.subplot(numRows,numCols,plotNum)
	plotNum += 1
	plt.plot(phi,sigma)
	plt.title('sigma')
	plt.xlabel('$\phi$')
	plt.xlim(my_xlim)

	plt.subplot(numRows,numCols,plotNum)
	plotNum += 1
	plt.plot(phi,elongation,label='XY plane')
	plt.plot(phi,elongation_in_Rz_plane,label='Rz plane')
	plt.title('elongation')
	plt.legend(loc=0,fontsize=6)
	plt.xlabel('$\phi$')
	plt.xlim(my_xlim)

	plt.subplot(numRows,numCols,plotNum)
	plotNum += 1
	plt.plot(phi,modBinv_sqrt_half_grad_B_colon_grad_B)
	#plt.title('modBinv_sqrt_half_grad_B_colon_grad_B')
	plt.title('Frobenius')
	plt.xlabel('$\phi$')
	plt.xlim(my_xlim)

	if order_r1_compute_B2:
		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,B20)
		plt.title('B20')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,B2s_array)
		plt.title('B2s_array')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,B2c_array)
		plt.title('B2c_array')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,B02)
		plt.title('B_0^{(2)}')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

	if order_r_squared:
		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,X20)
		plt.title('X20')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,X2s)
		plt.title('X2s')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,X2c)
		plt.title('X2c')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,Y20)
		plt.title('Y20')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,Y2s)
		plt.title('Y2s')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,Y2c)
		plt.title('Y2c')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,Z20)
		plt.title('Z20')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,Z2s)
		plt.title('Z2s')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,Z2c)
		plt.title('Z2c')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,B20)
		plt.title('B20')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,r_singularity_vs_zeta,label='refined')
		plt.plot(phi,r_singularity_basic_vs_zeta,':k',label='unrefined')
		plt.legend(fontsize=6,loc=0)
		plt.title('r_singularity_vs_zeta')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

	if order_r_cubed:
		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,X3s1)
		plt.title('X3s1')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,X3c1)
		plt.title('X3c1')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,X3s3)
		plt.title('X3s3')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,X3c3)
		plt.title('X3c3')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,Y3s1)
		plt.title('Y3s1')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,Y3c1)
		plt.title('Y3c1')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,Y3s3)
		plt.title('Y3s3')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,Y3c3)
		plt.title('Y3c3')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,Z3s1)
		plt.title('Z3s1')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,Z3c1)
		plt.title('Z3c1')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,Z3s3)
		plt.title('Z3s3')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,Z3c3)
		plt.title('Z3c3')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,B3s1)
		plt.title('B3s1')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,B3c1)
		plt.title('B3c1')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,B3s3)
		plt.title('B3s3')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

		plt.subplot(numRows,numCols,plotNum)
		plotNum += 1
		plt.plot(phi,B3c3)
		plt.title('B3c3')
		plt.xlabel('$\phi$')
		plt.xlim(my_xlim)

	plt.tight_layout()

	titleString = "Plot generated by "+ os.path.abspath(file)
	plt.figtext(0.5,0.004,titleString,horizontalalignment='center',verticalalignment='bottom',fontsize=8)
	plt.figtext(0.5,0.997,'File = '+os.path.abspath(filename),horizontalalignment='center',verticalalignment='top',fontsize=8)
	plt.savefig(file+'_params.pdf', bbox_inches = 'tight', pad_inches = 0)

	###################################################################3

	fig = plt.figure(figsize=(16,7))
	fig.patch.set_facecolor('white')

	numRows = 3
	numCols = 3

	plt.subplot(numRows,numCols,1)
	for jphi in range(N_phi):
		plt.plot(R[:,jphi],z[:,jphi])
	plt.xlabel('R')
	plt.ylabel('z')
	plt.gca().set_aspect('equal',adjustable='box')
	Rfig=R
	Zfig=z
	R0fig=R0
	Z0fig=z0

	for jphi in range(N_phi):
		plt.subplot(numRows,numCols,jphi+2)
		for k in range(N_theta):
			plt.plot([R0[jphi],R[k,jphi]], [z0[jphi],z[k,jphi]],'-')
		plt.plot(R[:,jphi],z[:,jphi])
		plt.plot(R0[jphi],z0[jphi],'o')
		plt.xlabel('R')
		plt.ylabel('z')
		plt.gca().set_aspect('equal',adjustable='box')


	plt.tight_layout()

	titleString = "Plot generated by "+ os.path.abspath(__file__)
	plt.figtext(0.5,0.004,titleString,horizontalalignment='center',verticalalignment='bottom',fontsize=8)
	plt.figtext(0.5,0.997,'File = '+os.path.abspath(filename),horizontalalignment='center',verticalalignment='top',fontsize=8)

	########################################################
	# Now make 3D surface plot
	########################################################

	from mpl_toolkits.mplot3d import Axes3D
	from matplotlib import cm
	from numpy.matlib import repmat

	fig = plt.figure(figsize=(16,7))
	fig.patch.set_facecolor('white')

	# Zoom in:
	#factor = 1
	factor = 0.2
	fig.subplots_adjust(bottom=-factor+0.0,top=0.95+factor)

	N_theta = 40
	N_phi = 150

	theta1D = np.linspace(0,2*np.pi,N_theta)
	phi1D = np.linspace(0,2*np.pi,N_phi)

	phi2D,theta2D = np.meshgrid(phi1D,theta1D)
	#print "theta2D:"
	#print theta2D

	R = np.zeros((N_theta,N_phi))
	z = np.zeros((N_theta,N_phi))
	for m in range(mpol+1):
		for jn in range(ntor*2+1):
			n = jn-ntor
			angle = m * theta2D - nfp * n * phi2D
			sinangle = np.sin(angle)
			cosangle = np.cos(angle)
			R += RBC[m,jn] * cosangle + RBS[m,jn] * sinangle
			z += ZBC[m,jn] * cosangle + ZBS[m,jn] * sinangle

	x = R * np.cos(phi2D)
	y = R * np.sin(phi2D)
	B = B0 * (1 + r * eta_bar * np.cos(theta2D)) 
	if order_r_squared:
		#B20_2D = repmat(B20,N_theta,1)
		B20_interpolated = np.interp(phi1D,phi,B20,period=2*np.pi/nfp)
		#print "B20:",B20
		#print "B20_interpolated:",B20_interpolated
		B20_2D = repmat(B20_interpolated,N_theta,1)
		B += r * r * (B2s * np.sin(2*theta2D) + B2c * np.cos(2*theta2D) + B20_2D)
	# Rescale to lie in [0,1]:
	B_rescaled = (B - B.min()) / (B.max() - B.min())

	ax = fig.gca(projection='3d')
	#ax.set_aspect('equal')
	ax.plot_surface(x, y, z, facecolors = cm.viridis(B_rescaled), rstride=1, cstride=1, antialiased=False)

	max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0

	mid_x = (x.max()+x.min()) * 0.5
	mid_y = (y.max()+y.min()) * 0.5
	mid_z = (z.max()+z.min()) * 0.5
	ax.set_xlim(mid_x - max_range, mid_x + max_range)
	ax.set_ylim(mid_y - max_range, mid_y + max_range)
	ax.set_zlim(mid_z - max_range, mid_z + max_range)

	plt.show()

	print("Save PDF")
	fig = plt.figure(figsize=(14,7))
	fig.patch.set_facecolor('white')
	for jphi in range(Nfig):
		plt.plot(Rfig[:,jphi],Zfig[:,jphi])
		plt.plot(R0fig[jphi],Z0fig[jphi])
	plt.xlabel('R')
	plt.ylabel('Z')
	plt.gca().set_aspect('equal',adjustable='box')
	plt.savefig(file+'_surfs.pdf', bbox_inches = 'tight', pad_inches = 0)
