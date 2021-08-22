#!/usr/bin/env python3

def main(file,qvfilename):

	import matplotlib.pyplot as plt

	filename = file

	token = open('neo_out.'+qvfilename,'r')
	linestoken=token.readlines()
	column_number = 1
	eps_eff=[]
	s_radial=[]
	for x in linestoken:
	    s_radial.append(float(x.split()[0])/251)
	    eps_eff.append(float(x.split()[1]))
	token.close()
	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.plot(s_radial,eps_eff)
	ax.set_yscale('log')
	plt.xlabel(r'$s=\psi/\psi_a$', fontsize=12)
	plt.ylabel(r'$\epsilon_{eff}$', fontsize=12)
	fig.savefig('neo_out_'+qvfilename+'.pdf', dpi=fig.dpi)
