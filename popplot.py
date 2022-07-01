import matplotlib.pyplot as plt
import numpy as np 
from matplotlib import patches
from matplotlib.collections import PatchCollection
from matplotlib.path import Path
import os, subprocess
from tqdm import tqdm


def angletodeg(angle):
	return angle/np.pi*180

def plotbac(posx, posy, angle, length, radius, color = 'orange', force = None):	

	ax = plt.gca()
	cornerx = posx - length/2*np.cos(angle) + radius*np.sin(angle) # corner of the rectangle after rotation
	cornery = posy - length/2*np.sin(angle) - radius*np.cos(angle)

	vecnorm = np.array([-2*radius*np.sin(angle),2*radius*np.cos(angle)])
	vecpar =  np.array([length*np.cos(angle),length*np.sin(angle)])

	#corners of the central rectangle of the bacteria
	corner1 = np.array([cornerx,cornery])
	corner2 = corner1 + vecnorm
	corner3 = corner2 + vecpar
	corner4 = corner3 - vecnorm

	#edges of the bacterium
	vecedge = np.array([np.sqrt(2)*radius*np.cos(angle+np.pi*3/4),np.sqrt(2)*radius*np.sin(angle+np.pi*3/4)])
	edge1 = corner1 + vecedge
	edge2 = corner3 - vecedge

	#bezier points
	vecbez = np.array([radius*np.cos(angle+np.pi),radius*np.sin(angle+np.pi)])
	control1 = corner1 + vecbez
	control2 = corner2 + vecbez
	control3 = corner3 - vecbez
	control4 = corner4 - vecbez

	# plt.plot([corner1[0],corner2[0],corner3[0],corner4[0]],[corner1[1],corner2[1],corner3[1],corner4[1]],'ob')
	# plt.plot([corner1[0],corner1[0]+vecnorm[0]],[corner1[1],corner1[1]+vecnorm[1]],'-r')
	# plt.plot([corner1[0],corner1[0]+vecpar[0]],[corner1[1],corner1[1]+vecpar[1]],'-g')

	verts = [
	   corner1,    # P0
	   control1,    # P1
	   edge1,    # P2
	   control2,  #these are the vertices for the 
	   corner2,  #Bezier curve
	   corner3,  # P3
	   control3,
	   edge2,
	   control4,
	   corner4,
	   corner1     #and back to P0
	]

	codes = [
	    Path.MOVETO,  # move to corner 1
	    Path.CURVE3,  #line to P2
	    Path.CURVE3,  #control point one for the Bezier curve
	    Path.CURVE3,  #control point two
	    Path.CURVE3,  #end point for the Bezier curve
	    Path.LINETO,
	    Path.CURVE3,  #line to P2
	    Path.CURVE3,  #control point one for the Bezier curve
	    Path.CURVE3,  #control point two
	    Path.CURVE3,  #end point for the Bezier curve
	    Path.LINETO	      #and back to P0
	 ]

	path = Path(verts, codes)
	patch = patches.PathPatch(path, facecolor='orange', edgecolor = 'k', lw=2, zorder = 0)
	ax.add_patch(patch)

	# plt.plot([posx],[posy],'bo')

	if force:
		plt.plot([],[])

foldername_output = 'output/'
foldername_colony = foldername_output+'colony/'
foldername_diffusible = foldername_output+'diffusible/'
files_colony = os.listdir(foldername_colony)
files_diffusible = os.listdir(foldername_colony)
timelist = []
for filename in files_colony:
#	if filename[0:10] == 'population'
	if ((filename[0:6] == 'colony') and (filename[-4:]=='.out')):
		timelist.append(filename[7:-4])
print('Unsorted list',timelist)
timelist = sorted(timelist,key = lambda x: float(x))
# print('Sorted list',timelist)
	
print("Creating Frames...")
for it,timepoint in enumerate(tqdm(timelist)):
	file_bac = np.loadtxt(foldername_colony+'colony_'+timepoint+'.out')
	for bac in file_bac:
		plotbac(bac[2],bac[3],bac[4],bac[5],0.2)
	ax = plt.gca()

	vmin = 0
	vmax = 5
	maxalpha = 0.8 # most opaque value of the morphogen
	file_diff = np.loadtxt(foldername_diffusible+'diffusible_0_'+timepoint+'.out')
	alphas_diff = file_diff
	alphas_diff[alphas_diff>vmax] = vmax
	alphas_diff[alphas_diff<vmin] = vmin
	alphas_diff = (alphas_diff-vmin)/(vmax-vmin)*maxalpha
	plt.imshow(file_diff,cmap = 'Blues',vmin = 0, vmax = 5,alpha = alphas_diff,
			   zorder = 10, extent = [-5,5,-5,5])

	ax.set_aspect('equal')
	plt.xlim(-5,5)
	plt.ylim(-5,5)

	plt.title('t = {}'.format(timepoint))


	plt.savefig(foldername_output + "/file%02d.png" % it)


	plt.clf()
os.chdir(foldername_output)
subprocess.call([
	'ffmpeg', '-framerate', '8', '-i', 'file%02d.png', '-r', '30', '-pix_fmt', 'yuv420p',
	'output.mp4'])
for filename in os.listdir():
	if(filename[-4:]=='.png'):
		os.remove(filename)
	# plt.show()	

