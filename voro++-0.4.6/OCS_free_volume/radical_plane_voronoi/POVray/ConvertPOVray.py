import numpy
import sys
import time

#------------------------------------------------------------------------------#


def GetCoordinates(lineData):
	# remove coordinates from string object 
	startIdx = lineData.find('<')
	endIdx = lineData.find('>')
	coords = lineData[startIdx+1 : endIdx-1]
	coords = coords.split(',')
	coords = map(float,coords)
	return coords

def GetCoordsCylinder(lineData):
	startIdx = lineData.find('<')
	endIdx = lineData.find('>')
	coords1 = lineData[startIdx+1 : endIdx-1]
	coords1 = coords1.split(',')
	coords1 = map(float,coords1)

	startIdx2 = lineData.find('<', endIdx+1)
	endIdx2 = lineData.find('>', endIdx+1)
	coords2 = lineData[startIdx2+1 : endIdx2-1]
	coords2 = coords2.split(',')
	coords2 = map(float,coords2)

	return coords1,coords2

def GetRadius(lineData):
	# remove radius from string object
	idx = lineData.find('>')
	radius = float(lineData[idx+2:-2])
	return radius

def GetLineID(lineData):
	idx = lineData.rfind(' ')
	ID = lineData[idx+1:]
	ID = int(ID.strip())
	return ID

def ParseParticleString(lineData):
	# remove coordinates from string object 
	startIdx = lineData.find('<')
	endIdx = lineData.find('>')
	coords = lineData[startIdx+1 : endIdx-1]
	coords = coords.split(',')
	coords = map(float,coords)

	# remove radius from string object
	raduis = float(lineData[endIdx+2:-2])

	return coords,raduis

def FindCenterPoint(dataPoints):
	xmin = min(dataPoints[:,0])
	xmax = max(dataPoints[:,0])

	ymin = min(dataPoints[:,1])
	ymax = max(dataPoints[:,1])

	zmin = min(dataPoints[:,2])
	zmax = max(dataPoints[:,2])

	return ((xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2)

def TranslateCoordsList(center,dataPoints):
	for i in range(len(dataPoints)):
		dataPoints[i][0] -= center[0]
		dataPoints[i][1] -= center[1]
		dataPoints[i][2] -= center[2]
	return dataPoints

def TranslateCoords(center,dataPoint):
	dataPoint[0] -= center[0]
	dataPoint[1] -= center[1]
	dataPoint[2] -= center[2]
	return dataPoint

def TransferPointData(inputfile,outputfile):
	ID = 0
	# load atom coordinates and find center point to translate simulation cell to origin
	dataPoints = []	
	fin = open(inputfile)
	for line in fin:
		if line[0] != '/':
			dataPoints.append(GetCoordinates(line))
	fin.close()

	# find the center of the simutlation cell and translate
	center = FindCenterPoint(numpy.array(dataPoints))
	dataPoints = TranslateCoordsList(center,dataPoints)

	# write out new atom coordinates, radius and color
	n = 0
	fin = open(inputfile)
	fout = open(outputfile, 'w')
	n = 0
	for line in fin:
		if line[0] != '/':
			coords = dataPoints[n]
			n += 1
			radius = GetRadius(line)

			if radius == 1.52:
				fout.write('sphere{ <%.4f, %.4f, %.4f> %.3f texture{tO}}\n' % (coords[0],coords[1],coords[2],radius*.75))
			elif radius == 1.7:
				fout.write('sphere{ <%.4f, %.4f, %.4f> %.3f texture{tC}}\n' % (coords[0],coords[1],coords[2],radius*.75))
			elif radius == 2.1: 
				fout.write('sphere{ <%.4f, %.4f, %.4f> %.3f texture{tSi}}\n' % (coords[0],coords[1],coords[2],radius*.75))

			# ''' conditions for voro++ example data'''
			# if radius < 0.4:
			# 	fout.write('sphere{ <%.4f, %.4f, %.4f> %.3f texture{tO}}\n' % (coords[0],coords[1],coords[2],radius))
			# elif radius > 0.4 and radius < 0.6:
			# 	fout.write('sphere{ <%.4f, %.4f, %.4f> %.3f texture{tC}}\n' % (coords[0],coords[1],coords[2],radius))
			# elif radius > 0.6: 
			# 	fout.write('sphere{ <%.4f, %.4f, %.4f> %.3f texture{tSi}}\n' % (coords[0],coords[1],coords[2],radius))

		else:
			# n += 1
			# if n > 300:
			# 	''' only write out 150 atoms '''
			# 	ID = GetLineID(line)
			# 	break
			''' else just write the line'''
			fout.write('%s'%line)

	fin.close()
	fout.close()
	return center, ID 

def TransferCellData(inputfile,outputfile,center,ID):
	fin = open(inputfile)
	fout = open(outputfile,'w')
	for line in fin:
		if line[0] != '/':
			if line[0][0] == 's':
				''' case for circle '''
				coords = GetCoordinates(line)
				coords = TranslateCoords(center,coords)
				fout.write('sphere{ <%.4f, %.4f, %.4f> r}\n' % (coords[0],coords[1],coords[2]))

			else:
				'''case for cylinder'''
				coords1,coords2 = GetCoordsCylinder(line)
				coords1 = TranslateCoords(center,coords1)
				coords2 = TranslateCoords(center,coords2)
				fout.write('cylinder{ <%.4f, %.4f, %.4f>, <%.4f, %.4f, %.4f>, r}\n' % 
					(coords1[0],coords1[1],coords1[2],coords2[0],coords2[1],coords2[2]))
		
		else:
			# ''' stop writing once last atom ID is reached '''
			# if GetLineID(line) == ID:
			# 	break
			fout.write('%s'%line)


#------------------------------------------------------------------------------#

if len(sys.argv) < 3:
	print 'Usage:'
	print '  python %s <particle file> <voronoi cell file>'%sys.argv[0]
	exit()

t0 = time.time()

pfile = sys.argv[1]
vfile = sys.argv[2]

outputfile_p = pfile[:-4]+'out.pov'
center,ID = TransferPointData(pfile,outputfile_p)

outputfile_v = vfile[:-4]+'out.pov'
TransferCellData(vfile,outputfile_v,center,ID)


print 'Time elapsed %.4f seconds' %(time.time() - t0)
