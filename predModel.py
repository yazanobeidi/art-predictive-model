#!/usr/bin/python
# This file contains class definitions, parsing routines and prediction handling 
# for the Kalman Filter based Predictive Model of ART Data.
# This file depends on input from plan.roi and a file named artplan.roi to dump 
# model output. Also for debugging, a file named output.md is created.

# For debugging purposes, opens output filestream to dump point data out
fout = open("output.md", 'w')

def kalmanf(zk, A, xk, Pk, B, uk, wk, Q, R, H, output, t):
	# Kalman Filter Function
	# Inputs:
	# * zk : raw data
	# * A : state transition coefficient matrix...list object for no reason
	# * xk : state estimate vector
	# * Pk : prediction covariance matrix
	# * B : control coefficient matrix
	# * uk : control signal vector
	# * wk : process noise (has covariance Q)
	# * Q : process covariance
	# * R : measurement covariance
	# * H :  transformation matrix
	# * t : initially size of zk, decrements as data is processed
	# Recursively outputs to output
	zk_index = len(zk) - t + 1
	# Time update (prediction)
	xk = float(A[0])*xk + B*uk + wk; # Predict state
	zk_pred = H*xk; # Predict measurement (comment out if alt est update used)
	Pk = float(A[0])*Pk*float(A[0]) + Q; # Predict error covariance
	# Measurement update  (correction)
	vk = (float(zk[zk_index-1]) - zk_pred); # Measurement residual
	S = (H*Pk*float(H) + R); # Prediction covariance
	Kk = Pk*float(H) / S ; # Compute Kalman Gain
	xk = xk + Kk * (vk); # Update estimate with gain * residual
	Pk = (1 - Kk*H)*Pk; # Update error covariance
	t -= 1
	output.append(xk)
	if t <= 0:
		return
	else:
		kalmanf(zk, A, xk, Pk, B, uk, wk, Q, R, H, output, t)

class Model:
# Model class is an abstracted container which holds a series of 'fractions'	
	def __init__(self):
		# initiate list of fraction objects
		self.fraction = []
	def addFraction(self, name, location, num):
		# each new fraction gets added to list
		self.fraction.append(Fraction(name, location, num)) 
	def getNumFractions(self):
		return len(self.fraction)
	def printAll(self):
		# This function runs through each fraction and returns name, then all points
		for i in range(0, self.getNumFractions()):
			print self.fraction[i].name + ':'
			print str(self.fraction[i].returnAll())
			# fout only works if fout = open("output.md", 'w') is declared at top
			fout.write(self.fraction[i].name + ':\n')
			fout.write(str(self.fraction[i].returnAll()) + '\n')
	def checkSize(self):
		total = int()
		for i in range(0, len(self.fraction)):
			total += int(self.fraction[i].numPoints)
		if total == len(self.fraction[0].x) * len(self.fraction): return True
		else: return False
	def predict(self, predictiveModel):
		# This member function applies Kalman Filter on each dimension of each point
		# through all available fractions. The output is passed onto
		#   Logic:
		#   let Pn = (xn,yn,zn) for all n > 0
		#   let Fm = [Pn] for all m > 0
		#   Fm = F(m-1) +- delta 
		#   Define Tn(Pn,Dd) = [F0(P0,D0);... Fm(Pn,Dn)] for all DdCPn, d = (1,2,3)
		#   Algorithm:
		#   output_n_d = kalman(Tn(Pn,Dd))
		#   model(n,d) = output_n_d(end) where model is a set of n points 
		# First initialize predictiveModel with fractions
		#for m in range(0, len(self.fraction)):
			#predictiveModel.addFraction(self.fraction[m].name, self.fraction[m].location, self.fraction[m].num)	
		A = [1]
		B = 0
		uk = 0
		Pi = 1
		wk = 0
		Q = 0 # covariance of wk
		R = 0.1
		H = 1
		# Through all points, through all fractions, append to time series
		for n in range(0, self.fraction[0].numPoints):
			Tn_x = []
			Tn_y = []
			Tn_z = []
			output_x = []
			output_y = []
			output_z = []
			for m in range(0, len(self.fraction)):
				Tn_x.append(self.fraction[m].x[n])
				Tn_y.append(self.fraction[m].y[n])
				Tn_z.append(self.fraction[m].z[n])
			t = len(Tn_x)
			# Apply Kalman Filter on time series
			kalmanf(Tn_x, A, Tn_x[0], Pi, B, uk, wk, Q, R, H, output_x, t)
			kalmanf(Tn_y, A, Tn_y[0], Pi, B, uk, wk, Q, R, H, output_y, t)
			kalmanf(Tn_z, A, Tn_z[0], Pi, B, uk, wk, Q, R, H, output_z, t)
			predictiveModel.fraction[0].add(str(output_x[-1])+' '+str(output_y[-1])+' '+str(output_z[-1]))
			
class Fraction:
	# Each fraction contains lists x, y, z, concatenates upon return
	def __init__(self, name, location, num):
		self.name = name
		self.location = location
		self.lineNum = num
		self.numPoints = 0
		self.x = []
		self.y = []
		self.z = []	
	def add(self, data): 
		# Adds data by appending onto appropriate sublist
		self.x.append(float(data.split()[0]))
		self.y.append(float(data.split()[1]))
		self.z.append(float(data.split()[2]))
	def setNumPoints(self, points):
		self.numPoints = points
	def getLocation(self):
		return self.location
	def getNumPoints(self):
		return self.numPoints
	def returnAll(self): 
		# This method goes down sublists, appends and returns x+y+z in a string
		self.all = ""
		for i in range(0,self.numPoints):
			self.all += str(self.x[i])+' '+str(self.y[i])+' '+str(self.z[i])+ '\n'
		return self.all
	def returnPoint(self, index):
		return str(str(self.x[index])+' '+str(self.y[index])+' '+str(self.z[index])+'\n')

# Create model object to hold data
model = Model()

# This section of code runs through plan.roi until it reaches CTV12 contours at which
# point it creates a new fraction in Model object and records the file location
with open('plan.roi') as plan: # with statement does not work in python 2.4.6, use plan = open('plan.roi')
	for num, line in enumerate(plan, 1):
		if "CTV12" in line and "org" not in line and "//" not in line and "SC" not in line and "IC" not in line:
			print "CONTOUR:" + line							#eventually remove this
			print "Found on line: " + str(num)		#eventually remove this
			model.addFraction(line.split()[1], plan.tell(), num)
			
# This second pass starts at each fraction file location and runs through plan.roi
# until it reaches 'number_of_vertices'. The value is recorded.
foundNumPoints = False
with open('plan.roi') as plan:
	# Iterate across all fractions
	for f_index in range(0, model.getNumFractions()):
		# Start at fractional location in plan
		plan.seek( model.fraction[f_index].getLocation(), 0) # 0=seek from beginning
		# Iterate through plan
		#for line in plan: ... removed as mixing iteration and readline causes errors
		while True:
			line = plan.readline()
			if not line: break
			if 'number_of_vertices' in line and not foundNumPoints:
				model.fraction[f_index].setNumPoints(int(line.split()[2][:-1]))
				foundNumPoints = True
			if 'vertices={' in line:
				foundNumPoints = False
				break
		# Iterate for however number of points in fraction and consume line
		for point_index in range(0,model.fraction[f_index].getNumPoints()):
			model.fraction[f_index].add(plan.readline())

# At this point, all points are loaded into model Model object
# To print points:
# (To adjust print to file/console, change object method in class def)
# model.printAll()

# create new model instance and add one fraction to hold predictive model output
predictiveModel = Model()
predictiveModel.addFraction('Model', 0, 1)
# Adjust number of points
predictiveModel.fraction[0].setNumPoints(model.fraction[0].getNumPoints())
# apply predictive model with parsed input
if model.checkSize(): model.predict(predictiveModel)
else: print 'Size mismatch'

# At this point, predictive model is now loaded into predictiveModel Model object
# Can take a look with:
#print predictiveModel.fraction[0].x
#print predictiveModel.fraction[0].y
#print predictiveModel.fraction[0].z
predictiveModel.printAll()

# save predictive_model to file by copying a contour, 
# predictive_model.toPinnacle()
with open('artplan.roi', 'r+') as plan:
	newfile = ''
	atPoints = 0
	point_index = 0
	firstTime = True
	for num, line in enumerate(plan,1):
		if atPoints > 0: 
			newfile += predictiveModel.fraction[0].returnPoint(point_index)
			point_index += 1
			atPoints -= 1
		else:
			newfile += line
		if 'vertices={' in line and firstTime:
			atPoints = predictiveModel.fraction[0].getNumPoints()
			firstTime = False
	newplan = open('newartplan.roi', 'w')
	newplan.write(newfile)
	newplan.close()

# At this point newartplan.roi contains kalman filter output contour

fout.close() # close debugging filestream
