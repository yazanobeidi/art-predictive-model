#!/usr/bin/python
""" Target Volume Predictive Model for Adaptive Radation Therapy
This file contains class definitions, parsing routines and prediction handling 
for the Kalman Filter based Predictive Model of ART Data.

This file depends on input from plan.roi and a file named artplan.roi to dump 
model output. Also for debugging, a file named output.txt is created.
"""

def kalmanf(zk, A, xk, Pk, B, uk, wk, Q, R, H, output, t):
    """Kalman Filter Function
    :param: zk: raw data
    :param: A: state transition coefficient matrix...list object for no reason
    :param: xk: state estimate vector
    :param: Pk: prediction covariance matrix
    :param: B: control coefficient matrix
    :param: uk: control signal vector
    :param: wk: process noise (has covariance Q)
    :param: Q: process covariance
    :param: R: measurement covariance
    :param: H:  transformation matrix
    :param: t: initially size of zk, decrements as data is processed
    """
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

class Model(object):
    """
    Model class is a container which holds a series of 'fractions'
    """
    def __init__(self):
        self.fractions = []

    def add_fraction(self, name, location, num):
        # each new fraction gets added to list
        self.fractions.append(Fraction(name, location, num)) 

    def get_num_fractions(self):
        return len(self.fractions)

    def printAll(self):
        """
        This function runs through each fraction and returns name, 
        followed by all its points
        """
        # no with statements to fit ART requirements
        fout = open("output.txt", 'w+')
        for fraction in self.fractions:
            print fraction.name + ':'
            print str(fraction.returnAll())
            fout.write(fraction.name + ':\n')
            fout.write(str(fraction.returnAll()) + '\n')
        fout.close()

    def checkSize(self):
        total = sum([fraction.numPoints for fraction in self.fractions])
        return True if \
            total == len(self.fractions[0].x) * len(self.fractions) else False

    def predict(self, predictiveModel):
        """
        This member function applies Kalman Filter on each dimension of each 
        point through all available fractions.
        """
        A = [1]
        B = 0
        uk = 0
        Pi = 1
        wk = 0
        Q = 0 # covariance of wk
        R = 0.1
        H = 1
        # Through all points, through all fractions, append to time series
        for n in range(0, self.fractions[0].numPoints):
            Tn_x = []
            Tn_y = []
            Tn_z = []
            output_x = []
            output_y = []
            output_z = []
            for m in range(0, len(self.fractions)):
                Tn_x.append(self.fractions[m].x[n])
                Tn_y.append(self.fractions[m].y[n])
                Tn_z.append(self.fractions[m].z[n])
            t = len(Tn_x)
            # Apply Kalman Filter on time series
            kalmanf(Tn_x, A, Tn_x[0], Pi, B, uk, wk, Q, R, H, output_x, t)
            kalmanf(Tn_y, A, Tn_y[0], Pi, B, uk, wk, Q, R, H, output_y, t)
            kalmanf(Tn_z, A, Tn_z[0], Pi, B, uk, wk, Q, R, H, output_z, t)
            predictiveModel.fractions[0].add(str(output_x[-1])+' '+str(output_y[-1])+' '+str(output_z[-1]))
            
class Fraction(object):
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
with open('data/plan.roi') as plan: # with statement does not work in python 2.4.6, use plan = open('plan.roi')
    for num, line in enumerate(plan, 1):
        if "CTV12" in line and "org" not in line and "//" not in line and "SC" not in line and "IC" not in line:
            print "CONTOUR:" + line                         #debugging
            print "Found on line: " + str(num)      #debugging
            model.add_fraction(line.split()[1], plan.tell(), num)
            
# This second pass starts at each fraction file location and runs through plan.roi
# until it reaches 'number_of_vertices'. The value is recorded.
foundNumPoints = False
with open('data/plan.roi') as plan:
    # Iterate across all fractions
    for f_index in range(0, model.get_num_fractions()):
        # Start at fractional location in plan
        plan.seek( model.fractions[f_index].getLocation(), 0) # 0=seek from beginning
        # Iterate through plan
        #for line in plan: ... removed as mixing iteration and readline causes errors
        while True:
            line = plan.readline()
            if not line: break
            if 'number_of_vertices' in line and not foundNumPoints:
                model.fractions[f_index].setNumPoints(int(line.split()[2][:-1]))
                foundNumPoints = True
            if 'vertices={' in line:
                foundNumPoints = False
                break
        # Iterate for however number of points in fraction and consume line
        for point_index in range(0,model.fractions[f_index].getNumPoints()):
            model.fractions[f_index].add(plan.readline())

# At this point, all points are loaded into model Model object
# To print points:
# (To adjust print to file/console, change object method in class def)
# model.printAll()

# create new model instance and add one fraction to hold predictive model output
predictiveModel = Model()
predictiveModel.add_fraction('Model', 0, 1)
# Adjust number of points
predictiveModel.fractions[0].setNumPoints(model.fractions[0].getNumPoints())
# apply predictive model with parsed input
if model.checkSize(): model.predict(predictiveModel)
else: print 'Size mismatch'

# At this point, predictive model is now loaded into predictiveModel Model object
# Can take a look with:
#print predictiveModel.fractions[0].x
#print predictiveModel.fractions[0].y
#print predictiveModel.fractions[0].z
predictiveModel.printAll()

# ETL predictive_model to file by copying a contour, 
# predictive_model.toPinnacle()
with open('data/plan.roi', 'r+') as plan:
    newfile = ''
    atPoints = 0
    point_index = 0
    firstTime = True
    for num, line in enumerate(plan,1):
        if atPoints > 0: 
            newfile += predictiveModel.fractions[0].returnPoint(point_index)
            point_index += 1
            atPoints -= 1
        else:
            newfile += line
        if 'vertices={' in line and firstTime:
            atPoints = predictiveModel.fractions[0].getNumPoints()
            firstTime = False
    newplan = open('newartplan.roi', 'w+')
    newplan.write(newfile)
    newplan.close()

# At this point newartplan.roi contains kalman filter output contour