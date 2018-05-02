# -*- coding: utf-8 -*-
"""
PyLocalCsys was written by Boris Burgarella
the objective of this module is to generate local coordinate system in the finite element soft abaqus
This is usually useful to make composite calculations.

  Copyright (C) <2018>  <Boris Burgarella>
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
   
"""


from abaqus import *
from abaqusConstants import *
import regionToolset
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import sketch
import part
import math
import numpy as np
import sys
import os
from abaqus import getInputs
from abaqus import getWarningReply, YES, NO





##############################
#		User options		 #
##############################

#Number of part to be treated
NbParts = float(getInput('How many parts would you like to generate local csys on ?'))
#Table of the part names:
PartList = getInput('Please enter the name of the parts separated by a comma "," ')
TblParts = PartList.split(',')

##################################
#		class definitions		 #
##################################

os.system('clear')
print "Start"
runstate = 1

class CSYS:
	"""
	This class is used to define a coordinate system in abaqus, it is initialized by a set of coordinates:
	P0 corresponds to the origin, P1, the end of the X axis, P2, the end of the Y axis, and P3 the end of the Z axis.
	
	Thanks to the CalculateAxes method, only P0 and P1 are required to generate a full coordinate system.
	This is subject to some hypotheses:
	
	- The abaqus mesh shall be generated using the Bottom-up technique in order to make sure that the nodes position match with the one used to validate this program.
	- The rotation around the X axis are not possible to take into accound with this program, This limitation is due to the fact that it is impossible to evaluate such a rotation without knowing more than just
	the node position (except for perfectly parrallellepipedic elements)
	
	This python script shall not be used in any other case, if you do, please validate the local csys by using the abaqus option "Plot -> Material Orientation" in the result tab of CAE
	
	"""
	def ___init___(self,P0 = [0,0,0],P1 = [0,0,0],P2=[0,0,0],P3=[0,0,0]):

		self.P0 = P0
		self.P1 = P1
		self.P2 = P2
		self.P3 = P3
		#self.Ax = [0,0,0]
		#self.Ay = [0,0,0]
		#self.Az = [0,0,0]
		
	def CalculateAxes(self):
		"""
		The axes are defined by take P0 and P1 to define the X axis direction, The vector given by [P0-P1] is then normed
		to have the first axis of the base.
		
		Then, the displacement between the end of usual X axis [1,0,0] and the end of new one is calculated, and allows the calculation
		of the Y axis by vector additionalRotationField
		
		Finally, the Z axis is calculated thanks to the vectorial product.
		
		In order to avoid compulational errors, the found Y and Z axes are then normed (they usually have a magnitude of 0.99 or 1.01 due to float point issues)
		
		"""
		# definition of the X axis
		self.Ax = Norm(TupleSoustraction(self.P1,self.P0))
		Diff = TupleSoustraction([1,0,0],self.Ax,)
		# definition of the Y axis
		self.Ay = Norm(TupleAddition([0,1,0],Diff))
		#definition of the Z axis
		self.Az = Norm(VectorialProduct(self.Ax,self.Ay))
		self.Ay = Norm(VectorialProduct(self.Az,self.Ax))
		
		#print "-----------------------------"
		#print self.Ax
		#print self.Ay
		#print self.Az
		if self.ValidOrthoNormal() == 0:
			print "Error in Csys Generation, Local Csys is not an orthonormal base"
			return 0
		else:
			return 1
		
		
	def ValidOrthoNormal(self):
		"""
		This function simply checks that each vector in the Csys has a magnitude of 1 and that they are all orthogonal vectors
		"""
	
		Valid = 0

		if round(NormVal(self.Ax),5) != 1.0:
			Valid = 0
			print "Ax norm is not 1"
			print "Ax norm is "+str(NormVal(self.Ax))
			return Valid
		if  round(NormVal(self.Ay),5) != 1.0:
			Valid = 0
			print "Ay Norm is not 1"
			print NormVal(self.Ay)
			return Valid
		if round(NormVal(self.Az),5) != 1.0:
			Valid = 0
			print NormVal(self.Az)
			print "Az Norm is not 1"
			return Valid
		if round(abs(ScalarProduct(self.Ax,self.Ay)),5) >= 0.0001:
			Valid = 0
			print "Ax and Ay are not orthogonal vectors"
			print round(abs(ScalarProduct(self.Ax,self.Ay)),5)
			return Valid	
		if round(abs(ScalarProduct(self.Ax,self.Az)),5) >= 0.0001:
			Valid = 0
			print "Ax and Az are not orthogonal vectors"
			print round(abs(ScalarProduct(self.Ax,self.Az)),5)
			return Valid				
		if round(abs(ScalarProduct(self.Ay,self.Az)),5) >= 0.0001:
			Valid = 0
			print "Ay and Az are not orthogonal vectors"
			print round(abs(ScalarProduct(self.Ay,self.Az)),5)
			return Valid	
		else:
			#print "Ok - The local Csys is orthonormal"
			return 1
			
class node:
	"""
	A simple class used to store nodes,
	 - Label -> the abaqus node label / number
	 - coordinates -> A 3x1 vector used to store the node coordinates
	 
	"""

	def __init__(self,label):
		self.label = label
		self.coordinates = [0,0,0]
		
	def QueryCoordinates(self,PartName):
		self.coordinates = mdb.models['Model-1'].parts[PartName].nodes[self.label].coordinates
		
class element:
	"""
	
	A simple class used to store the information about the elements
	nodetable -> List of node objects that gathers all the nodes connected to the elements
	Label -> Abaqus Element number / label
	LocalCsys -> Csys associated with this element, generated by the GenerateLocalCsys method
	
	"""
	def __init__(self,ElementLabel,connectivity,PartName,ElementType="C3D8R"):
		self.nodetable = []
		self.Label = ElementLabel
		self.LocalCsys = CSYS()
		for i in connectivity:
			self.nodetable.append(node(i))
			self.nodetable[-1].QueryCoordinates(PartName)
			
	def GenerateLocalCsys(self):
		global runstate
		self.LocalCsys.P0 = self.nodetable[0].coordinates
		self.LocalCsys.P1 = self.nodetable[4].coordinates
		self.LocalCsys.P2 = self.nodetable[1].coordinates
		self.LocalCsys.P3 = self.nodetable[3].coordinates
		runstate = self.LocalCsys.CalculateAxes()

		
class Mesh:
	"""
	
	This class is used to gather all the information and methods applied on the whole mesh
	PartName -> Gives the name of the part the mesh is associated to
	ElemAbaqusFormLis -> List of the Raw elements directly imported  from abaqus
	ElemPythonList -> List of Element objects converted from the ElemAbaqusFormLis
	
	The runstate variable is used here to check if everything runs fine, you should not have issues with this variable if you don't modify this program
	
	"""
	def __init__(self,PartName):
		self.PartName = PartName
		self.ElemAbaqusFormLis = mdb.models['Model-1'].parts[PartName].elements
		self.ElemPythonList = []
		for i in self.ElemAbaqusFormLis:
			if runstate == 0:
				break
			self.ElemPythonList.append(element(i.label,i.connectivity,self.PartName))
			self.ElemPythonList[-1].GenerateLocalCsys()
		if runstate == 1:
			print "All the Csys were generated and validated"
			
	def ExportLocalCsys(self):
		"""
		This method simply calls the abaqus command to generate the table of local csys 
		Please refer to the abaqus documentation about discrete fields for more information
		
		"""
		
		strlabels = "("
		for i in range(len(self.ElemPythonList)):
			if i == len(self.ElemPythonList)-1:
				strlabels += str(self.ElemPythonList[i].Label)+")"
			else:
				strlabels += str(self.ElemPythonList[i].Label)+", "
		#print strlabels
		strcoord = "("
		for i in range(len(self.ElemPythonList)):
			if i == len(self.ElemPythonList)-1:
				strcoord += str(self.ElemPythonList[i].LocalCsys.Ax[0])+","+str(self.ElemPythonList[i].LocalCsys.Ax[1])+","+str(self.ElemPythonList[i].LocalCsys.Ax[2])+","+str(self.ElemPythonList[i].LocalCsys.Ay[0])+","+str(self.ElemPythonList[i].LocalCsys.Ay[1])+","+str(self.ElemPythonList[i].LocalCsys.Ay[2])+")"
			else:
				strcoord += str(self.ElemPythonList[i].LocalCsys.Ax[0])+","+str(self.ElemPythonList[i].LocalCsys.Ax[1])+","+str(self.ElemPythonList[i].LocalCsys.Ax[2])+","+str(self.ElemPythonList[i].LocalCsys.Ay[0])+","+str(self.ElemPythonList[i].LocalCsys.Ay[1])+","+str(self.ElemPythonList[i].LocalCsys.Ay[2])+","			
		#print strcoord
		StrToEvalutate1 = "mdb.models['Model-1'].DiscreteField(data=(('', 6, "+strlabels+", "+strcoord+"), ), dataWidth=6, defaultValues=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0),description='', fieldType=ORIENTATION, location=ELEMENTS, name='DiscField_"+str(self.PartName)+"', orientationType=CARTESIAN, partLevelOrientation=True)"
		eval(StrToEvalutate1)
		DiscFieldName = "DiscField_"+str(self.PartName)
		mdb.models['Model-1'].parts[self.PartName].MaterialOrientation(
			additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
			, axis=AXIS_1, fieldName=DiscFieldName, localCsys=None, 
			orientationType=FIELD, region=Region(
			cells=mdb.models['Model-1'].parts[self.PartName].cells.getSequenceFromMask(
			mask=('[#1 ]', ), )), stackDirection=STACK_3)
	
#####################################
#		function definitions 		#
#####################################	
		
def TupleSoustraction(T1,T2):
	"""
	Simple tuple soustraction function
	"""

	return [T1[0]-T2[0],T1[1]-T2[1],T1[2]-T2[2]]

def TupleAddition(T1,T2):
	"""
	Simple tuple addition function
	"""

	return [T1[0]+T2[0],T1[1]+T2[1],T1[2]+T2[2]]
	
def ScalarProduct(V1,V2):
	"""
	Simple scalar product function
	"""

	return V1[0]*V2[0]+V1[1]*V2[1]+V1[2]*V2[2]
	
def AngleVector(V1,V2):
	"""
	Calculate the angle between two vectors, not used anymore
	"""
	return math.acos((ScalarProduct(V1,V2))/(NormVal(V1)*NormVal(V2)))

def VectorialProduct(V1,V2):
	"""
	Simple Vectorial Product function
	"""
	return [V1[1]*V2[2]-V1[2]*V2[1],V1[2]*V2[0]-V1[0]*V2[2],V1[0]*V2[1]-V1[1]*V2[0]]	
	
def NormVal(Vect):
	"""
	Calculate the norm of a vector, returns a scalar
	"""
	X = Vect[0]
	Y = Vect[1]
	Z = Vect[2]
	return np.sqrt((X**2)+(Y**2)+(Z**2))
	
def Norm(Vect):
	"""
	Calculate the norm of a vector and then returns a normalized version of it
	"""
	X = Vect[0]
	Y = Vect[1]
	Z = Vect[2]
	NormValTemp = np.sqrt((X**2)+(Y**2)+(Z**2))
	return [X/NormValTemp,Y/NormValTemp,Z/NormValTemp]

	
#####################
#		Main		#
#####################	
	
if __name__=='__main__':
	if len(TblParts) != NbParts:
		print "Error, The number of part specified doesn't match the number of part names given !"
	else:
		for i in TblParts:
			Part = Mesh(i)	
			Part.ExportLocalCsys()
			#print i+" -  Done"
			stringmessage = "Part - "+str(i)+" - done, continue ?"
			print stringmessage
			reply = getWarningReply(message=stringmessage, buttons=(YES,NO))
			if reply == NO:	
				break
