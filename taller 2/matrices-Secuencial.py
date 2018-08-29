import threading
import timeit
import random
import numpy as np

class matrices:
	def __init__(self, matrixA=[], matrixB=[]):
		self.matrixA=matrixA
		self.matrixArows=len(matrixA)
		self.matrixAcols=len(matrixA[0])
		self.matrixA=self.matrixA.flatten()
		self.matrixB=matrixB
		self.matrixBrows=len(matrixB)
		self.matrixBcols=len(matrixB[0])
		self.matrixB=self.matrixB.flatten()
		self.matrixResult=range(self.matrixArows*self.matrixBcols)
		self.matrixResult=np.array(self.matrixResult)
		self.posOp=0

	def operarPOS(self, pos):
		result=0
		filaPos=int(pos/self.matrixArows)# find row to result
		columPos=pos%self.matrixBcols# find column to result
		for i in range(self.matrixArows):
			result+=self.matrixA[filaPos*self.matrixArows+i]*self.matrixB[i*self.matrixBcols+columPos]
		
		self.matrixResult[filaPos*self.matrixArows+columPos]=result

	def operateAll(self):
		for i in range(self.matrixArows):# per row
			
			for j in range(self.matrixBcols): #per columns 
				result=0
				for k in range(self.matrixArows):
					result+=self.matrixA[i*self.matrixArows+k]*self.matrixB[k*self.matrixBcols+j]
				self.matrixResult[self.matrixArows*i+j]=result
				

if __name__ == "__main__":
	# matrix
	#print("execute secuential")
	timeit.default_timer()
	
	fileRead=open("input.txt", "r").read()
	matrixA, matrixB=fileRead.split("\n")


	matrixA=eval(matrixA)
	matrixB=eval(matrixB)

	matrixA=np.array(matrixA, np.int32)
	matrixB=np.array(matrixB, np.int32)
	### show input data Matrix A
	#print("Matrix A")
	#for i in matrixA:
	#	print(i)
	#print("Matrix B")
	#for i in matrixB:
	#	print(i)
	### show input data Matrix B
	opMatrix=matrices(matrixA, matrixB)
	startTime=timeit.default_timer()
	opMatrix.operateAll()
	endTime=timeit.default_timer()
	#print("execute time", endTime-startTime)

	fileSave=open("times-Secuencial.txt", "a")
	fileSave.write(str(len(matrixA))+"	"+str(endTime-startTime).replace(".", ",")+"\n")
	fileSave.close()
	#print("END execute secuential")
	#show result
	#print("result")
	#print(opMatrix.matrixResult)