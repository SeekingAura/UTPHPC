import threading
import timeit
import random
import numpy as np

class matrices:
	def __init__(self, matrixA=[], matrixB=[]):
		self.matrixA=matrixA
		self.matrixB=matrixB
		self.matrixResult=[]
		for i in range(len(matrixA)):
			self.matrixResult.append([])			
			for j in range(len(matrixB)):
				self.matrixResult[i].append(0)
		self.matrixResult=np.array(self.matrixResult)
		self.posOp=0

	def operatePOS(self, pos):
		result=0
		filaPos=int(pos/len(self.matrixA))# find row to result
		columPos=pos%len(self.matrixA)# find column to result
		for row, i in enumerate(self.matrixA[filaPos]):
			# print("doing thread", pos)
			result+=i*self.matrixB[row][columPos]
		self.matrixResult[filaPos][columPos]=result

	def operateAll(self):
		for i in range(len(self.matrixA)):# per row
			
			for j in range(len(self.matrixB)): #per columns 
				result=0
				for enum, k in enumerate(self.matrixB):
					result+=self.matrixA[i][enum]*k[j]
				self.matrixResult[i][j]=result


if __name__ == "__main__":
	# matrix
	#print("execute paralelo")
	timeit.default_timer()
	fileRead=open("input.txt", "r").read()
	matrixA, matrixB=fileRead.split("\n")

	# Converts string to list
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
	#startTimeCreate=timeit.default_timer()
	hilos=[]
	for i in range(len(matrixA)*len(matrixA[0])):
		hilo1=threading.Thread(target=opMatrix.operatePOS, args=([i]))
		hilos.append(hilo1)
	#endTimeCreate=timeit.default_timer()
	#startTimeExecute=timeit.default_timer()
	for i in hilos:
		i.start()

	for i in hilos:
		i.join()
	#endTimeExecute=timeit.default_timer()
	endTime=timeit.default_timer()
	#print("tiempo creando hilos", #endTimeCreate-startTimeCreate)
	#print("tiempo execute ", endTimeExecute-startTimeExecute)
	#print("execute time", endTime-startTime)

	fileSave=open("times-Paralelo.txt", "a")
	fileSave.write(str(len(matrixA))+"	"+str(endTime-startTime).replace(".", ",")+"\n")
	fileSave.close()
	#print("END execute paralelo")

	#show result
	#print(opMatrix.matrixResult)