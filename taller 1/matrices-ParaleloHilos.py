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
		self.matrixResult=[]
		self.matrixResult=range(self.matrixArows*self.matrixBcols)
		self.matrixResult=np.array(self.matrixResult)
		self.posOp=0

	def operatePOS(self, pos, lock):
		result=0
		filaPos=int(pos/self.matrixArows)# find row to result
		columPos=pos%self.matrixBcols# find column to result
		for i in range(self.matrixArows):
			result+=self.matrixA[filaPos*self.matrixArows+i]*self.matrixB[i*self.matrixBcols+columPos]
		lock.acquire()
		self.matrixResult[filaPos*self.matrixArows+columPos]=result
		lock.release()

	def operateRow(self, row, lock):
		result=0
		for i in range(self.matrixAcols):
			for j in range(self.matrixBcols):
				result+=self.matrixA[row*self.matrixArows+j]*self.matrixB[j*self.matrixBcols+i]
			lock.acquire()
			self.matrixResult[row*self.matrixArows+i]=result
			lock.release()
			result=0

	def operateAll(self):
		for i in range(self.matrixArows):# per row
			
			for j in range(self.matrixBcols): #per columns 
				result=0
				for k in range(self.matrixArows):
					result+=self.matrixA[i*self.matrixArows+k]*self.matrixB[k*self.matrixBcols+j]
				self.matrixResult[self.matrixArows*i+j]=result


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
	
	#instnaciando objeto
	opMatrix=matrices(matrixA, matrixB)
	
	startTime=timeit.default_timer()
	hilos=[]
	# for i in range(opMatrix.matrixArows*opMatrix.matrixBcols):
	# 	hilo1=threading.Thread(target=opMatrix.operartePOS, args=([i, threading.Lock()]))
	# 	hilos.append(hilo1)

	for i in range(opMatrix.matrixArows):
		hilo1=threading.Thread(target=opMatrix.operateRow, args=([i, threading.Lock()]))
		hilos.append(hilo1)

	for i in hilos:
		i.start()

	for i in hilos:
		i.join()
	endTime=timeit.default_timer()

	fileSave=open("times-Threads.txt", "a")
	fileSave.write(str(len(matrixA))+"	"+str(endTime-startTime).replace(".", ",")+"\n")
	fileSave.close()

	#show result
	#print(opMatrix.matrixResult)