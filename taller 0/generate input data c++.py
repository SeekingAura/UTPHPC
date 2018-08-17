import random
import sys
if __name__ == "__main__":
	# matrix
	#cant=int(input("indique el tamaño n de las matrices a crear -> "))
	if(len(sys.argv)!=2):#Verifica la cantidad de argumentos a la hora de compilar si no son 2. "py 'fichero.py' 'archivo'"
		sys.stderr.write('Usage: "{}" number_MaxN\n'.format(sys.argv[0]))#permite que al al compilar indique que debe de darse el archivo de la forma python.exe "fichero.py" "Archivo a abrir, como un simple print"
		raise SystemExit(1)#termina el programa
	negativeMAX=-100#-2147483647
	positiveMAX=100#2147483647
	cant=int(sys.argv[1])
	matrixA=[]
	for i in range(cant):
		matrixA.append([])			
		for j in range(cant):
			matrixA[i].append(random.randint(negativeMAX, positiveMAX))
	matrixB=[]
	for i in range(cant):
		matrixB.append([])			
		for j in range(cant):
			matrixB[i].append(random.randint(negativeMAX, positiveMAX))

	fileSave=open("inputc++.txt", "w")
	fileSave.write(str(len(matrixA))+"\n")
	fileSave.write(str(len(matrixA[0]))+"\n")
	for i in matrixA:
		fileSave.write(str(i).replace("[", "").replace("]", "").replace(" ", "")+"\n")
	fileSave.write(str(len(matrixB))+"\n")
	fileSave.write(str(len(matrixB[0]))+"\n")
	for i in matrixB:
		fileSave.write(str(i).replace("[", "").replace("]", "").replace(" ", "")+"\n")
	fileSave.close()