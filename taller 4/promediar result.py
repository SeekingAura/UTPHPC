import sys
if __name__ == "__main__":
	if(len(sys.argv)!=2):#Verifica la cantidad de argumentos a la hora de compilar si no son 2. "py 'fichero.py' 'archivo'"
		sys.stderr.write('Usage: "{}" "filename"\n'.format(sys.argv[0]))#permite que al al compilar indique que debe de darse el archivo de la forma python.exe "fichero.py" "Archivo a abrir, como un simple print"
		raise SystemExit(1)#termina el programa
	fileReadParalelo=open(sys.argv[1], "r").read()
	data=fileReadParalelo.split("\n")
	temp=0
	count=1
	result=0
	resultDict={}
	for i in data:
		try:
			tam, resultInput=i.split("	")
		except:
			# print("se hizo except {}".format(i))
			continue
		# print("value tam {} value result{}".format(tam, resultInput))
		if(tam in resultDict):
			resultDict.get(tam).append(float(resultInput))
		else:
			resultDict[tam]=[]
			resultDict.get(tam).append(float(resultInput))
		
	fileWriteParalelo=open(sys.argv[1].replace(".txt"," result.txt"), "w")
	times={}
	for i in [str(y) for y in sorted([int(x) for x in resultDict])]:# order by number
		result=0
		for j in resultDict.get(i):
			result+=j
			# print("variable {}, total {}".format(i, j))
		fileWriteParalelo.write(str(i)+"	"+str(result/len(resultDict.get(i))).replace(".", ",")+"\n")
		times[i]=len(resultDict.get(i))
	fileWriteParalelo.write("\n"+"#"*50+"\n")
	for time in times:
		fileWriteParalelo.write(str(time)+" -> "+str(times.get(time))+"\n")
	fileWriteParalelo.close()
			






	#fileReadSecuencial=open("input.txt", "r").read()