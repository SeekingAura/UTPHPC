if __name__ == "__main__":
	fileReadParalelo=open("timesc++Parallel.txt", "r").read()
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
		
	fileWriteParalelo=open("times-result.txt", "w")
	for i in resultDict:
		result=0
		for j in resultDict.get(i):
			result+=j
			# print("variable {}, total {}".format(i, j))
		fileWriteParalelo.write(str(i)+"	"+str(result/len(resultDict.get(i))).replace(".", ",")+"\n")
	fileWriteParalelo.close()
			






	#fileReadSecuencial=open("input.txt", "r").read()