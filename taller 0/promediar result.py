if __name__ == "__main__":
	fileReadParalelo=open("times-Paralelo.txt", "r").read()
	data=fileReadParalelo.split("\n")
	temp=0
	count=1
	result=0
	resultList=[]
	for i in data:
		tam, resultInput=i.split("	")
		if(temp!=tam):
			if(temp!=0):
				print("append")
				resultList.append(str(int(tam)-1)+"	"+str(result/count)+"\n")
			result=float(resultInput)
			temp=tam
			count=1
		else:
			result+=float(resultInput)
			count+=1
	fileWriteParalelo=open("times-resultParalelo.txt", "w")
	for i in resultList:
		fileWriteParalelo.write(i.replace(".", ","))
	fileWriteParalelo.close()
			






	#fileReadSecuencial=open("input.txt", "r").read()