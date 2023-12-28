#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import json

def separete(path,filename,pathw,filenamew):
	#items =[it1[],it2[],it3[],it4[],it5[],it6[],it7[],it8[],it9[],it10[],it11[],it12[],\
	#it13[],it14[],it15[],it16[],it17[],it18[],it18[],it19[],it20[],it21[],it22[],itX[],itY[]]
	items = []
	for i in range(25): 
		abc=[]
		items.append(abc)

	batch = []
	for i in range(25): 
		batch.append(0)

	for line in open(os.path.join(path,filename)).readlines():	
		split_list = line.split()
		temp=split_list[0]
		#print temp
		if "chrX" == temp:
			items[22].append(line)
			batch[22]= batch[22]+1
			if batch[22]==1000:
				file=open(os.path.join(pathw,filenamew[22]), "a" )  
				file.writelines(items[22])
				#print line 
				file.close() 
				batch[22]=0
				items[22][:]=[]
		elif temp == "chrY":
			items[23].append(line)
			batch[23] +=1
			if batch[23]==1000:
				file=open(os.path.join(pathw,filenamew[23]), "a" )  
				file.writelines(items[23])
					#print items 
				file.close() 
				batch[23]=0
				items[23][:]=[]

		elif temp.find('chr') == 0:
			try:
				indextemp=int(temp.split('r')[1])-1
				items[indextemp].append(line)
				batch[indextemp] +=1
				if batch[indextemp]==1000:
					file=open(os.path.join(pathw,filenamew[indextemp]), "a" )  
					file.writelines(items[indextemp])
				#print items 
					file.close() 
					batch[indextemp]=0
					items[indextemp][:]=[]
			except:
				items[24].append(line)
				batch[24] +=1
				if batch[24]==1000:
					file=open(os.path.join(pathw,filenamew[24]), "a" )  
					file.writelines(items[24])
				#print items 
					file.close() 
					batch[24]=0
					items[24][:]=[]	
		else:
			items[24].append(line)
			batch[24] +=1
			if batch[24]==1000:
				file=open(os.path.join(pathw,filenamew[24]), "a" )  
				file.writelines(items[24])
				#print items 
				file.close() 
				batch[24]=0
				items[24][:]=[]
				

	i=0
	for i in range(25):
		file=open(os.path.join(pathw,filenamew[i]), "a" )  
		file.writelines(items[i])
		#print items 
		file.close() 
				
	


def main():
	path="data/WGBS/"
	filename1="ENCFF179VKR.bed"
	pathw="data/WGBS/SK-N-SH_Separete_GRCh38"
	filenamew=["chr1out.bed","chr2out.bed","chr3out.bed","chr4out.bed","chr5out.bed",\
	"chr6out.bed","chr7out.bed","chr8out.bed","chr9out.bed","chr10out.bed","chr11out.bed","chr12out.bed",\
	"chr13out.bed","chr14out.bed","chr15out.bed","chr16out.bed","chr17out.bed","chr18out.bed","chr19out.bed","chr20out.bed",\
	"chr21out.bed","chr22out.bed","chrXout.bed","chrYout.bed","unfit.bed"]
	filenamew1=filenamew
	for i in range(25):
		filenamew1[i]="ENCFF179VKR"+filenamew[i]
	print(filenamew1[0])
	separete(path,filename1,pathw,filenamew1)

	filename1="ENCFF940XWW.bed"
	filenamew=["chr1out.bed","chr2out.bed","chr3out.bed","chr4out.bed","chr5out.bed",\
	"chr6out.bed","chr7out.bed","chr8out.bed","chr9out.bed","chr10out.bed","chr11out.bed","chr12out.bed",\
	"chr13out.bed","chr14out.bed","chr15out.bed","chr16out.bed","chr17out.bed","chr18out.bed","chr19out.bed","chr20out.bed",\
	"chr21out.bed","chr22out.bed","chrXout.bed","chrYout.bed","unfit.bed"]
	filenamew1=filenamew
	for i in range(25):
		filenamew1[i]="ENCFF940XWW"+filenamew[i]
	print(filenamew1[0])
	separete(path,filename1,pathw,filenamew1)

	path="data/WGBS/"
	filename1="ENCFF005TID.bed"	
	pathw="data/WGBS/A549_Separete_GRCh38"
	filenamew=["chr1out.bed","chr2out.bed","chr3out.bed","chr4out.bed","chr5out.bed",\
	"chr6out.bed","chr7out.bed","chr8out.bed","chr9out.bed","chr10out.bed","chr11out.bed","chr12out.bed",\
	"chr13out.bed","chr14out.bed","chr15out.bed","chr16out.bed","chr17out.bed","chr18out.bed","chr19out.bed","chr20out.bed",\
	"chr21out.bed","chr22out.bed","chrXout.bed","chrYout.bed","unfit.bed"]
	filenamew1=filenamew
	for i in range(25):
		filenamew1[i]="ENCFF005TID"+filenamew[i]
	print(filenamew1[0])
	separete(path,filename1,pathw,filenamew1)

	filename1="ENCFF003JVR.bed"
	filenamew=["chr1out.bed","chr2out.bed","chr3out.bed","chr4out.bed","chr5out.bed",\
	"chr6out.bed","chr7out.bed","chr8out.bed","chr9out.bed","chr10out.bed","chr11out.bed","chr12out.bed",\
	"chr13out.bed","chr14out.bed","chr15out.bed","chr16out.bed","chr17out.bed","chr18out.bed","chr19out.bed","chr20out.bed",\
	"chr21out.bed","chr22out.bed","chrXout.bed","chrYout.bed","unfit.bed"]
	filenamew1=filenamew
	for i in range(25):
		filenamew1[i]="ENCFF003JVR"+filenamew[i]
	print(filenamew1[0])
	separete(path,filename1,pathw,filenamew1)






if __name__ == '__main__':
	main()

	#use admin
	#db.shutdownServer()
