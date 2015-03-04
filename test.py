import cr_detection
import matplotlib.pyplot as plt
import numpy.random as rnd
import numpy
import time as tm
from scipy.integrate import simps
import os

frame_t=1
rdoutnoise=2
gain=1
cr_class=cr_detection.cr_detection_class(rdoutnoise,gain,frame_t)
rnseed=220022

rnd.seed(rnseed)

def makeramp(nframes,csmcmag,rmpslope,cr_class,rnseed,cosmdist='f'):
	numhits=0

	readN=numpy.zeros(nframes)
	cosmagarray=[0]*len(readN)


	#add a number of cosmics based on poisson distribution
	#size of cosmic taken from normal distribution
	#magnitude csmcmag and SD of 5. No cosmic in last frame. (because an increase
	#in last frame means cosmic was in the nframes-1st frame.
	for j in range(nframes-1):
		numhits=rnd.poisson(0.1)
                if (numhits>0):
                        readN[j]=j
                        
			if cosmdist=='n':
				magadd=numpy.abs(rnd.normal(csmcmag,5))
			elif cosmdist=='f':
				magadd=numpy.abs(rnd.uniform(1,csmcmag*2))
                        cosmagarray[j]=magadd*numhits
			numhits=0

	time,cnts,cnterr = cr_class.mkramp(rmpslope,1,nframes,cosmagarray,readN)
	
	return (time,cnts,cnterr,cosmagarray)


def detectCR(time,cnts,cnterr,method,threshold,cr_class):
	if (method=='2pt'):
		detected=cr_class.findCRs_2ptdiff(cnts,threshold)
	if (method=='devfit'):
		detected=cr_class.findCRs_devfit(time,cnts,cnterr,threshold)
	if (method=='yint'):
		detected=cr_class.findCRs_yint(time,cnts,cnterr,threshold)	
        if (method=='qt1'):
                detected=cr_class.findCRs_qmeth(cnts,threshold)

	return detected


def runNramps(numramps,rmpslope,nframes,cosmag=0,method='2pt',thresh=1.0):
	t_arr=[]
	c_arr=[]
	e_arr=[]
	m_arr=[]
	d_arr=[]
	t_arr.append([])
        c_arr.append([])
        e_arr.append([])
        m_arr.append([])
        d_arr.append([])

	for rampn in range(0,numramps):
		t,c,e,m=makeramp(nframes,cosmag,rmpslope,cr_class,rnseed)
		d=detectCR(t,c,e,method,thresh,cr_class)
		t_arr[0].append(t)
		c_arr[0].append(c)
		e_arr[0].append(e)
		m_arr[0].append(m)
		d_arr[0].append(d)

	return t_arr[0],c_arr[0],e_arr[0],m_arr[0],d_arr[0]


def missfalsedets(time,mags,dets):
	false=[]
	miss=[]
	truedets=[]
	truemiss=[]
	cosmics=[]

	for rampn in range(0,len(time)):
		falseramp=0
		missramp=0
		detramp=0
		Tmissramp=0
		cosmramp=sum(1 for i in mags[rampn] if i)
		k=0
		for val in dets[rampn]:
			if(mags[rampn][val-1]>0):
				detramp+=1
				print('D'),
			else:
				print('F'),
				falseramp+=1

		for val in mags[rampn]:
			k+=1
			if(val>0):
				if (((k) in dets[rampn])):
					pass
				else:
					print('M'),
					missramp+=1
			else:
				if(not((k+1) in dets[rampn])):
					Tmissramp+=1
					#print('X'),
			

		false.append(falseramp)
		miss.append(missramp)
		truedets.append(detramp)
		truemiss.append(Tmissramp)
		cosmics.append(cosmramp)

	return false,truedets,miss,truemiss,cosmics


def rampsiterate(numramps,method='2pt',iterarray=['thresh',0.1,10,0.1],rmpslope=50,nframes=20,cosmag=20,thresh=3.0):
        iterable=iterarray[0]
	if(method=='2pt'):
                print "iterating over "+iterable+" using 2-point difference method for a sample size of "+str(numramps)+" ramps per point"
	elif(method=='yint'):
                print "iterating over "+iterable+" using y-intercept method for a sample size of "+str(numramps)+" ramps per point"
	elif(method=='devfit'):
		print "iterating over "+iterable+" using deviation from fit method for a sample size of "+str(numramps)+" ramps per point"
        elif(method=='qt1'):
                print "iterating over "+iterable+" using First order Q test for a sample size of "+str(numramps)+" ramps per point"

	print " " 

	false=numpy.zeros(len(numpy.arange(iterarray[1],iterarray[2],iterarray[3]))+1)
	truedets=numpy.zeros(len(numpy.arange(iterarray[1],iterarray[2],iterarray[3]))+1)
	miss=numpy.zeros(len(numpy.arange(iterarray[1],iterarray[2],iterarray[3]))+1)
	truemiss=numpy.zeros(len(numpy.arange(iterarray[1],iterarray[2],iterarray[3]))+1)
	dets=numpy.zeros(len(numpy.arange(iterarray[1],iterarray[2],iterarray[3]))+1)
	
	i=0
	for iter1 in numpy.arange(iterarray[1],iterarray[2],iterarray[3]):
		exec(iterarray[0] + "= iter1")
		
		print iterarray[0]		
                print iter1
		
		t1,c1,e1,m1,d1=runNramps(int(numramps),float(rmpslope),int(nframes),float(cosmag),method,float(thresh))
		f1,td1,M1,tm1,cn1=missfalsedets(t1,m1,d1)
		
		false[i]=numpy.sum(f1)
		truedets[i]=numpy.sum(td1)
		miss[i]=numpy.sum(M1)
		truemiss[i]=numpy.sum(tm1)
		dets[i]=numpy.sum(cn1)
		i+=1
	

	return false,truedets,miss,truemiss,dets


def makeROC(numramps,method='2pt',iterarray=['thresh',0.1,10,0.1],rmpslope=20,nframes=20,cosmag=50,thresh=3.0 ):
	f1,td1,M1,tm1,cm1=rampsiterate(numramps,method,iterarray,rmpslope,nframes,cosmag,thresh)
	truepos=[td1[i]/cm1[i] for i in range(len(cm1)-1)]
	falsepos=[f1[i]/((nframes-1)*numramps-cm1[i]) for i in range(len(cm1)-1)]
	x=numpy.arange(0,1.01,0.01)
	
	plt.scatter(falsepos,truepos,marker="+",c='k')
	plt.plot(x,x,linestyle='--',c='r')
	plt.title("ROC curve for "+method+" detection method varying the "+iterarray[0])
	plt.ylabel("true positive rate (detections)/(cosmics)")
	plt.xlabel('false positive rate (false detections)/(total points-cosmics)')
	plt.grid(b=True,which='major',color='b',linestyle=':')
	plt.xlim(-0.01,1.01)
	plt.ylim(-0.01,1.01)

	plt.show()
	fp1=[falsepos[i] for i in range(len(truepos)-1) if falsepos[i]>0]
	tr1=[truepos[i] for i in range(len(truepos)-1) if falsepos[i]>0]

	ROCAUC=-simps(tr1,fp1)
	print "Area Under Curve= "+str(ROCAUC)
	return ROCAUC, [truepos, falsepos]


def plot_falsemiss(numramps,method='2pt',iterarray=['thresh',0.1,10,0.1],rmpslope=20,nframes=20,cosmag=50,thresh=3.0):
	f1,td1,M1,tm1,cm1=rampsiterate(numramps,method,iterarray,rmpslope,nframes,cosmag,thresh)

	filename="plots_"+str(numramps)+"ramps_"+method+"_rmpslope="+str(rmpslope)+"_nframes="+str(nframes)+"_cosmag="+str(cosmag)+"_thresh="+str(thresh)+"_VAR_"+iterarray[0]+str(iterarray[1])+"->"+str(iterarray[2])
	
	falsepos=[f1[i]/((nframes-1)*numramps-cm1[i]) for i in range(len(cm1)-1)]
	miss=[M1[i]/(cm1[i]) for i in range(len(cm1)-1)]

	fd,=plt.plot(numpy.arange(iterarray[1],iterarray[2],iterarray[3]),falsepos,c='r',label='false detections')
	ms,=plt.plot(numpy.arange(iterarray[1],iterarray[2],iterarray[3]),miss,c='k',label='missed cosmics')

	plt.legend([fd, ms],["False detections","Missed Cosmic Rays"],loc=4,frameon=False)
	plt.title("Missed cosmics and False detections for "+method+" method with "+str(numramps)+" runs per point")
	plt.xlabel(str(iterarray[0])+" varied with a resolution of "+str(iterarray[3]))
	plt.ylabel("Fraction of missed cosmics or False detections")	

	fig = plt.gcf()
	fig.set_size_inches(18.5,10.5)
	fig.savefig(method+"_False_miss_.pdf",dpi=100)
	plt.show()
	
	return [fd,ms]


def iterate_falsemiss(numramps,method='2pt',iterarray=['thresh',0.1,10,0.1],iterarray2=['rmpslope',1,100,1],rmpslope=20,nframes=20,cosmag=50,thresh=3.0):
	
	
	
	return 0








