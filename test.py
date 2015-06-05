import cr_detection
import matplotlib.pyplot as plt
import numpy.random as rnd
import numpy
import time as tm
from scipy.integrate import simps
import os
import multiprocessing as mp


frame_t=27.75		#default 27.7 (long exposure
rdoutnoise=16/numpy.sqrt(8) #default 16/sqrt8 
gain=5
cr_class=cr_detection.cr_detection_class(rdoutnoise,gain,frame_t)
rnseed=213448622

rnd.seed(rnseed)

def makeramp(nframes,csmcmag,rmpslope,cr_class,rnseed,cosmdist='n',nlinpwr=0):
	numhits=0

	readN=numpy.zeros(nframes)
	cosmagarray=[0]*len(readN)


	#add a number of cosmics based on poisson distribution
	#size of cosmic taken from normal distribution
	#magnitude csmcmag and SD of 5. No cosmic in last frame. (because an increase
	#in last frame means cosmic was in the nframes-1st frame.
#	for j in range(nframes-1):
#		numhits=rnd.poisson(0.008)  #default 0.008 gives one cosmic per proper time period
#                if (numhits>0):
#                        readN[j]=j
#                        
#			if cosmdist=='n':
#				magadd=numpy.abs(rnd.normal(csmcmag,csmcmag/10.))
#			elif cosmdist=='f':
#				magadd=numpy.abs(rnd.uniform(1,csmcmag*2))
#                        cosmagarray[j]=magadd*numhits
#			numhits=0



	for j in range(1):
		wherehit=int(numpy.floor(numpy.random.uniform(0,nframes-1)))
		readN[wherehit]=wherehit
		cosmagarray[wherehit]=numpy.abs(rnd.normal(csmcmag,csmcmag/5.))
	

	time,cnts,cnterr = cr_class.mkramp(rmpslope,2400,nframes,cosmagarray,readN)
	if (nlinpwr)!=0:
		cnts=[c[i]*(float((i+3)/float(27.75))**(float(-nlinpwr))) for i in range(len(c))]
		max_cnts=max(cnts)
		c=[c[i]*(cnts[i]/max_cnts) for i in range(len(cnts))]

#	time,cnts,cnterr = cr_class.mkramp(rmpslope,1,nframes,cosmagarray,readN)
	
	return (time,cnts,cnterr,cosmagarray)


def detectCR(time,cnts,cnterr,method,threshold,cr_class):
	if (method=='2pt'):
		detected=cr_class.findCRs_2ptdiff(cnts,threshold)
	elif (method=='devfit'):
		detected=cr_class.findCRs_devfit(time,cnts,cnterr,threshold)
	elif (method=='yint'):
		detected=cr_class.findCRs_yint(time,cnts,cnterr,threshold)	
        elif (method=='qt1'):
                detected=cr_class.findCRs_qmeth(cnts,threshold)
	elif (method=='qt1m'):
                detected=cr_class.findCRs_qmeth_mod(cnts,threshold)
        elif (method=='qt2'):
                detected=cr_class.findCRs_qmeth2(cnts,threshold)
	elif (method=='GESD'):
                detected=cr_class.findCRs_GESD(cnts,threshold)
	elif (method=="IQR"):
		detected=cr_class.findCRs_IQR(time,cnts,threshold)
	else:
		print "Invalid detection method"
		return 0

	return detected


def runNramps(numramps,rmpslope,nframes,cosmag=0,method='2pt',thresh=1.0,dist='n'):
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
#	print 'creating ',numramps,' ramps'
	for rampn in range(0,numramps):
		t,c,e,m=makeramp(nframes,cosmag,rmpslope,cr_class,rnseed,dist)
		d=detectCR(t,c,e,method,thresh,cr_class)
		t_arr[0].append(t)
		c_arr[0].append(c)
		e_arr[0].append(e)
		m_arr[0].append(m)
		d_arr[0].append(d)
		#print rampn
	return t_arr[0],c_arr[0],e_arr[0],m_arr[0],d_arr[0]


def runPramps(numramps,rmpslope,nframes,cosmag=0,method='2pt',thresh=1.0):

#	output = mp.Queue()
#	processes = [mp.Process(target=runNramps, args=(numramps,rmpslope,nframes,cosmag,method,thresh)) for x in range(4)]
#
#	for p in processes:
#		p.start()
#		print "started", p
#	output.get(False)
#	for p in processes:
#		p.join()
#		print 'ended',p


#	results = output.get(False)
#	print(results)

	pool = mp.Pool()
#	results = [pool.apply(runNramps, args=(100,rmpslope,nframes,cosmag,method,thresh)) 

	for x in range(1,(numramps)):
		pool.apply(runNramps, args=(100,rmpslope,nframes,cosmag,method,thresh))


	pool.close()
	pool.join()
	#r=[results[i].get(timeout=525600*60) for i in range(1000-1)]
	print 'done main'	
	
	return results	





def missfalsedets(time,mags,dets):
	false=[]
	miss=[]
	truedets=[]
	truemiss=[]
	cosmics=[]
	print 'counting misses'
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
#				print('D'),
			else:
#				print('F'),
				falseramp+=1

		for val in mags[rampn]:
			k+=1
			if(val>0):
				if (((k) in dets[rampn])):
					pass
				else:
#					print('M'),
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


def rampsiterate(numramps,method='2pt',iterarray=['thresh',0.1,10,0.1],rmpslope=50,nframes=20,cosmag=20,thresh=3.0,dist='n'):
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
		
		t1,c1,e1,m1,d1=runNramps(int(numramps),float(rmpslope),int(nframes),float(cosmag),method,float(thresh),dist)
		f1,td1,M1,tm1,cn1=missfalsedets(t1,m1,d1)
		
		false[i]=numpy.sum(f1)
		truedets[i]=numpy.sum(td1)
		miss[i]=numpy.sum(M1)
		truemiss[i]=numpy.sum(tm1)
		dets[i]=numpy.sum(cn1)
		i+=1
	

	return false,truedets,miss,truemiss,dets


def makeROC(numramps,method='2pt',iterarray=['thresh',0.1,10,0.1],rmpslope=70,nframes=40,cosmag=10,thresh=3.0 ):
	f1,td1,M1,tm1,cm1=rampsiterate(numramps,method,iterarray,rmpslope,nframes,cosmag,thresh)
	truepos=[td1[i]/cm1[i] for i in range(len(cm1)-1)]
	falsepos=[f1[i]/((nframes-1)*numramps-cm1[i]) for i in range(len(cm1)-1)]
	x=numpy.arange(0,1.01,0.01)
	
	plt.scatter(falsepos,truepos,marker="+",c='k')
	plt.plot(x,x,linestyle='--',c='r')
	plt.title("ROC curve for "+method+" detection method varying the "+iterarray[0]+" crmag=10")
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


def plot_falsemiss(numramps,method='2pt',iterarray=['thresh',0.1,10,0.1],rmpslope=20,nframes=20,cosmag=50,thresh=3.0,dist='n'):
	f1,td1,M1,tm1,cm1=rampsiterate(numramps,method,iterarray,rmpslope,nframes,cosmag,thresh,dist)

	filename="plots_"+str(numramps)+"ramps_"+method+"_rmpslope="+str(rmpslope)+"_nframes="+str(nframes)+"_cosmag="+str(cosmag)+"_thresh="+str(thresh)+"_VAR_"+iterarray[0]+str(iterarray[1])+"->"+str(iterarray[2])
	
	falsepos=[f1[i]/float((nframes-1)*numramps-cm1[i]) for i in range(len(cm1)-1)]
	miss=[M1[i]/float(cm1[i]) for i in range(len(cm1)-1)]

	fd,=plt.plot(numpy.arange(iterarray[1],iterarray[2],iterarray[3]),falsepos,c='r',label='false detections')
	ms,=plt.plot(numpy.arange(iterarray[1],iterarray[2],iterarray[3]),miss,c='k',label='missed cosmics')

	plt.legend([fd, ms],["False detections","Missed Cosmic Rays"],loc=4,frameon=False)
	plt.title("Missed cosmics and False detections for "+method+" method with "+str(numramps)+" runs per point")
	plt.xlabel(str(iterarray[0])+" varied with a resolution of "+str(iterarray[3]))
	plt.ylabel("Fraction of missed cosmics or False detections")	

	fig = plt.gcf()
	fig.set_size_inches(18.5,10.5)
#	fig.savefig(method+str(rmpslope)+"_False_miss_.pdf",dpi=100)
#	plt.show()
	plt.close()	
	return [falsepos,miss]


def iterate_falsemiss(numramps,method='2pt',iterarray=['thresh',0.1,10,0.1],iterarray2=['rmpslope',1,100,10],rmpslope=20,nframes=20,cosmag=50,thresh=3.0):
	
	
	
	return 0

def plot_rampfalsemiss(numramps, method='2pt', rmpslope=20,nframes=20,cosmag=50,thresh=2.0,dist='n'):
        false=numpy.zeros(nframes)
        miss=numpy.zeros(nframes)
	cosmics=numpy.zeros(nframes)
	det=numpy.zeros(nframes)
        t1,c1,e1,mags,dets=runNramps(int(numramps),float(rmpslope),int(nframes),float(cosmag),method,float(thresh),dist) 
	for rampm in range(len(t1)):
	        for rampn in range(0,len(t1[rampm])):
	                k=0
	                for val in dets[rampm]:
	                        if(mags[rampm][val-1]>0):
					det[val]+=1
#	                                print('D'),
	                        else:
#	                                print('F'),
	                                false[val]+=1

        	        for val in mags[rampm]:
                	        k+=1
                        	if(val>0):
					cosmics[k]+=1
					#print " EEEEEEEEE", k, dets[rampm]
	                                if (((k) in dets[rampm])):
	                                        pass
	                                else:
#	                                        print('M'),
	                                        miss[k]+=1
	                        else:
	                                if(not((k+1) in dets[rampm])):
	                                        pass
	                                        #print('X'),


	
	
        import matplotlib.pyplot as plt
        xred=[i-0.2 for i in range(0,len(false))]
        xdet=[i-0.4 for i in range(0,len(det))]
        #xcosm=[i-0.6 for i in range(0,len(cosmics))]
        xmiss=[i-0.8 for i in range(0,len(miss))]
        miss=[miss[i]/cosmics[i] for i in range(len(miss))]
        false=[false[i]/cosmics[i] for i in range(len(false))]
        det=[det[i]/cosmics[i] for i in range(len(det))]
        mb=plt.bar(xmiss,miss,0.3,color='black')
        fb=plt.bar(xred,false,0.2,color='red')
        db=plt.bar(xdet,det,0.2,color='blue')
        #cb=plt.bar(xcosm,cosmics,0.2,color='green')
        plt.legend([mb[0],fb[0],db[0]],["missed","false","dets"])
        plt.title("number of detected,missed,false and cosmics versus ramp position")
        plt.xlabel("Frame number")
        plt.show()
        return false,miss,det,cosmics
	
	
	
def ROC_track_1(numramps,method='2pt',iterarray=['thresh',0.1,10,0.1],iterarray2=['cosmag',1,250,1],rmpslope=70,nframes=40):
	tp=0
	fp=0
	truepos=[]
	falsepos=[]
	thresharray=[]
	miss=[]
	thresh=4

	for cosmag in range(1,250,5):
		for thresh in numpy.arange(4,-0.0001,-0.2):
			
			t1,c1,e1,m1,d1=runNramps(int(numramps),float(rmpslope),int(nframes),float(cosmag),method,float(thresh))
			f1,td1,M1,tm1,cn1=missfalsedets(t1,m1,d1)
	        	tp=numpy.sum(td1)/float(numpy.sum(cn1))
	        	fp=numpy.sum(f1)/float(((nframes-1)*numramps-numpy.sum(cn1)))
			ms=numpy.sum(M1)/float(numpy.sum(cn1))
			print '\n', 'fp= ',fp," cosmag= ",cosmag
			if (fp>0.05):
				break

		truepos.append(tp)
		falsepos.append(fp)
		thresharray.append(thresh)
		miss.append(ms)

        x=numpy.arange(0,1.01,0.01)

        plt.scatter(falsepos,truepos,marker="+",c='k')
        plt.plot(x,x,linestyle='--',c='r')
        plt.title("ROC track 1 for "+method+" detection method and with cr_mag 0-250")
        plt.ylabel("true positive rate (detections)/(cosmics)")
        plt.xlabel('false positive rate (false detections)/(total points-cosmics)')
        plt.grid(b=True,which='major',color='b',linestyle=':')
        plt.xlim(-0.01,1.01)
        plt.ylim(-0.01,1.01)
        
        plt.show()

			

	
	return truepos,falsepos,miss,thresharray	


def ROC_track_2(numramps,method='2pt',iterarray=['thresh',0.1,10,0.1],iterarray2=['cosmag',1,250,1],rmpslope=70,nframes=40):
        tp=0
        fp=0
        truepos=[]
        falsepos=[]
        thresharray=[]
	miss=[]
        thresh=4


        for cosmag in range(1,250,25):

                for thresh in numpy.arange(4.0,-0.00001,-0.1):

                        t1,c1,e1,m1,d1=runNramps(int(numramps),float(rmpslope),int(nframes),float(cosmag),method,float(thresh))
                        f1,td1,M1,tm1,cn1=missfalsedets(t1,m1,d1)
                        tp=numpy.sum(td1)/float(numpy.sum(cn1))
                        fp=numpy.sum(f1)/float(((nframes-1)*numramps-numpy.sum(cn1)))
                        ms=numpy.sum(M1)/float(numpy.sum(cn1))
                        print '\n', 'tp,1-fp= ',tp,(1-fp)," cosmag= ",cosmag,thresh


			if (tp>(1-fp)):
                                break


                truepos.append(tp)
                falsepos.append(fp)
                thresharray.append(thresh)
	        miss.append(ms)

        x=numpy.arange(0,1.01,0.01)

        plt.scatter(falsepos,truepos,marker="+",c='k')
        plt.plot(x,x,linestyle='--',c='r')
        plt.title("ROC track2 for "+method+" detection method with cr_mag 0-250")
        plt.ylabel("true positive rate (detections)/(cosmics)")
        plt.xlabel('false positive rate (false detections)/(total points-cosmics)')
        plt.grid(b=True,which='major',color='b',linestyle=':')
        plt.xlim(-0.01,1.01)
        plt.ylim(-0.01,1.01)

        plt.show()




        return truepos,falsepos,miss,thresharray



