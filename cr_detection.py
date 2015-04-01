#! /usr/bin/env python

''' 
ABOUT: 
 This script defines a class used to make ramps with/without cosmic rays for IR detectors,
 and includes algorithms to detect cosmic rays using methods described in
 Anderson, R., & Gordon, K. 2011, PASP, 132, 1237

DEPENDS: 
Python 2.7 

AUTHOR: 
Rachel E. Anderson, 2012

''' 
from math import sqrt
import random
import numpy as np
import operator

class cr_detection_class():

    def __init__(self, readnoise, gain, frame_t, frames_per_group=1):

        self.frame_t = frame_t  #seconds
        self.gain = gain  #e- per DN
        self.readnoise = readnoise  #DN per read

    def findCRs_qmeth2(self, counts, rej_thr=0.2):
        cr_readNs = []
        # diff takes 1st-order differences between the elements of counts
        diff = np.diff(counts) 
        rn2 = self.readnoise * self.readnoise
        max_crs = len(diff)   # Upper limit, there won't be this many though...
        cr_free=np.copy(counts)
        cr_loc=[]
        
        # See if there is a CR in all of these segments
        for j in xrange(max_crs):
            flag = 0
            # Remove outliers we already found from diff by setting them to close to the minimum value of diff
            diff = np.diff(cr_free)
       	    diff2= [cr_free[i]-cr_free[i-2] for i in range(2,len(cr_free)-1,1)]
#	    diff2= diff2[::-1]
#   	    print len(diff2) 
            if (len(cr_readNs) != 0):
                cr_loc=[cr_readNs[i]-1 for i in range(len(cr_readNs))]
                for i in range(len(cr_loc)):
#			print cr_loc[i]
                        diff[cr_loc[i]]=min(diff) #replaces 'cr' with minimum
           		diff2[cr_loc[i]-2]=min(diff2) #replaces 'cr' with minimum
 
            Qrange=(max(diff)-min(diff))
	    Qrange2=(max(diff2)-min(diff2))
            
            if (Qrange<=0):
                break
            
            Q1=np.diff(np.sort(diff))/Qrange
	    Q2=np.diff(np.sort(diff2))/Qrange2
#	    print Q2
            diffarg=np.argsort(diff)	
	    diffarg2=np.argsort(diff2)
            
#           print Q1
#           print diffarg
            # Calculate ratio to compare to rejection threshold
            max_index, max_value = max(enumerate(np.argsort(diff)), key=operator.itemgetter(1))          
	    max_index2, max_value2 = max(enumerate(np.argsort(diff2)), key=operator.itemgetter(1))            
            if ((max(Q1) >  rej_thr) and (max(Q2)> (rej_thr/2.))): # use normal rejection threshold to rejecm
                    if((diffarg2[max_value2]+1) in cr_readNs): break #sometimes triggers on already found CR     
		    print(diffarg2[max_value2]+1)	
                    cr_readNs.append(diffarg2[max_value2]+1)   # This is the frame the CR first appears in         
                    flag = 1
            else: break 
            if (flag == 0): break   # We didn't find any outliers, so we are done
	
	return cr_readNs


    def findCRs_qmeth(self, counts, rej_thr=0.2):
        cr_readNs = []
        # diff takes 1st-order differences between the elements of counts
        diff = np.diff(counts)
        rn2 = self.readnoise * self.readnoise
        max_crs = len(diff)   # Upper limit, there won't be this many though...
        cr_free=np.copy(counts)
        cr_loc=[]

        # See if there is a CR in all of these segments
        for j in xrange(max_crs):  
            flag = 0
            # Remove outliers we already found from diff by setting them to close to the minimum value of diff
            diff = np.diff(cr_free)
	    
            if (len(cr_readNs) != 0):
                cr_loc=[cr_readNs[i]-1 for i in range(len(cr_readNs))]
                for i in range(len(cr_loc)):
                        diff[cr_loc[i]]=min(diff) #replaces 'cr' with minimum

            Qrange=(max(diff)-min(diff))
#Qrange=[max(numpy.delete(diff,i))-min(numpy.delete(diff,i))  for i in range(len(diff))]  sdiff=np.sort(diff)   Q1d=[(sdiff[i+1]-sdiff[i])/Qrange[i] for i in range(len(sdiff)-1)]

            if (Qrange<=0):
                break 

            Q1=np.diff(np.sort(diff))/Qrange
	    diffarg=np.argsort(diff)
	    
#	    print Q1
#	    print diffarg
            # Calculate ratio to compare to rejection threshold
            max_index, max_value = max(enumerate(np.argsort(diff)), key=operator.itemgetter(1))
	    
            if (Q1[-1]) >  rej_thr: # use normal rejection threshold to rejecm
                    if((diffarg[max_value]+1) in cr_readNs): break #sometimes triggers on already found CR
		    cr_readNs.append(diffarg[max_value]+1)   # This is the frame the CR first appears in
                    flag = 1
            else: break
            if (flag == 0): break   # We didn't find any outliers, so we are done



        return cr_readNs




    def findCRs_2ptdiff(self, counts, rej_thr=3.0): 
        cr_readNs = []
        # diff takes 1st-order differences between the elements of counts
        diff = np.diff(counts)   
        rn2 = self.readnoise * self.readnoise
        max_crs = len(diff)   # Upper limit, there won't be this many though...
    
        # See if there is a CR in all of these segments
        for j in xrange(max_crs): 
            flag = 0
            cr_free = np.copy(diff)
        
            # Remove outliers we already found from cr_free
            if (len(cr_readNs) != 0):
                cr_free = np.delete(cr_free,cr_readNs)

            # Calculate ratio to compare to rejection threshold
            meddiff = np.median(cr_free)
            pnoise = sqrt(abs(meddiff) * self.gain) / self.gain
            diff_unc = sqrt(pnoise*pnoise + rn2 + rn2)   # We have two read's            
            ratio = np.abs(diff - meddiff) / diff_unc
            sortindx = np.argsort(ratio)[::-1]  # Indices that would reverse-sort an array

            # Remove outliers
            for i in sortindx:  # Largest first
                if ((i+1) in cr_readNs): continue   # Skip those already tagged as CRs
                if ratio[i] > rej_thr:  
                    cr_readNs.append(i+1)   # This is the frame the CR first appears in
                    flag = 1
                break   # If we find one break, if the largest isn't one, break.
            if (flag == 0): break   # We didn't find any outliers, so we are done
        return cr_readNs
#	return readN

    def devfitramp(self, time_in, counts_in, count_errors_in, rej_thr=3.0):
#        import ipdb
#        ipdb.set_trace()
	if(len(time_in)<3.):
		return 0
        fits,fitserr,fity,fityerr=self.fitline(time_in,counts_in, random_only=0)
	pnoise = np.array([sqrt(abs(fits)*self.frame_t*self.gain)/self.gain for k in xrange(len(counts_in))])

	count_errors_in=[np.sqrt(pnoise[i]*pnoise[i]+ self.readnoise*self.readnoise) for i in xrange(len(counts_in))]
	dev=[(counts_in[i]-(fits*time_in[i]+fity))/(count_errors_in[i]) for i in range(len(time_in))]
#        print dev
	diffdev=[np.absolute(dev[i]-dev[i-1]) for i in range(1,len(dev))]
#	print(diffdev)
        (mx,ind)=max((v,j) for j,v in enumerate(diffdev))
	if (mx>rej_thr):
		return (ind+1)
	else:
		return 0


    def findCRs_yint(self,time_in,counts_in,count_errors_in,rej_thr=3.0):
        crcnt = 0
        goodflags = np.empty(0)
        cr_readNs = np.array([len(counts_in)-1])  # For the sake of dealing with the last semi-ramp
                                                  #I add this as a marker (I will remove later)
        for iteration in xrange(len(counts_in)-1):   # For each possible semi-ramp
            sortedreadNs = np.sort(cr_readNs)  # I could have added a readN
	                                       #at the end which belongs between.
            for i in xrange(len(sortedreadNs)):  # I would loop over cr_readNs, but
                                                 #it is a growing array, instead
                                                 #I just have to break this once
                                                 #we loop over all semi-ramps
                if sortedreadNs[i] in goodflags: 
			#print(goodflags)	
			continue  # This semi-ramp is fine

		#Pick out CR free ramp segment
                start = end = None
                if (crcnt != 0):
                    if (i == 0):  # If this is the first CR..
                        start = 0
                    else:
                        start = sortedreadNs[i-1]+1 
                    end = sortedreadNs[i]+1     
                else:
                    start = 0
                    end = len(counts_in) 
		time = time_in[start:end]
                counts = counts_in[start:end]
                count_errors = count_errors_in[start:end]

                # If ramp is too short, maybe start is a cosmic not detected! 
		#(low threshold won't catch this CR if we 'continue')
		# ALSO, CRs don't seem to be detected at begin/end of semiramps!
                if (end - start) < 3:
		   b=0
                   goodflags = np.append(goodflags,sortedreadNs[i])
		   #print "Gii= ",i,'startend=',start,end,'g',goodflags
		   if(start in sortedreadNs):continue
		   if(start == 0): continue
		   t2pread=self.findCRs_2ptdiff(counts_in[(start-1):end],rej_thr)
		   #print (t2pread+start)
		   if (len(t2pread)==2):
		   	if not((start+t2pread[0]-1) in sortedreadNs):
				cr_readNs = np.append(cr_readNs,start + t2pread[0]-1)
		                crcnt += 1
				b=1
			if (not((start +t2pread[1]-1) in sortedreadNs) and not(t2pread[1]==len(time_in)) ):
                        	cr_readNs = np.append(cr_readNs,start + t2pread[1]-1)
	                        crcnt += 1
				b=1
                        if b : break
		   continue

                # Now search for CRs
                n_possible_crs = len(counts)-1
                slp = np.empty((2,n_possible_crs))
                slp_err = np.empty((2,n_possible_crs))
                yint = np.empty((2,n_possible_crs))
                yint_err = np.empty((2,n_possible_crs))
                yint_err_sqr = np.empty((2,n_possible_crs))
                for rn in xrange(len(counts)-1):  # Iterate over possible cr_readNs
                                                  #(forget first two and last three,
                                                  #as they must be in line, and remember
                                                  #that last is not included)
                    if rn == 0:  # Shift the x-axis to be '0' at single sample
                        slp[0][rn],slp_err[0][rn],yint[0][rn],yint_err[0][rn] = 0,99999,counts[rn],self.readnoise  
                        slp[1][rn],slp_err[1][rn],yint[1][rn],yint_err[1][rn] = self.fitline(time[rn+1:]-time[rn],counts[rn+1:],random_only=1)
                    elif rn == len(counts)-2: # Again, shift the x-axis to be '0 at single sample
                        slp[0][rn],slp_err[0][rn],yint[0][rn],yint_err[0][rn] = self.fitline(time[:rn+1]-time[rn+1],counts[:rn+1],random_only=1)  
                        slp[1][rn],slp_err[1][rn],yint[1][rn],yint_err[1][rn] = 0,99999,counts[rn+1],self.readnoise   
                    else:
                        slp[0][rn],slp_err[0][rn],yint[0][rn],yint_err[0][rn] = self.fitline(time[:rn+1]-time[rn],counts[:rn+1],random_only=1) 
                        slp[1][rn],slp_err[1][rn],yint[1][rn],yint_err[1][rn] = self.fitline(time[rn+1:]-time[rn],counts[rn+1:],random_only=1)

                # For each possible CR, calculate the average slope for the Poisson noise (photon noise)
                slpavg = np.empty(n_possible_crs)
                for k in xrange(n_possible_crs):
                    if slp_err[0][k] == 99999:
                        slpavg[k] = slp[1][k]
                    elif slp_err[1][k] == 99999:
                        slpavg[k] = slp[0][k]
                    else:
                        slpavg[k] = (slp[0][k]/(slp_err[0][k]*slp_err[0][k]) + slp[1][k]/(slp_err[1][k]*slp_err[1][k]))/(1/(slp_err[0][k]*slp_err[0][k]) + 1/(slp_err[1][k]*slp_err[1][k]))
                pnoise = np.array([sqrt(abs(slpavg[k])*self.frame_t*self.gain)/self.gain for k in xrange(len(slpavg))])

                #Add photon noise and readnoise in quadrature to get expected uncertainty in the y-intercept
                yerr_exp = np.array([sqrt(pnoise[k]*pnoise[k] + yint_err[0][k]*yint_err[0][k] + yint_err[1][k]*yint_err[1][k]) for k in xrange(len(pnoise))])   #we have two read's
            
                ydiff = np.abs(yint[1] - yint[0])  #random_only=1 doesn't mess with how yints are calculated.
                ratio = ydiff/yerr_exp
                rn2add = ratio.argmax()
                if ratio[rn2add] > rej_thr:
                    cr_readNs = np.append(cr_readNs,start + rn2add)  
                    crcnt += 1
		    #print "i=",i,'startend=',start,end,'cr',start+rn2add
                    break   # Every time a cr is found, go on to next iteration
                else:
                    goodflags = np.append(goodflags,sortedreadNs[i])  #the semi-ramp ending at this read number is good to go!
                    if len(goodflags) == crcnt+1: break  

        #remember, we put the last frame in cr_readNs so we need to remove it
        to_rm = np.where(cr_readNs == len(counts_in)-1)[0][0]
        if cr_readNs[to_rm] != len(counts_in)-1:  #it does not equal the last read
            print 'The last item in list cr_readNs:',cr_readNs[to_rm],' is not the final read: ',len(counts)-1,'.  Do not remove it!'
            sys.exit()
        cr_readNs = np.delete(cr_readNs,to_rm)  
        return cr_readNs + 1  # Jump is detected between frame before CR and frame after CR


    def findCRs_devfit(self,time_in,counts_in,count_errors_in,rej_thr=3.0):
        crcnt = 0
        goodflags = np.empty(0)
        cr_readNs = np.array([len(counts_in)-1])  # For the sake of dealing with the last semi-ramp
                                                  #I add this as a marker (I will remove later)
        for iteration in xrange(len(counts_in)-1):   # For each possible semi-ramp
            sortedreadNs = np.sort(cr_readNs)  # I could have added a readN
	                                       #at the end which belongs between.
            for i in xrange(len(sortedreadNs)):  # I would loop over cr_readNs, but
                                                 #it is a growing array, instead
                                                 #I just have to break this once
                                                 #we loop over all semi-ramps
                if sortedreadNs[i] in goodflags: 
# 			print(goodflags)	
			continue  # This semi-ramp is fine

		#Pick out CR free ramp segment
                start = end = None
                if (crcnt != 0):
                    if (i == 0):  # If this is the first CR..
                        start = 0
                    else:
                        start = sortedreadNs[i-1]+1 
                    end = sortedreadNs[i]+1     
                else:
                    start = 0
                    end = len(counts_in) 
		time = time_in[start:end]
                counts = counts_in[start:end]
                count_errors = count_errors_in[start:end]
		
                # If ramp is too short, maybe start is a cosmic not detected! 
		#(low threshold won't catch this CR if we 'continue')
		# ALSO, CRs don't seem to be detected at begin/end of semiramps!
                if (end - start) < 3:
		   b=0
                   goodflags = np.append(goodflags,sortedreadNs[i])
#		   print "Gii= ",i,'startend=',start,end,'g',goodflags
		   if(start in sortedreadNs):continue
		   if(start == 0): continue
		   t2pread=self.findCRs_2ptdiff(counts_in[(start-1):end],rej_thr)
# 	 	   print (t2pread+start)
		   if (len(t2pread)==2):
		   	if not((start + t2pread[0]-1) in sortedreadNs):
				cr_readNs = np.append(cr_readNs,start + t2pread[0]-1)
		                crcnt += 1
				b=1
			if not((start + t2pread[1]-1) in sortedreadNs):
                        	cr_readNs = np.append(cr_readNs,start + t2pread[1]-1)
	                        crcnt += 1
				b=1
                        if b : break

		   continue

                # Now search for CRs

#	print 'start',start
#		pnoise = np.array([sqrt(abs(slpavg[k])*self.frame_t*self.gain)/self.gain for k in xrange(len(slpavg))])

                #Add photon noise and readnoise in quadrature to get expected uncertainty in the y-intercept
#                err_exp = np.array([sqrt(pnoise[k]*pnoise[k] + yint_err[0][k]*yint_err[0][k] + yint_err[1][k]*yint_err[1][k]) for k in xrange(len(pnoise))])   #we have two read's
                cosmindex=self.devfitramp(time, counts,count_errors_in,rej_thr)
                if cosmindex:
                    cr_readNs = np.append(cr_readNs, start + cosmindex-1)  
                    crcnt += 1
	            #print cr_readNs
		    #print "i=",i,'startend=',start,end,'cr',start+rn2add
                    break   # Every time a cr is found, go on to next iteration
                else:
                    goodflags = np.append(goodflags,sortedreadNs[i])  #the semi-ramp ending at this read number is good to go!
                    if len(goodflags) == crcnt+1: break  

        #remember, we put the last frame in cr_readNs so we need to remove it
        to_rm = np.where(cr_readNs == len(counts_in)-1)[0][0]
        if cr_readNs[to_rm] != len(counts_in)-1:  #it does not equal the last read
            print 'The last item in list cr_readNs:',cr_readNs[to_rm],' is not the final read: ',len(counts)-1,'.  Do not remove it!'
            sys.exit()
        cr_readNs = np.delete(cr_readNs,to_rm)  
        return cr_readNs + 1  # Jump is detected between frame before CR and frame after CR


    # Returns slope, slope error, y-intercept, y-intercept error.
    # Calculates the fit using a covariance matrix taking into account
    # correlated and uncorrelated errors.
    # Note that random_only == 1 if it is yint method and we want only random error
    # for the y-intercept error
    def fitline(self, xx, yy, random_only=0):
        #initial estimate of the slope not considering correlated and uncorrelated errors
        slope,yint = np.polyfit(xx,yy,1) 
        photon_noise = sqrt(abs(slope) * self.frame_t * self.gain) / self.gain  
        nn = len(xx)
        YY = np.asmatrix(yy).T
        AA = np.asmatrix([np.ones(nn),xx]).T
        guts = [np.zeros(nn)] * nn
        PP = np.asmatrix(guts)

        #photon noise (which is correlated)
        pn2 = photon_noise * photon_noise   
        pn0 = pn2   #this is the noise we will count in the first read.  
        for j in xrange(nn):
            for i in xrange(j,nn): 
                PP[i,j] = PP[j,i] = pn2 * (j + 1)   
                
        #now add readnoise on the diagonal
        RR = np.asmatrix(np.identity(nn) * self.readnoise * self.readnoise)  
        CC = PP + RR

        #calculate slope and y-intercept
        part1 = AA.T * (CC.I * AA)
        part2 = AA.T * (CC.I * YY)
        XX = part1.I * part2
        bb = XX[0,0]
        mm = XX[1,0]

        #calculate slope and y-intercept error
        mm_err = sqrt(part1.I[1,1])
        if random_only == 1:
            part1 = AA.T * (RR.I * AA)   #random error only
        bb_err = sqrt(part1.I[0,0])
        #chi2 = (YY - AA * XX).T * CC.I * (YY - AA * XX)
        return mm,mm_err,bb,bb_err
 

    # readN is the first read the CR is detected in, so the first read AFTER the CR hit.
    # Can make a CR-free ramp if you don't include crmag and readN
    # crmag is in DN.  crmag and readN's should be arrays
    def mkramp(self, slope, yint, nframes, crmag=None, readN=None):  
        # Calculate time
        int_t = (nframes-1)*self.frame_t  
        time = np.linspace(self.frame_t, int_t + self.frame_t, nframes)

        # Get counts
        cnts2add = slope*self.frame_t
        cnts2add_electrons = cnts2add * self.gain
        counts = np.empty(nframes)   #not actually empty, just has whatever garbage this memory contains.
        counts[0] = yint + np.random.poisson(lam=cnts2add_electrons) / self.gain  
        for frame in xrange(1,nframes):
            #included correlated noise
            counts[frame] = counts[frame-1] +  np.random.poisson(lam=cnts2add_electrons) / self.gain
            
            # Add CR if necessary
            if crmag and ((frame-1) in readN):   # Occured between frame n and n+1
                counts[frame] += crmag[np.where(readN == (frame-1))[0][0]] 
                
        #now add un-correlated noise
        counts += [random.normalvariate(0,self.readnoise) for frame in xrange(nframes)]

        #Calculate errors
        photonNoise = np.array([sqrt(abs(counts[i])*self.gain)/self.gain for i in xrange(nframes)])
	#photonNoise =np.array([sqrt(abs(counts[i])*0)/self.gain for i in xrange(nframes)])
        rn2 = self.readnoise*self.readnoise
        count_errors = np.array([sqrt(photonNoise[i]*photonNoise[i] + rn2) for i in xrange(nframes)])
        return time, counts, count_errors
    

