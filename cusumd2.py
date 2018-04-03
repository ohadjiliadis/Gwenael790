from numpy import *
from scipy.io import write_array
from scipy.io import read_array

def cusumd( col1, d, start, end, angle=45, thresh=100 ):
    detectO=-1
    rd=[0]
    for n in range(start,end):
        if (col1[n]!=-10):
            val=max(0,rd[size(rd)-1]-(col1[n]-angle))
            d[n]=val
            rd.append(val)
        else:
            d[n]=-10
            
        if d[n]>thresh: 
            detectO=n
            break;
    return detectO 




def cusume( col1, e,  start, end, angle=45, thresh=100 ):
    detectO=-1
    rd=[0]
    for n in range(start,end):
        if (col1[n]!=-10):
            val=max(0,rd[size(rd)-1]+(col1[n]-angle))
            e[n]=val
            rd.append(val)
        else:
            e[n]=-10
            
        if e[n]>thresh: 
            detectO=n
            break;
    return detectO

#
    
#HELPER, DO NOT WORRY ABOUT IT NOW
def invalidate_missing(col1,object, start, end):
    for n in range(start,end):
        if (col1[n]==-10):
            object[n]=-10
            
            
    
#detect sequence low to high to low to high ...
#need to start from low
def cusumall( col1, d, e, object, angled=45, threshd=20, anglee=30, threshe=15):
    running=True
    start=1
    end=size(col1)
    while running:
        startObject=cusumd( col1, d,  start, end, angled, threshd )
        #print "startObject: ",startObject
        if (startObject==-1):
            running=False
        else:
            endObject=cusume( col1, e,  startObject+1, end, anglee, threshe )
            #print "endObject: ",endObject
            if (endObject==-1):
                object[startObject:end]=10
                #invalidate_missing(col1,object,startObject,end);
                running=False
            else:
                object[startObject:endObject]=10
                #invalidate_missing(col1,object,startObject,endObject);
                start=endObject
                
    invalidate_missing(col1, object, 0, end)
    

#Performs cusumall for all columns
#NO NEED FOR IT NOW
#data is a set of columns, etc.
#d,e,object = zeros(shape(data))
def cusumallColumns( data, d, e, object, angled=45, threshd=20, anglee=30, threshe=15):
    r,c=shape(data)
    for i in range(0,c):
        cusumall( data[:,i],d[:,i],e[:,i],object[:,i] )
    


#change: take into account missing data
def trimmedVar( col1 ):
    sz=size(col1)
    realcol1=[];
    for i in range(0,sz):
        if (col1[i]!=-10):
            realcol1.append( col1[i] );
    rsz=size(realcol1);
    #print "real size is ", rsz
    
    if (rsz<6):
        return -1
        
    else:
        
        return var(realcol1[2:rsz-2])
        
        
    
def trimmedMean( col1 ):
    sz=size(col1)
    realcol1=[];
    for i in range(0,sz):
        if (col1[i]!=-10):
            realcol1.append( col1[i] );
    rsz=size(realcol1);
    #print "real size is ", rsz
    
    if (rsz<6):
        return -1
        
    else:
        
        return mean(realcol1[2:rsz-2])
    
 




def trimmedMeanVarPerObject2( col1, object ):   
    sz=size(object)
    mobj=zeros(sz);vobj=zeros(sz);mlist=[];vlist=[];slist=[];mlist0=[];vlist0=[];mlist1=[];vlist1=[]
    running=True
    index=1
    endObject=0
    avgM0=0
    avgV0=0    
    avgM1=0
    avgV1=0
    while running: 
        if (endObject==sz):
            break
        startObject=endObject
        endObject=0
        while endObject==0:
            if (index==sz-1): 
                   endObject=sz
            elif (object[index]==object[startObject]): 
                   index=index+1
            else: 
                   endObject=index      
        vobj  [startObject:endObject]=trimmedVar ( col1[startObject:endObject] )
        mobj  [startObject:endObject]=trimmedMean( col1[startObject:endObject] )
        if ((object[startObject]==0) & (mobj[startObject]!=-1)):
            mlist0.append( mobj  [startObject] )
            vlist0.append( vobj  [startObject] )
        elif ((object[startObject]!=0) & (mobj[startObject]!=-1)):
            mlist1.append( mobj  [startObject] )
            vlist1.append( vobj  [startObject] )
        
        mlist.append( mobj[startObject] )
        vlist.append( vobj[startObject] )
        slist.append( (startObject,endObject) )
        
        #print "Start: ",startObject," End: ",endObject,"Mean : ", mobj[startObject], " Var: ", vobj[startObject], "Sigma :", sqrt(vobj[startObject]), "Obj :", object[startObject]/10
        
    #compute weighted (based on size of region) average mean/variation per region
    divBy0=0;
    divBy1=0;
    print mlist
    for i in range (0,size(mlist)):
        stObj=slist[i][0]
        enObj=slist[i][1]
        if mlist[i]!=-1:
            if object[stObj]==0:
                avgM0=avgM0+(enObj-stObj)*mlist[i]                
                avgV0=avgV0+(enObj-stObj)*vlist[i]
                divBy0=divBy0+(enObj-stObj)
            else:
                avgM1=avgM1+(enObj-stObj)*mlist[i]
                avgV1=avgV1+(enObj-stObj)*vlist[i]
                divBy1=divBy1+(enObj-stObj)
    print avgM0,divBy0
    return double(avgM0)/double(divBy0), double(avgV0)/double(divBy0), double(avgM1)/double(divBy1), double(avgV1)/double(divBy1)
         

    


            
            
#sequence of distances, detection of changes in the variance
#detects low to high
def cusumeDistVar( col1, e, start, end, insigma0s=2.0, insigma1s=100.0 ):
    #sz=size(col1)
    rd=[0.0]
    sigma0s=float(insigma0s)
    sigma1s=float(insigma1s)
    
    for j in range(start, end):
        if (col1[j]!=-1):
            #print (1.0/(2.0*sigma0s) - 1.0/(2.0*sigma1s))
            #print 0.5*log(sigma1s/sigma0s)
            #print col1[j]*col1[j]*(1.0/(2.0*sigma0s) - 1.0/(2.0*sigma1s))
            #debugv[j]=col1[j]*col1[j]*(1.0/(2.0*sigma0s) - 1.0/(2.0*sigma1s))
            val=max(0.0,rd[size(rd)-1] + col1[j]*col1[j]*(1.0/(2.0*sigma0s) - 1.0/(2.0*sigma1s))-0.5*log(sigma1s/sigma0s))
            rd.append(val)
            e[j]=val
        else:
            e[j]=-20
            
            
  

#sequence of distances, detection of changes in the variance
#return position of change from low to high variance
def cusumeDistVarThresh( col1, e, start, end, insigma0s=2.0, insigma1s=25, thresh=3 ):
    #sz=size(col1)
    rd=[0.0]
    sigma0s=float(insigma0s)
    sigma1s=float(insigma1s)
    
    for j in range(start, end):
        if (col1[j]!=-10):
            #print (1.0/(2.0*sigma0s) - 1.0/(2.0*sigma1s))
            #print 0.5*log(sigma1s/sigma0s)
            #print col1[j]*col1[j]*(1.0/(2.0*sigma0s) - 1.0/(2.0*sigma1s))
            #debugv[j]=col1[j]*col1[j]*(1.0/(2.0*sigma0s) - 1.0/(2.0*sigma1s))
            val=max(0.0,rd[size(rd)-1] + col1[j]*col1[j]*(1.0/(2.0*sigma0s) - 1.0/(2.0*sigma1s))-0.5*log(sigma1s/sigma0s))
            rd.append(val)
            e[j]=val
        else:
            e[j]=-3
                        
        if (e[j] > thresh):
             return j
    
    return -1    


#sequence of distances, detection of changes in the variance
#return position of change from high to low
def cusumdDistVarThresh( col1, d, start, end, insigma0s=2.0, insigma1s=25, thresh=12 ):
    #sz=size(col1)
    rd=[0.0]
    sigma0s=float(insigma0s)
    sigma1s=float(insigma1s)
    
    for j in range(start, end):
        if (col1[j]!=-10):
            #print (1.0/(2.0*sigma0s) - 1.0/(2.0*sigma1s))
            #print 0.5*log(sigma1s/sigma0s)
            #print col1[j]*col1[j]*(1.0/(2.0*sigma0s) - 1.0/(2.0*sigma1s))
            #debugv[j]=col1[j]*col1[j]*(1.0/(2.0*sigma0s) - 1.0/(2.0*sigma1s))
            val=max(0.0,rd[size(rd)-1] - col1[j]*col1[j]*(1.0/(2.0*sigma0s) - 1.0/(2.0*sigma1s))+0.5*log(sigma1s/sigma0s))
            rd.append(val)
            d[j]=val
        else:
            d[j]=-3
                        
        if (d[j] > thresh):
             return j
    
    return -1   
    
    
    
#detects low to high to low ... variance
#need to start from low
def cusumallDist( col1, d, e, object, insigma0se=0.1, insigma1se=25, insigma0sd=0.1, insigma1sd=25, threshe=2.5, threshd=2.5):
    running=True
    start=1
    end=size(col1)
    while running:
        startObject=cusumeDistVarThresh( col1, d, start, end, insigma0se, insigma1se, threshe)
        if (startObject==-1):
            running=False
        else:
            endObject=cusumdDistVarThresh( col1, e, startObject+1, end, insigma0sd, insigma1sd, threshd)
            if (endObject==-1):
                object[startObject:end]=10
                invalidate_missing(col1,object,startObject,end);
                running=False
            else:
                object[startObject:endObject]=10
                invalidate_missing(col1,object,startObject,endObject);
                start=endObject
        
    invalidate_missing(col1,object,0,end)
    
    
    
    
#For ALL    
#data is a set of columns, etc.
#d,e,object = zeros(shape(data))
def cusumallDistColumns( data, d, e, object, insigma0se=2.0, insigma1se=100.0, insigma0sd=2.0, insigma1sd=100.0, threshe=2.5, threshd=12):
    r,c=shape(data)
    for i in range(0,c):
        cusumallDist( data[:,i],d[:,i],e[:,i],object[:,i], insigma0se, insigma1se, insigma0sd,insigma1sd, threshe,threshd)



#def cusumdVar( col1, d, start, end, sigma0s=2.0, sigma1s=100.0 ):
#    for n in range(start,end):
        #print sigma0s, sigma1s, sigma1s/sigma0s
        #print 0.5*log(sigma1s/sigma0s)
        #print d[n-1]+col1[n]*col1[n]*(1.0/(2.0*sigma0s)-1.0/(2.0*sigma1s))-0.5*log(sigma1s/sigma0s)
#        d[n]=max(0,d[n-1]+col1[n]*col1[n]*(1.0/(2.0*sigma0s)-1.0/(2.0*sigma1s))-0.5*log(sigma1s/sigma0s))
        
    


     
    

#converts a signed angle vector to +1,-1,-10
def convert_signedangles_to_signs(col):
    scol=zeros(size(col))
    for i in range(0,size(col)):
        #missing becomes -10 (to be used by next routine)
        if (col[i]==360):
            scol[i]=-10
        elif ( (sign(col[i])==1) | (sign(col[i])==0)  ):
            scol[i]=1
        else:
            scol[i]=-1
            
    return scol

#removes -10 from col
def delete_missing(col):
    scol=[]
    for i in range(0,size(col)):
        if col[i]!=-10:
            scol.append(col[i])
    return scol


    
#removes 360 from col (for signed angles)
def delete_missing2(col):
    scol=[]
    for i in range(0,size(col)):
        if col[i]!=360:
            scol.append(col[i])
    return scol


            
 
def comp_missing2(col):
    nm=0
    for i in range(0,size(col)):
        if col[i]==360:
            nm+=1
    return nm

#Three StateHMM----->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#simple gaussian
def gaussian_simple( Xt, Mvalue, sigma=10.0 ):
    #sigma=10.0
    sigma2=sigma*sigma;
    return longdouble(1.0/(sigma*sqrt(2.0*pi)))*exp( -(Xt-Mvalue)*(Xt-Mvalue) / (2.0*sigma2) )


    
def transformLLR(LLRin,LLRout, maprt):
    N=size(maprt)
    for rt in range(0,N):
        if maprt[rt]!=-10:
            LLRout[ rt ] = LLRin[ maprt[rt] ]
        else:
            LLRout[ rt ] = -3


#initialize apH, apK, at time t, objervation vector col
#ps[i]=1/N for i=0..N-1 (1/2 if N==2, 1/3 if N==3)...

def initialize_apH_apK_Ns( t, col, apH, apK,  Means, Sigmas, ps, N):
    for i in range(0,N):        
        a=ps[i]*gaussian_simple(col[t], Means[i], Sigmas[i])
        apH[t][i]=a
        apK[t][i]=a
    
    #print "Initialization at time ", t, "apH is ", apH[t], "apK is ", apK[t]

#update ap (H or K) at time t from previous value (time t-1)
#col is vector of objervations
# p: probability on->off
# q: probability off->on
# works only for 3 states now
def update_ap_Ns( t, col, ap, pqr, Means, Sigmas, Nstates ):
    tempB=0
    for i in range(0,Nstates):
        tempB += longdouble(ap[t-1][i]);
    #print("A")
    #print pqr[0],pqr[1],pqr[2]
    tempA=longdouble(ap[t-1][0]*pqr[0]) + longdouble(ap[t-1][1]*pqr[1]) + longdouble(ap[t-1][2]*pqr[2])
    tempA*=longdouble(gaussian_simple(col[t],Means[0],Sigmas[0]))
    ap[t][0]=longdouble(tempA)/longdouble(tempB)
    #print("B")
    #print pqr[3],pqr[4],pqr[5]
    tempA=longdouble(ap[t-1][0]*pqr[3]) + longdouble(ap[t-1][1]*pqr[4]) + longdouble(ap[t-1][2]*pqr[5])
    tempA*=longdouble(gaussian_simple(col[t],Means[1],Sigmas[1]))
    #print("C")
    ap[t][1]=longdouble(tempA)/longdouble(tempB)
    
    np1=longdouble(1.0-pqr[0]-pqr[3])
    np2=longdouble(1.0-pqr[1]-pqr[4])
    np3=longdouble(1.0-pqr[2]-pqr[5])
    
    #print "npi:", np1, np2, np3
    tempA=longdouble(ap[t-1][0]*np1) + longdouble(ap[t-1][1]*np2) + longdouble(ap[t-1][2]*np2)
    tempA*=longdouble(gaussian_simple(col[t],Means[2],Sigmas[2]))
    
    ap[t][2]=longdouble(tempA)/longdouble(tempB)
    
    #print "ap:", ap[t][0],ap[t][1],ap[t][2]
    
    



#calculate the ln ration at t, given apH and apK
def calc_ln_ratio_Ns(t, apH, apK,Nstates):
    scale=longdouble(1.0);
    tempA=0
    tempB=0
    for i in range(0,Nstates):
        tempA += longdouble(scale*apK[t][i])
        tempB += longdouble(scale*apH[t][i])
    tempA = longdouble( abs(tempA) )
    tempB = longdouble( abs(tempB) )
    #print tempA, tempB
    #print tempA/tempB
    #print log( longdouble(tempA/tempB) )
    return log( tempA / tempB  )
    


    
    
    
    
                
#col1 is an array of measurements
#(p1,q1) are the transition probabilities under hypothesis H (p1: on (object) -> off (ground)
#                                                            (q1: off(ground) -> on (object)
#(p2,q2) are the transition probabilities under hypothesis K (p2: on -> off, q2: off -> on)
#stop if exceed threshold h
#Return LLR as well (should be zeros(size(col)))
#it skips missing data
#LLR gets the value -3 in missing data (just for plotting-it does not have any effect on the algorithm)
#
#The last parameter is called reverse: If it is False, the usual technique is used (default)
# otherwise a reverse situation is implemented

#Three States
  #prqH=[p1,q1,r1,p2,q2,r2] for H 
#prqK=      -//-          for K
def N_state_CUSUM_HMM_w_missing( col, LLRout, h=100, Means=[90.0,15.0,-90.0], Sigmas=[12.0,2.0,2.0], pqrH=[0.9,0.1,0.0,0.1,0.9,0.0], pqrK=[1.0/3.0,1.0/3.0,1.0/3.0,1.0/3.0,1.0/3.0,1.0/3.0],Nstates=3):
    #two states
    #two hypotheses: H, K
    #N observations
    
    N=size(col)
    
    LLR=zeros( N );
    #Null hypothesis
    
    apH=zeros((N,Nstates))
    apK=zeros((N,Nstates))
    init_probs=[]
    for i in range(0,Nstates):
        init_probs.append(1.0/double(Nstates))    
    needInit=True
    #time
    t=-1
    #map from rt to t
    maprt=zeros(N);
    #map from t to rt
    mapt=[]; scol=delete_missing2(col)
    for rt in range(0,N):
        #if col[rt]==360.0:
            #skipping missing
            #print "col rt ", rt, "is ", col[rt]; maprt[rt]=-10
        if col[rt]!=360.0:
            t=t+1
            maprt[rt]=t
            mapt.append(rt)
            #print "scol ", t , " is ", scol[t]
            if needInit:
                initialize_apH_apK_Ns( t, scol, apH, apK, Means, Sigmas,init_probs, Nstates )
                needInit=False
            
            if t==0:
                LLR[0]=calc_ln_ratio_Ns(0, apH, apK, Nstates);
                #debugA[0]=LLR[0]
                print(LLR[0])
            else:
                update_ap_Ns(t, scol, apH, pqrH, Means, Sigmas, Nstates)
                update_ap_Ns(t, scol, apK, pqrK, Means, Sigmas, Nstates)
                tempVal=longdouble( calc_ln_ratio_Ns(t, apH, apK, Nstates) )
                #print tempVal
                #debugA[t]=tempVal;
            
                if isnan(tempVal):

                    print "Hit nan:"
                    transformLLR(LLR,LLRout,maprt);
                    #return t
                    return rt
                
                LLR[t]=LLR[t-1]+tempVal
                #print "LLR: ", LLR[t]
                #debugA[t]=LLR[t];
            
            if LLR[t]>h:
                transformLLR(LLR,LLRout,maprt)
                #print "real t ", rt, "wout missing t ", t
                return rt
            
            elif LLR[t]<0.0:
                LLR[t]=0
                #initialize_apH_apK(t, col, apH, apK )
                needInit=True
    
    #It would return the end of the region
    transformLLR(LLR,LLRout,maprt);
    #return t
    return rt
  


def N_state_SPRT_w_missing( col, LLRout, hp=20, hn=-20, Means=[90.0,15.0,-90.0], Sigmas=[12.0,2.0,2.0], pqrH=[0.9,0.1,0.0,0.1,0.9,0], pqrK=[1.0/3.0,1.0/3.0,1.0/3.0,1.0/3.0,1.0/3.0,1.0/3.0],Nstates=3):
    #Nstates states
    #two hypotheses: H, K
    #N observations
    
    N=size(col)
    
    LLR=zeros( N );
    #Null hypothesis
    
    apH=zeros((N,Nstates))
    apK=zeros((N,Nstates))
    init_probs=[]
    for i in range(0,Nstates):
        init_probs.append(1.0/double(Nstates))    
    needInit=True
    #time
    t=-1
    #map from rt to t
    maprt=zeros(N);
    #map from t to rt
    mapt=[]; scol=delete_missing2(col)
    for rt in range(0,N):
        if col[rt]==360:
            #skipping missing
            maprt[rt]=-10
            continue
        else:
            t=t+1
            maprt[rt]=t
            mapt.append(rt)
         
        if needInit:
            initialize_apH_apK_Ns( t, scol, apH, apK, Means, Sigmas,init_probs, Nstates )
            needInit=False
            
        if t==0:
            LLR[0]=calc_ln_ratio_Ns(0, apH, apK, Nstates);
            #debugA[0]=LLR[0]
            #print(LLR[0])
        else:
            update_ap_Ns(t, scol, apH, pqrH, Means, Sigmas, Nstates)
            update_ap_Ns(t, scol, apK, pqrK, Means, Sigmas, Nstates)
            tempVal=longdouble( calc_ln_ratio_Ns(t, apH, apK, Nstates) )
            #print tempVal
            #debugA[t]=tempVal;
            
            if isnan(tempVal):

                print "Hit nan:"
                transformLLR(LLR,LLRout,maprt);
                #return t
                return rt
                
            LLR[t]=LLR[t-1]+tempVal
            #print "LLR: ", LLR[t]
            #debugA[t]=LLR[t];
            
        if LLR[t]>hp:
            transformLLR(LLR,LLRout,maprt)
            #return t
            #stop due to positive threshold
            #print "real t ", rt, "wout missing t ", t
            return 1,rt
            
        elif LLR[t]<hn:
            #stop due to negative threshold
            transformLLR(LLR,LLRout,maprt);
            #print "real t ", rt, "wout missing t ", t
            return -1,rt
    
    #It would return the end of the region
    transformLLR(LLR,LLRout,maprt);
    #stop due to end of data
    return 0,rt
    
    
def N_state_CUSUM_HMM_followed_by_SPRT_w_missing( col,  h=10, hp=10, hn=-5, Means=[90.0,15.0,-90.0], Sigmas=[12.0,2.0,2.0], pqrH=[0.9,0.1,0.0,0.1,0.9,0], pqrK=[1.0/3.0,1.0/3.0,1.0/3.0,1.0/3.0,1.0/3.0,1.0/3.0],Nstates=3):
    LLRout=zeros(size(col))
    N=size(col)
    t0=0
    continueLoop=True
    a=-1
    t=-1
    rcol=col
    while continueLoop:
        if size(rcol)==0:
            print "Uknown condition (4)"
            return -1,-1; #a,size(col)-1
            
        a=N_state_CUSUM_HMM_w_missing( rcol, LLRout, h, Means, Sigmas, pqrH, pqrK,Nstates)
        #print "CUSUM HMM goes off at time ", a+t0
        if a==N-1:
            print "No tree"
            return -1,-1
     
        #check whether SPRT becomes negative
        rcol=rcol[a+1:]
        t0 += a+1
        #run SPRT
        
        if size(rcol)==0:
            print "Uknown condition (5)"; return -1, -1
        val,t=N_state_SPRT_w_missing(rcol,LLRout,hp,hn,Means, Sigmas, pqrH, pqrK,Nstates)
        #print "SPRT after CUSUM HMM: ", val, "at time ", t+t0
        if val==0:
            #print "Unknown condition"
            return -1,-1; #a,t+t0
            
        elif val==-1:
             
            #print "Possible obstacle between", t0-1, " and", t+t0
            rcol=rcol[t+1:]
            t0 += t+1
        else:
            a=t0-1
            #print "tree starts at", a
            #print "exit first loop..."
            continueLoop=False
    
    while t<N:
        rcol=rcol[t+1:]
        #plot(rcol,hold=False)
        #print "rcol has size", size(rcol)
        if size(rcol)==0:
            #print "Uknown condition (3)"
            return a, t+t0
            
        t0 += t+1
        val,t=N_state_SPRT_w_missing(rcol,LLRout,hp,hn,Means, Sigmas, pqrH, pqrK,Nstates)
        #print "SPRT after CUSUM HMM (2): ", val, "at time ", t + t0
        if val==0:
            #print "Uknown condition (2)"
            return a,t+t0
            
        elif val==-1:
            #print "Tree ends at", t+t0
            return a,t+t0
        
        #else: print "restart at", t+t0
                    


def detect_trees_using_N_CUSUM_HMM(col):
    object=zeros(size(col))
    N=size(col)
    rcol=col
    continueLoop=True
    t0=0
    prevEndOfTree=-1
    while continueLoop:
        if size(rcol)==0:
            #print "rcol is empty"
            continueLoop=False
        else:
            #print "rcol has size ",size(rcol)
            #plot(rcol,hold=False)
            a,b=N_state_CUSUM_HMM_followed_by_SPRT_w_missing(rcol)
            #print "N_state_CUSUM_HMM returned: ", a, b
            if a==-1:
                continueLoop=False
            else:
            #tree between [a, b]
                #print "check for first negative from ", a, " to ", b, " size(rcol) is ", size(rcol)
                c=check_for_first_negative(rcol, a, b)
                #print "first negative at ", c
                if c!=-1: b=c+1
                object[a+t0:b+t0+1]=1
                #print "pETree is", prevEndOfTree, "a+t0 is ", a+t0
                if prevEndOfTree!=-1:
                    if a+t0-prevEndOfTree<5:
                        object[prevEndOfTree:a+t0]=1
                prevEndOfTree=b+t0+1
                
                if b+t0<N:
                    rcol=rcol[b+1:]
                    t0=t0+b+1
                else:
                    continueLoop=False
    
    return object



#return first negative in b, b-1, ..., a+1
#if not found returns -1
def check_for_first_negative(col, a, b):
    for i in range(-b,-a):
        j=int(-i)
        if col[j] < 0.0: return j
    return -1
    
    
    
def detect_trees_all_columns(data,limit=0):
    r,c=shape(data)
    object=zeros((r,c))
    if limit==0: limit=c
    for i in range(0,min(limit,c)):
        print "Column ", i
        object[:,i]=detect_trees_using_N_CUSUM_HMM( data[:,i] )
    return object


#trees is the array returned by detect_trees_all_columns
#it has ones in tree regions
#and 0 in non-tree regions
def calculate_vert_horizontal_all(data, trees, limit=0):
    r,c=shape(data)
    object=zeros((r,c))
    if limit==0: limit=c
    for i in range(0,min(limit,c)):
        print "Column ", i
        object[:,i]=calculate_vert_horizontal_scan_line(data[:,i],trees[:,i])
    return object
    

def modify_obj_for_plotting(obj):
    r,c=shape(obj)
    for i in range(0,c):
        for j in range(0,r):
            if obj[j,i]==1: obj[j,i]=-5
            if obj[j,i]==10: obj[j,i]=-10

    
def calculate_vert_horizontal_scan_line(col, trees):
    object=zeros(size(col))
    N=size(col)
    startObject=0
    endObject=0
    ncol=zeros(size(col))
    for i in range(0,N):
        if ncol[i]==360: ncol[i]=-10
        else: ncol[i]=abs(col[i])
    print "start: ", startObject
    continueLoop=True
    while continueLoop:
        for i in range(startObject,N):
            if trees[i]==0: continue
            else: break
        #here i is 1
        if trees[i]==0: endObject=N; continueLoop=False
        else: endObject=i
        print "end :", endObject
        if startObject<endObject:
            col1=ncol[startObject:endObject]
            d=zeros(size(col1))
            e=zeros(size(col1))
            print "cusumall from ", startObject, " to ", endObject
            cusumall( col1, d, e, object[startObject:endObject] )
        if continueLoop:
            for j in range(endObject,N):
                if trees[j]==1: object[j]=1; continue
                else: break
            if trees[j]==1: continueLoop=False
            else: 
                startObject=j
                print "start ", startObject
    
    return object
    
    

#Simple union-find algorithm
def map_sets(union_array):
    set_array=zeros(size(union_array))
    set_num=1
    for i in range(0,size(union_array)):
        if union_array[i] < 0:
            set_array[i]=set_num
            set_num+=1
    return set_array
    

def find_final( union_array, set_array, A):
    if union_array[A] < 0: return set_array[A]
    else: return find_final( union_array, set_array, union_array[A] )

def find( union_array, A):
    if union_array[A] < 0: return A
    else: return find( union_array, union_array[A] )

    
def unite(union_array, A, B):
    #print "Unite: ", A, " with ", B
    setA=find( union_array, A)
    setB=find( union_array, B)
    if setA!=setB:
        C=min(setA,setB)
        D=max(setA,setB)
        union_array[C] += union_array[D]
        union_array[D]=C
        
        
#region_growing
#sAngles: signed angles
#distT  : distance tilda
#dist   : distance (regular)
#Tthresh : thrshold of distT
#Ttrhesh2 : threshold for dist
#vht: 0 if horizontal, 10 if vertical, 1 if tree

def region_growing_vert_horiz(sAngles, distT, dist, vht,limit=0,Tthresh=2.0,Tthresh2=2.0,sClust=20):
    r,c=shape(sAngles)
    obj=zeros((r,c))
    union_array=[-1]
    if limit==0: limit=c
    regionNum=0
    #first column
    #print "A"
    for j in range(0,r):
        Aj=vht[j,0]; Dj=sAngles[j,0]
        if ( (Aj!=1) & (Dj!=360) ):
            #not a tree and not missing data
            if j==0:
                regionNum+=1
                obj[j,0]=regionNum
                union_array.append(-1)
            else: 
                if ( (Aj==vht[j-1,0]) & (sAngles[j-1,0]!=360) & (obj[j-1,0]!=0) & (dist[j-1,0]<Tthresh2) ): obj[j,0]=obj[j-1,0]
                else:
                    regionNum+=1
                    obj[j,0]=regionNum
                    union_array.append(-1)
    #remaining columns
    for i in range(1,min(limit,c)):
        print "Column : ", i
        for j in range(0,r):
            Aji=vht[j,i]; Dji=sAngles[j,i]; 
            if ( (Aji!=1) & (Dji!=360) ):
                #not a tree and not missing data
                if j==0:
                    #look left
                    if Aji==vht[j,i-1]:
                        #check Dtilda
                        dT=distT[j,i-1]
                        if ((dT<Tthresh) & (dT!=-10) & (obj[j,i-1]!=0) ):
                            obj[j,i]=obj[j,i-1]
                        else:
                            regionNum+=1
                            obj[j,i]=regionNum
                            union_array.append(-1)
                    else:
                        regionNum+=1
                        obj[j,i]=regionNum
                        union_array.append(-1)
                else:
                    dT=distT[j,i-1]
                    if ((dT>=Tthresh) | (dT==-10)):
                        #use only previous
                        if ( (Aji==vht[j-1,i]) & (sAngles[j-1,i]!=360) & (obj[j-1,i]!=0) & (dist[j-1,i]<Tthresh2) ): obj[j,i]=obj[j-1,i]
                        else: 
                            regionNum+=1
                            obj[j,i]=regionNum
                            union_array.append(-1)
                    else:
                        #here we have two pieces of information
                        #previous in same scanline
                        Oprev=obj[j-1,i]
                        Aprev=vht[j-1,i]
                        #left
                        Oleft=obj[j,i-1]
                        Aleft=vht[j,i-1]
                        candPrev= ((Oprev!=0) & ( sAngles[j-1,i]!=360 ) & (Aprev==Aji) & (dist[j-1,i]<Tthresh2) )
                        candLeft= ((Oleft!=0) & ( sAngles[j,i-1]!=360 ) & (Aleft==Aji) & (distT[j,i-1]<Tthresh) )
                        if ((candPrev==False) & (candLeft==False)):
                            regionNum+=1
                            obj[j,i]=regionNum
                            union_array.append(-1)
                        elif ( (candPrev==False) & (candLeft==True) ):
                            #if ( (Aji==Aleft) & (sAngles[j,i-1]!=360) ):
                            if (Oleft==0): print "Error Oleft"
                            obj[j,i]=Oleft
                            #else:
                            #    regionNum+=1
                            #    obj[j,i]=regionNum
                            #    union_array.append(-1)
                        elif ( (candLeft==False) & (candPrev==True) ):
                            #if ( (Aji==Aprev) & (sAngles[j-1,i]!=360) ):
                            if (Oprev==0): print "Error Oprev"
                            obj[j,i]=Oprev
                            #else:
                            #    regionNum+=1
                            #    obj[j,i]=regionNum
                            #    union_array.append(-1)
                        else:
                            #both True
                            #if Oprev==Oleft: obj[j,i]=Oprev
                            #else:
                            #obj[j,i]=min(Oleft,Oprev)
                            #unite(union_array, int(Oleft), int(Oprev))
                            if Aprev != Aleft:
                                #regions of different types
                                #merge with previous
                                if (obj[j-1,i]==0): print "Error obj_j1"
                                obj[j,i]=obj[j-1,i]
                            else:
                                #same types
                                if min(Oleft,Oprev)==0: print "Error min"
                                obj[j,i]=min(Oleft,Oprev)
                                unite(union_array, int(Oleft), int(Oprev))
#connect united
    set_array=map_sets(union_array)
    #sizes=zeros((r,c))
    max_region=0
    for i in range(0,min(limit,c)):
        for j in range(0,r):
            if obj[j,i]!=0:
                #if (int(abs(find(union_array,int(obj[j,i])))) < 30): obj[j,i]=1
                #else:
                #sizes[j,i]=union_array[abs(find(union_array,int(obj[j,i])))]
                obj[j,i]=find_final(union_array, set_array, int(obj[j,i]))
                if obj[j,i]>max_region: max_region=obj[j,i]
            if ( (vht[j,i]==1) & (obj[j,i]!=0) ): print "problem "
    print "Max Region: ", max_region
    sizes_array=zeros(max_region+1)
    for i in range(0,min(limit,c)):
        for j in range(0,r):
                sizes_array[int(obj[j,i])] +=1
    #suppress small regions
    for i in range(0,min(limit,c)):
        for j in range(0,r):
            if obj[j,i]!=0:
                if sizes_array[int(obj[j,i])] < sClust: obj[j,i]=1
    return obj


#MERGEREGIONS-------------------------------------------------------------------------------
#region is a 2D array
#num is the column of region we are looking at
#maxRegion is the absolute value of last region
def mergeregions(prevcol, curcol, maxRegion,theRegions,colval):
    M=size(curcol)
    for N in range(0,size(prevcol)):
        if prevcol[N]==0: 
            break;
    #N is the start of the tree on the first scan
    count=0
    countM=0
    #sizethresh as a percentage (1%) of smallest region
    #print "previous size is ", M
    #print "current size is ", N
    
    #threshsize=int(0.05*min(N,M))
    
    #print "threshsize is ", threshsize
    
    similarity=False
    #if max(N,M) - min(N,M) < threshsize:
    #    print "passed threshsize "
    for i in range(0,min(N,M)):
        if ( (prevcol[i]<=-4) & (curcol[i]==0) ):
            count += 1
        elif ( (prevcol[i]>=4) & (curcol[i]==10) ):
            count += 1
        elif ( (prevcol[i]==-1) | (curcol[i]==-10) ):
            countM += 1
    #print "similarity count is ", count
    actualsize=min(N,M) - countM
    #print "actualsize is ", actualsize
    if count > int(0.9*actualsize):
        similarity=True
    
    #there is similarity
    #do simple assignment
    minimumRegionNum=-1
    if similarity==True:
        for i in range(0, min(N,M)):
            if ( (prevcol[i]<=-4) & (curcol[i]==0) ):
                curcol[i]=prevcol[i]
                if minimumRegionNum==-1:
                    minimumRegionNum=abs(curcol[i])
                theRegions[int(abs(curcol[i])-4)][3].append( colval[i] )
            elif ( (prevcol[i]>=4) & (curcol[i]==10) ):
                curcol[i]=prevcol[i]
                if minimumRegionNum==-1:
                    minimumRegionNum=abs(curcol[i])
                theRegions[int(abs(curcol[i])-4)][3].append( colval[i] )
            elif curcol[i]==-10:
                curcol[i]=-1
            else:
                curcol[i]=-2
        for i in range(min(N,M), M):
            curcol[i]=-2
        return maxRegion,minimumRegionNum
        
    else:
        
        return initial_assignment(curcol,maxRegion,theRegions,colval),maxRegion+1
                
                
    
def calculate_MeanVar_of_Regions( theRegions, start, end ):
    lowmeanT=0;highmeanT=0;lowvarT=0;highvarT=0;
    numTlow=0; numThigh=0;
    for i in range(start, end):
        reg=theRegions[i]
        
        reg[1]=trimmedMean( reg[3] )
        reg[2]=trimmedVar( reg[3] )
        if ( isnan(reg[1]) | isnan(reg[2]) ):
            reg[1]=-1
            reg[2]=-1
        if ( (reg[1]!=-1) & (reg[2]!=-1) ):
            sz=double(size(reg[3]))
            #sz=1
            if reg[0]==0:
                numThigh += sz
                highmeanT += sz*reg[1]
                highvarT  += sz*reg[2]
                #print "highMeanT is ", highMeanT, "highVarT is ", highVarT
            else:
                numTlow += sz
                lowmeanT += sz*reg[1]
                lowvarT  += sz*reg[2]
                #print "lowMeanT is ", lowMeanT, "lowVarT is ", lowVarT
    highmean=-1;highvar=-1;lowmean=-1;lowvar=-1;
    
    if numThigh>0:
        
        highmean=double(highmeanT)/double(numThigh)
        highvar  =double(highvarT)/double(numThigh)
        
    if numTlow>0:
        lowmean=double(lowmeanT)/double(numTlow)
        lowvar  =double(lowvarT)/double(numTlow)       
        
    return lowmean,lowvar,highmean,highvar

#provides assignment in scanline without merging
#gets previous max. region (absolute value)
#returns current max. region
#i.e. you need to increase maxRegion by one in the next iteration
def initial_assignment(curcol,maxRegion,theRegions,colval):#,theRegions,angleval):
    prevRegion=0
    maxRegion += 1
    for i in range(0,size(curcol)):
        if (curcol[i]==0):
            if prevRegion==0:
                prevRegion=-maxRegion
                curcol[i]  = -maxRegion
                #Create a new region
                theRegions.append( [0,-2,-2,[colval[i]]] )
            elif prevRegion<0:
                curcol[i]= prevRegion
                #print int(abs(curcol[i])-4)
                theRegions[int(abs(curcol[i])-4)][3].append( colval[i] )
            else: #prevRegion>0
                prevRegion = -prevRegion - 1
                curcol[i]=prevRegion
                theRegions.append([0,-2,-2,[ colval[i] ]])
        elif (curcol[i]==10):
            if prevRegion==0:
                curcol[i]  = maxRegion
                prevRegion = maxRegion
                theRegions.append( [1,-2,-2,[colval[i]]] )
            elif prevRegion>0:
                curcol[i]= prevRegion
                #print int(abs(curcol[i])-4)
                theRegions[int(abs(curcol[i])-4)][3].append( colval[i] )
            else: #prevRegion<0
                prevRegion = -prevRegion + 1
                curcol[i]=prevRegion
                theRegions.append([1,-2,-2,[ colval[i] ]])
        elif curcol[i]==-10:
            curcol[i]=-1
    if prevRegion==0:
        
        return maxRegion-1 
    
    else:
        
        return abs(prevRegion)



def ok_res(res):
    if ( (isnan(res[0])==False) & (isnan(res[1])==False) ):
        a=res[0];b=res[1]
    else:
        a=-1;b=-1;
    if ( (isnan(res[2])==False) & (isnan(res[3])==False) ):
        c=res[2];d=res[3]
    else:
        c=-1;d=-1;
        
    return a,b,c,d

#def calculate_means_var_from_regions( region, i ):
    

    
#main function
#merges regions and calculates means, and variances of them
def cusumallColumnsHMM_compMeanVar( data, end, update=False, angled=45, threshd=20, anglee=30, threshe=15, lowMean=15, highMean=90, lowSigma=12.0, highSigma=2.0):
    plMean=lowMean;phMean=highMean;plSigma=lowSigma;phSigma=highSigma;
    r,c=shape(data)
    m0=[];v0=[];m1=[];v1=[];
    limit=min(end,c)
    #holds the region information
    #-r -> low of region r-3 (r>=4), +r-> high of region r, 0: start of tree, -1:missing data, -2: unknown
    #region=ones(shape(data))
    #region *= -2
    region=zeros(shape(data))
    firstTime=True
    theRegions=[]
    #list of regions with info
    # for each region r the list should contain:
    #theRegions[ abs(r)-4 ] should contain: (0/1, mean, var, [list of data values for region])
        
    for i in range(0,limit):
        LLR=zeros(r)
        stop=two_state_HMM_w_missing(data[:,i],LLR,6,plMean,phMean,plSigma,phSigma)
        plot(LLR);
        print "stop in scanline ", i , " is ", stop
        if stop>10:
            col=data[0:stop,i]
            d=zeros(size(col))
            e=zeros(size(col))
            object=zeros(size(col))
            cusumall( col,d,e,object )
            if firstTime==True:
                #plot(object)
                maxRegion=initial_assignment( object, 3, theRegions, col )
                print "calculate_MeanVar_of_Regions from 0 to ", maxRegion+1-4
                res=calculate_MeanVar_of_Regions( theRegions, 0, maxRegion+1-4 )
                #print res
                    
                #print "maxRegion in scanline ", i , " is ", maxRegion
                #plot(object)
                #region[0:size(object),i]=object
                firstTime=False
            else:
                #print "Merge scanline ", i
                maxRegion,minRegion=mergeregions(region[:,i-1],object,maxRegion, theRegions, col)
                print "calculate_MeanVar_of_Regions from ", int(minRegion-4), "to ", int(maxRegion+1-4)
                res=calculate_MeanVar_of_Regions( theRegions, int(minRegion-4), int(maxRegion+1-4) )
                #print res
                #print "maxRegion in scanline ", i , " is ", maxRegion
            
            #update region
            region[0:size(object),i]=object
            #place marker for tree
            if stop<r:
                region[size(object),i]=0
                
            #update means/vars for HMM
            res=ok_res(res)
            #print res
            if (update==True):
                if ((res[0]!=-1) & (res[1]!=-1)):
                    plMean=res[0];plSigma=sqrt(res[1])
                if ((res[2]!=-1) & (res[3]!=-1)):
                    phMean=res[2];phSigma=sqrt(res[3])
                
            print "High-mean: ", phMean, " High-std : ", phSigma, "Var ( ", phSigma*phSigma, " )"
            print " Low-mean: ", plMean, " Low-std : ", plSigma, "Var ( ", plSigma*plSigma, " )"
            
    return region,theRegions
                    
                    