import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np


def credibleInterval(df,error):
    percentile=df.quantile(q=0.5+error/2)-df.quantile(q=0.5-error/2)
    p=np.abs(np.abs(df.quantile(q=0.5+error/2)) - np.abs(df.quantile(q=0.5-error/2))) == np.abs(percentile)

    return p.values[0]


def makeDfForPlot(modelName,imputeType,averageOrPeriod,adultAge):
    
    error = 0.95
    
    if averageOrPeriod == 'average':

        startDir = 'results/' + modelName + '/MultiLevel/'+imputeType+'/'+averageOrPeriod +'/'+adultAge +'/'

        l='djub'

        M=[]

        lag=sorted(os.listdir(startDir))
        lag=sorted([int(float(ll)) for ll in lag])
        lag=[str(ll) for ll in lag]
        for v in lag:
            VAR = os.listdir(startDir + v + '/')
            for var in VAR:
                BETA=os.listdir(startDir + v + '/' + var + '/')
                DEP=np.asarray(BETA)
                ii=[s.startswith('beta') for s in DEP]

                DEP=list(DEP[ii])
                DEPlab=[s[-3:] for s in DEP]

                for dep in DEP:
                    fn=startDir + v + '/' + var + '/' + dep
                    Z=pd.read_csv(fn,index_col=0)

                    m = Z.median().values[0]

                    sem = Z.std().values[0]

                    p = credibleInterval(Z,error)


                    M.append([l,int(v),var,dep,m,sem,p])


        N=pd.DataFrame(M,columns=['value','lag','targ','pred','M','sem','p'])

    
    if averageOrPeriod == 'period':
        startDir = 'results/' + modelName + '/MultiLevel/'+imputeType+'/'+averageOrPeriod +'/'+adultAge +'/'

        l='djub'

        M=[]

        PER=sorted(os.listdir(startDir))


        for per in PER: 


            lag=sorted(os.listdir(startDir + per + '/'))
            lag=sorted([int(float(ll)) for ll in lag])
            lag=[str(ll) for ll in lag]


            for v in lag:
                VAR = os.listdir(startDir + per + '/' + v + '/')
                for var in VAR:
                    BETA=os.listdir(startDir + per + '/' + v + '/' + var + '/')
                    DEP=np.asarray(BETA)
                    ii=[s.startswith('beta') for s in DEP]

                    DEP=list(DEP[ii])
                    DEPlab=[s[-3:] for s in DEP]

                    for dep in DEP:
                        fn=startDir + per + '/' + v + '/' + var + '/' + dep
                        Z=pd.read_csv(fn,index_col=0)

                        m = Z.median().values[0]

                        sem = Z.std().values[0]

                        p = credibleInterval(Z,error)


                        M.append([l,per,int(v),var,dep,m,sem,p])


        N=pd.DataFrame(M,columns=['value','period','lag','targ','pred','M','sem','p'])



    return N




def plotRow(av3,br2,pred,row,col,averageOrPeriod):
    minV = np.asarray([0])
    maxV = np.asarray([0])

    av4 = av3[av3.loc[:,'pred'] == pred]
    br3 = br2[br2.loc[:,'pred'] == pred]

    num_t = len(br3.loc[:,'lag'])
    pert = (np.random.random(num_t)-0.5)/3

    t=br3.loc[:,'lag'] + pert
    vals = br3.loc[:,'M']
    sem = br3.loc[:,'sem'] 

    minV = np.append(minV,vals-sem)
    maxV = np.append(maxV,vals+sem)

    if averageOrPeriod == 'period': 
        ax[row,col].errorbar(t,vals,yerr=sem, fmt='o', color='r',alpha=0.6)

    
    t=sorted(av4.loc[:,'lag'].values)

    
    vals = av4.loc[:,'M'].values
    sem = av4.loc[:,'sem'] .values

    minV = np.append(minV,vals-sem*1.1)
    maxV = np.append(maxV,vals+sem*1.1)


    p = av4.loc[:,'p'] .values
    

    for iii in range(len(t)):        
        if p[iii] == True: alphaMinus=0.
        if p[iii] == False: alphaMinus=0.45


        ax[row,col].bar(t[iii],vals[iii],width=0.9999, color='k', yerr=sem[iii],
                   error_kw=dict(ecolor='k', lw=2.5, capsize=2, capthick=1.,alpha=0.9-alphaMinus),alpha=0.7-alphaMinus)

        
    if len(av4.loc[:,'lag']) > 8:
        ax[row,col].set_xticks([5,15,25])
        ax[row,col].set_xticklabels(['5 yrs','15 yrs','25 yrs'],fontsize=11)
        
        
    if len(av4.loc[:,'lag']) < 8:
        ax[row,col].set_xticks([1,2,3])
        ax[row,col].set_xticklabels(['10 yrs','20 yrs','30 yrs'],fontsize=11)

    #set y ticks and limits
    ymin = round(np.min(minV),2) 
    ymax = round(np.max(maxV),2)

    ax[row,col].set_ylim([-0.75,0.75])

    return



            
def setRowColumnLabels(var,row,col,iftarget):

    var=var#lookup.loc[var]
    
    if iftarget==True:
        ax[row,col].text(0.5,0.5,r"$" + var + "_t = $", horizontalalignment='center', verticalalignment='center', fontsize=33)

    if iftarget==False: 
        ax[row,col].text(0.5,0.5,r"$" + var + "_{t-g}$", horizontalalignment='center', verticalalignment='center', fontsize=33)

    ax[row,col].set_yticks([])
    ax[row,col].set_xticks([])  

    return

modelName = 'CivicValuesDemocracyCivicCosmopolitan'

imputeType = 'linear' #'quad'
averageOrPeriod = 'average'  #'average'
adultAge='0'

#targets = ['DEM']
#preds = ['betaDEM','betaCRI','betaTRU','betaSUP','betaCOS','betaGDP','betaCxD']#, 'betaCxD']#None#['betaCPL', 'betaONE', 'betaGDP']
#'betaAUT',
targets = ['DEM']

preds = ['betaDEM','betaCOS','betaCON','betaTRU','betaSUP','betaGDP']
#preds = ['betaCPL','betaONE','betaGDP','betaDEM']

#targets = ['CPL','ONE','GDP']
#preds = ['betaCPL','betaONE','betaGDP']#,'betaDEM','betaCON','betaLRP']#, 'betaCxD']#None#['betaCPL', 'betaONE', 'betaGDP']
#preds = ['betaCPL','betaONE','betaGDP','betaLRP','betaCON']
# 'betaCON','betaLRP',


directory = 'plots/'+modelName + '/' + imputeType + '/' +  averageOrPeriod + '/' + adultAge + '/'



lookup=pd.Series(['AxO','D','GDP','R','C','CH','U','P','GINI','C','W','E','T','S','S','CON','E','C','G','Gold','AxP'
                  ,'C',
                  'CFL','E','L','SEC','ENG','I','VIO','OUT','I','U','UD','PC3','PC4','L','L','N','C'],
            index=['AxO','DEM','GDP','RAT','COS','CHU','URB','POP','IQZ','LRP','WAR','ONE','TRU','SUP','CON','COP',
                   'PAR','COM','GIN','INS','AXP','AUT','CFL','EDS','LEX','SEC','ENG','INV','INN','OUT','CPL','CRI',
                   'CxD','PC3','PC4','IS1','IS2','NAT','COF'])


if averageOrPeriod == 'average':
    av=makeDfForPlot(modelName,imputeType,averageOrPeriod,adultAge)
    br=av
    
    
if averageOrPeriod == 'period':
    av=makeDfForPlot(modelName,imputeType,'average',adultAge)
    br=makeDfForPlot(modelName,imputeType,averageOrPeriod,adultAge)    

    
    


if targets is None:
    targets = sorted(av.loc[:,'targ'].unique())
    print(targets)

if preds is None:
    preds = sorted(av.loc[:,'pred'].unique())
    print(preds)



av1 = av
br1 = br

num = len(targets)+1

print(targets)


#if len(preds) == len(targets): f,ax = plt.subplots(len(targets)+1,len(preds)+1,figsize=[4,10])
#if len(preds) > len(targets): f,ax = plt.subplots(len(targets)+1,len(preds)+1,figsize=[3*num,2.2*num])#8.5

    
f,ax = plt.subplots(len(targets)+1,len(preds)+1,figsize=[18,4])    

row=1
for targ in targets:

    av3 = av1[av1.loc[:,'targ'] == targ]
    br2=br1[br1.loc[:,'targ'] == targ]

    col=1

    setRowColumnLabels(targ,row,0,True)

    for pred in preds:#sorted(av3.loc[:,'pred'].unique()):
        plotRow(av3,br2,pred,row,col,averageOrPeriod)

        setRowColumnLabels(pred[4:],0,col,False)


        col+=1

    row+=1



plt.tight_layout()
ax[0,0].set_visible(False)

#if not os.path.exists(directory):
 #   os.makedirs(directory)

plt.savefig('RegressionResults.pdf')
#plt.savefig(directory + 'RegressionResults.png')

#plt.close()


plt.show()






