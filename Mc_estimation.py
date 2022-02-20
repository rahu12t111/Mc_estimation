#Python conversion of a R code to calculate Mc.
import scipy.stats as ss
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import seaborn as sns
from scipy.stats import norm
import math
import plotly.graph_objects as go
import plotly.figure_factory as ff
import subprocess
import statistics as stats


################## IMPORTANT:(data in file should be seperated by "|" and header of the column containing magnitude should be "magnitude") ###################
#Taking input from user.
file= input("Enter name of the file containing data : ")
country_name= input("Enter prefix to save figures : ")

#Extracting data from a csv file into a pandas dataframe.
data = pd.read_csv(file, sep='|')

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#FUNCTIONS

#0)function to bin the magnitude.
def binned_mag(mag,mbin):
	'''
	Input : take a list of magnitude and rational bin size as input
	Return : a list of binned magnitude. '''

	binned_mag=[]
	for i in mag:
		value=mbin*math.modf((mbin/2+i)/mbin)[1]
		binned_mag.append(value)
	return(binned_mag)


#1)#Function returning a list of unique_mag,counts,cumulative_count.
def fmd(mag,mbin):
	'''
	Create a frequency magnitude distribution for binned magnitude.
	Input : Magnitude list, binsize
	Output : Dictionary containing a list of unique values, their frequencies
			 and cumulative count list. '''

	res_unique = np.unique(mag, return_counts=True)

	#Creating unique list.
	unique_list=np.arange(min(res_unique[0]),max(res_unique[0])+mbin,mbin).tolist()

	#getting count for each element in unique list.
	count_list=[]

	for i in unique_list:
	    i=round(i,1)
	    rounded_unique_list=np.round_(res_unique[0],1)
	    if i in rounded_unique_list:
	        count_list.append(res_unique[1][(rounded_unique_list.tolist()).index(i)])
	    else:
	        count_list.append(0)

	#getting cumulative count list.
	reverse_count_list=count_list[::-1]
	rev_cumulative_count=np.cumsum(reverse_count_list)
	cumulative_count_list=(rev_cumulative_count[::-1]).tolist()

	res={'m':unique_list,'noncum':count_list,'cum':cumulative_count_list}
	return(res)

#2)bootstrap the list.
def bootstraplist(sample,nsample):
	'''
	Function to bootstrap samples from a list
	Input : sample -> list, nsample -> number of samples to be created.
	Output : returns a 2D array with each row as a bootstrapped sample. '''

	length=len(sample)
	value=np.random.choice(sample, size=(nsample, len(sample)))
	return(value)


#3)bootstrapping a sample from a probability distribution.
#Making bootstrap samples of a given frequency distribution for discrete observations.
#Function returns a 2D array with each row as a bootstrapped sample. 
def bootstrap_distribution(sample_obs,number_of_observations_per_sample,nsample,frequency_distribution_func):
	'''
	Function to return bootstrap sample from a frequency distribution.
	Input:-
	sample_obs : unique mag in sample.
	nsample : number of samples.
	frequency_distribution_function : returns frequency for a particular value of magnitude.

	Output:-
	returns a 2D array with nsample rows and each row containing a bootstrapped sample.'''

	length=len(sample_obs)
	frequency_list=[]
	for i in sample_obs:
		freq_i=frequency_distribution_func(i)
		frequency_list.append(freq_i)
	prob_list=[i/sum(frequency_list) for i in frequency_list]
	value=np.random.choice(sample_obs ,p=prob_list , size=(nsample, number_of_observations_per_sample))
	return(value)

#2)function to calculate value of Goodness of fit R.
def GFT_R(mag,mbin,FMD,Mc):
	''' Goodness parameter calculated using GFT. '''
	indmag=[]
	for j in range(len(mag)):
		if mag[j]>Mc-mbin/2:
			indmag.append(j)

	if indmag==[]:
		print('Dont have enough data')
		return()

	list1=[]
	[(list1.append(mag[i])) for i in indmag]
	var1=np.mean(list1)
	b=math.log10(math.exp(1))/((var1)-(Mc-mbin/2))
	a=math.log10(len(indmag))+b*Mc

	FMDcum_model=[]
	for j in FMD['m']:
		FMDcum_model.append(10**(a-b*j))

	indmi=[]
	for j in range(len(FMD['m'])):
		if FMD['m'][j]>=Mc:
			indmi.append(j)

	if indmi==[]:
		print('Dont have enough data')
		return()

	list2=[]
	list3=[]
	for j in indmi:
		list2.append(abs(FMD['cum'][j]-FMDcum_model[j]))
		list3.append(FMD['cum'][j])
	return((sum(list2)/sum(list3))*100)

# ---------------------------------------- PLOTTING -------------------------------------------
#3)function to plot GR_law.
def plot_GR(min_mag,max_mag,FMD,mag,mbin):
	'''
	Input:-
	min_mag : minimum magnitude in magnitude list.
	max_mag : maximum magnitude in magnitude list.
	FMD : frequency magnitude function calculated using "fmd" function.
	mag : list of magnitude.
	mbin : bin size

	Output:-
	returns a labled scatterplot

	NOTE:-Plots made using this function can be visualized by using plt.show() or saved using plt.savefig(path)
		  after applying this function, Don't forget to clear the already existing plots in your memory after saving
		  this. '''

	# here output will be log(cumulative count>=mag).
	def gr_law(mag,a,b):
	    return(a-b*(mag-min_mag))

	#fitting the data for magnitude (>Mc) and removing the magnitude with 0 events. 
	x_array=[FMD['m'][j] for j in range(len(FMD['m'])) if max_mag>=FMD['m'][j]>=min_mag and FMD['cum'][j]!=0]
	ind_x_array=[j for j in range(len(FMD['m'])) if FMD['m'][j] in x_array]
	y_array_test=[math.log(FMD['cum'][k]) for k in ind_x_array]
	param,correlation=curve_fit(gr_law,x_array,y_array_test)
	y_array=[gr_law(i,param[0],param[1]) for i in x_array]

	#Evaluating value of R.
	R=GFT_R(mag=mag,mbin=mbin,FMD=FMD,Mc=min_mag)

	#finding R^2 vatue for the fit.
	residuals=[(y_array_test[i]-y_array[i])**2 for i in range(len(x_array))]
	ss_res=sum(residuals)
	mean_diff_square=[(y_array_test[i]-np.mean(y_array_test))**2 for i in range(len(x_array))]
	ss_tot=sum(mean_diff_square)
	r_squared=1-(ss_res/ss_tot)

	#extrapolating the fit for all magnitudes.
	x_array1=[FMD['m'][j] for j in range(len(FMD['m'])) if FMD['cum'][j]!=0]
	y_array1=[gr_law(i,param[0],param[1]) for i in x_array1]

	#plotting the data on the graph. 
	plt.text(x=min(FMD['m'])+0.5, y=min(y_array1), s=r'$\log_{10}(N(m))$'+('={a:.2f}-{b:.2f}*(m-{Mc:.1f}) \n '+r'$r^2=$'+'{r_squared:.2f}').format(a=param[0],b=param[1],Mc=min_mag,r_squared=r_squared))
	plt.text(x=max(FMD['m'])*4/5, y=max(y_array1)*(8/10), s=('Mc={Mc:.1f},\nGFT R={R:.2f}').format(Mc=min_mag,R=R),fontweight='bold',bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10})
	plt.scatter(x=FMD['m'],y=[math.log(j) for j in FMD['cum']],s=1.5,color="black",alpha=1)
	plt.plot(x_array1,y_array1,'--',color='red')	 # extrapolated plot for complete magnitude range.
	#plt.plot(x_array,y_array,'--',color='red')      # plot of data of applicable magnitude range.
	plt.ylabel(r'$\log_{10}(Cumulative\ count)$',fontsize=12)
	plt.xlabel('Magnitude',fontsize=12)

#4)Function to plot cdf of normal distribution to our data.
def plot_EMR(xprocessed_data,yprocessed_data,res,FMD,mag,mbin):
	'''
	Input:-
	xprocessed_data : magnitude list correaponding to unique(mag).
	yprocessed_data : frequency for each magnitude, calculated using Ogata's distribution.
	FMD : frequency magnitude function calculated using "fmd" function.
	mag : list of magnitude.
	mbin : bin size

	Output:-
	returns a labled scatterplot
	
	NOTE:-Plots made using this function can be visualized by using plt.show() or saved using plt.savefig(path)
		  after applying this function, Don't forget to clear the already existing plots in your memory after saving
		  this. '''

	x_scatter=FMD['m']
	y_scatter=FMD['noncum']

	#finding R^2 vatue for the fit.
	residuals=[(yprocessed_data[i]-y_scatter[i])**2 for i in range(len(xprocessed_data))]
	ss_res=sum(residuals)
	mean_diff_square=[(y_scatter[i]-np.mean(y_scatter))**2 for i in range(len(x_scatter))]
	ss_tot=sum(mean_diff_square)
	r_squared=1-(ss_res/ss_tot)

	#Evaluating value of R.
	R=GFT_R(mag=mag,mbin=mbin,FMD=FMD,Mc=res['Mc'])

	#Plotting non-cumulative count.
	plt.text(x=max(FMD['m'])*(5.5/10), y=max(FMD['noncum'])-(max(FMD['noncum'])-min(FMD['noncum']))*(8/10), s=r'$N=Δ(10$'+('^({a:.2f}-{b:.2f}*(mag-{Mc:.1f})))').format(a=res['a'],b=res['b'],Mc=res['Mc'],r_squared=r_squared))
	plt.text(x=min(FMD['m']), y=min(FMD['noncum'])+(max(FMD['noncum'])-min(FMD['noncum']))*(1/5), s=('cdf.norm('+r'$μ=$'+'{mu:.2f}'+r'$,σ=$'+'{sigma:.2f})\n'+r'$r^2=$'+'{r_squared:.4f}').format(sigma=res['sigma'],mu=res['mu'],Mc=res['Mc'],r_squared=r_squared))
	plt.text(x=max(FMD['m'])*(7.5/10), y=max(FMD['noncum'])*(9.2/10), s=('Mc={Mc:.1f}').format(Mc=res['Mc']),fontweight='bold',bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10})
	plt.scatter(x=x_scatter,y=y_scatter,s=1.5,color="black",alpha=1)
	plt.plot(xprocessed_data,yprocessed_data,'--',color='red')
	plt.title('EMR non-cumulative distribution',fontsize=18, fontweight='bold')
	plt.ylabel('Count',fontsize=12)
	plt.xlabel('Magnitude',fontsize=12)
	plt.savefig(country_name+'EMR noncum')
	#plt.show()
	plt.clf()

	#Converting non-cumulative data to log(cumulative) data.
	x_scatter=FMD['m']
	y_scatter=[math.log10(i) for i in FMD['cum']]

	x_plot=xprocessed_data
	#calculating log cumulative count for processed data.
	y_plot=[]
	for i in range(len(yprocessed_data)):
		y=[j for j in yprocessed_data[i:]]
		y_plot.append(sum(y))

	#plotting the data on the graph. 
	plot_GR(min_mag=res['Mc'],max_mag=max(FMD['m']),FMD=FMD,mag=mag,mbin=mbin)
	plt.title('EMR',fontsize=18, fontweight='bold')
	plt.savefig(country_name+'EMR')
	#plt.show()
	plt.clf()

# --------------------------------------Functions for Mc Estimation-------------------------------------- #

#1) Maximum Curvature(MAXC)
def maxc(mag,mbin,plot=False):
	'''
	Input:-
	mag : binned magnitude list.
	mbin : binsize
	Output:-
	returns dictionary with vallue of Mc. '''

	FMD=fmd(mag,mbin)
	Mc=[]
	for i in range(len(FMD['noncum'])):
		if FMD['noncum'][i]==max(FMD['noncum']):
			Mc.append(round(FMD['m'][i],len(str(mbin).split('.')[1])))
	if plot:
		plot_GR(min_mag=Mc[0], max_mag=max(FMD['m']), FMD=FMD, mag=mag, mbin=mbin)
		plt.title('MAXC', fontsize=18, fontweight='bold')
		plt.savefig(country_name+'MAXC')
		#plt.show()
		plt.clf()
	#MAXC can give multiple magnitudes which have a same higheast number of events but we have to chose 1 of them and I chose the first one.
	return({'Mc':Mc[0]})     


#2) Goodness-of-fit test (GFT)
def gft(mag,mbin,plot=False):
	'''
	Input:-
	mag : list of binned magnitude
	mbins : bin size

	Output:-
	returns a dictionary containing Mc value, Confidence value('best'), R value . '''

	FMD=fmd(mag,mbin)
	McBound=maxc(mag,mbin)['Mc']  
	Mco=[]
	for i in np.arange(1,16,1):
		number_to_append=McBound-0.4+(i-1)/10
		Mco.append(number_to_append)

	R=[]
	[R.append(0) for i in range(15)]

	for i in range(15):
		indmag=[]
		for j in range(len(mag)):
			if mag[j]>Mco[i]-mbin/2:
				indmag.append(j)

		#Taking this to be the magnitude which can be observed.
		if len(indmag)==0:
			return({'Mc':Mco[i]})

		list1=[]
		[(list1.append(mag[i])) for i in indmag]
		var1=np.mean(list1)
		b=math.log10(math.exp(1))/((var1)-(Mco[i]-mbin/2))
		a=math.log10(len(indmag))+b*Mco[i]

		FMDcum_model=[]
		for j in FMD['m']:
			FMDcum_model.append(10**(a-b*j))

		indmi=[]
		for j in range(len(FMD['m'])):
			if FMD['m'][j]>=Mco[i]:
				indmi.append(j)

		if indmi==[]:
			return({'Mc':Mco[i]})

		list2=[]
		list3=[]
		for j in indmi:
			list2.append(abs(FMD['cum'][j]-FMDcum_model[j]))
			list3.append(FMD['cum'][j])
		R[i]=(sum(list2)/sum(list3))*100 


	indGFT=[j for j in range(len(R)) if R[j]<=5]   #95% confidence
	if len(indGFT)!=0:
		Mc=Mco[indGFT[0]]
		best='95%'
	else:
		indGFT=[j for j in range(len(R)) if R[j]<=10]   #90% confidence
		if len(indGFT)!=0:
			Mc=Mco[indGFT[0]]
			best='90%'
		else:
			Mc=McBound
			best='MAXC'

	if plot:
		plot_GR(min_mag=Mc,max_mag=max(FMD['m']),FMD=FMD,mag=mag,mbin=mbin)
		plt.title('GFT',fontsize=18, fontweight='bold')
		plt.savefig(country_name+'GFT')
		#plt.show()
		plt.clf()

	return({'Mc':Mc,'best':best,'Mco':Mco,'R':R})



#3)Mc by b-val Stability (MBS)
def mbs(mag,mbin,cutoff=0.03,plot=False):
	'''
	Input:-
	mag : list of binned magnitudes,
	mbin : binsize

	Output:-
	returns a dictionary containing value of Mc ans other parameters. '''

	McBound=maxc(mag,mbin)['Mc']
	Mco=[]
	for i in np.arange(1,21,1):
		number_to_append=McBound-0.7+(i-1)/10
		Mco.append(number_to_append)

	bi=[]
	unc=[]
	for i in range(20):
		indmag=[j for j in range(len(mag)) if mag[j]>Mco[i]-mbin/2]
		nbev=len(indmag)

		#Taking this to be the magnitude which can be observed.
		if nbev<2:
			return({'Mc':Mco[i]})

		list1=[]
		[(list1.append(mag[i])) for i in indmag]
		var1=np.mean(list1)
		bi.append(math.log10(math.exp(1))/(var1 -(Mco[i]-mbin/2)))
		list2=[]
		for j in list1:
			diff_square=(i-var1)**2
			list2.append(diff_square)
		unc.append(2.3*bi[i]**(2*(math.sqrt(sum(list2))/(nbev*(nbev-1)))))

	bave=[];[bave.append(np.mean(bi[i:i+6])) for i in range(15)]

	dbi_old=[];[dbi_old.append(abs(bi[i+1]-bi[i])) for i in range(len(bi)-1)]

	indMBS_old=[i for i in range(len(dbi_old)) if dbi_old[i]<=cutoff]
	dbi=[abs(bave[i]-bi[i]) for i in range(15)]
	indMBS=[i for i in range(15) if dbi[i]<=unc[i]]
	Mc=Mco[indMBS[0]]

	#If all data is complete and follows a gutenberg-richter law then value of mc predicted by the above algorithm can be less than min value of magnitude in our catalogue.
	minimum_magnitude=min(mag)
	if Mc<minimum_magnitude:
		Mc=minimum_magnitude

	if plot:
		FMD=fmd(mag,mbin)
		plot_GR(min_mag=Mc,max_mag=max(FMD['m']),FMD=FMD,mag=mag,mbin=mbin)
		plt.title('MBS',fontsize=18, fontweight='bold')
		plt.savefig(country_name+'MBS')
		#plt.show()
		plt.clf()

	return({'Mc':Mc,'Mco':Mco,'bi':bi,'unc':unc,'bave':bave})


#4)Entire Magnitude Range method (EMR)
def emr(mag,mbin,plot=False,add={'Mc':50}):
	'''
	Input:-
	mag : list of binned magnitudes.
	mbin : binsize (rational number).
	plot : boolean to specify plotting.
	add : Do not change this parameter.

	Output:
	returns a dictionary containing value for Mc, parameters a and b in G-R law,
	parameters of norm CDF, output of the model, Probablity of getting the obse-
	rvations from the model. '''

	FMD=fmd(mag,mbin)
	nbm=len(FMD['m'])

	McMAXC=maxc(mag,mbin)['Mc']
	mu=abs(McMAXC/2); sig=abs(McMAXC/4)
	if mu>1:
		mu=abs(McMAXC/10); sig=abs(McMAXC/20)
	McBound=McMAXC
	Mco=[]; [Mco.append(round(McBound-0.3+((i-1)/10),1)) for i in np.arange(1,10,1)]
	params=[];[params.append([0,0,0,0]) for i in range(9)]    #here each column corresponds to a,b,mu,sigma.
	list_len_nbm=[]; [list_len_nbm.append(0) for i in range(nbm)]
	savedmodel=[]; [savedmodel.append(list_len_nbm) for j in range(9)]
	prob=[]; [prob.append(list_len_nbm) for i in range(9)]
	for i in range(9):
		indmag=[j for j in range(len(mag)) if mag[j]>Mco[i]-mbin/2]
		nbev=len(indmag)

		#Taking this to be the magnitude which can be observed.
		if nbev<2:
			return({'Mc':Mco[i]})

		list1=[]
		[(list1.append(mag[i])) for i in indmag]
		var1=np.mean(list1)
		b=math.log10(math.exp(1))/((var1)-(Mco[i]-mbin/2))
		a=math.log10(len(indmag))+b*Mco[i]

		cum_N=[]
		for j in FMD['m']:
			cum_N.append(10**(a-b*j))

		params[i][0]=a ; params[i][1]=b

		cumNtmp=10**(a-b*(max(FMD['m'])+mbin))
		cumNtmp=cum_N[:]+[cumNtmp]

		N=[]
		for j in range(len(cumNtmp)-1):
			abs_diff=abs(cumNtmp[j]-cumNtmp[j+1])
			N.append(abs_diff)

		data=pd.DataFrame({'N':N,'m':FMD['m'],'Nd':FMD['noncum']})
		indLow=[j for j in range(len(FMD['m'])) if FMD['m'][j]<Mco[i]]
		indHigh=[j for j in range(len(FMD['m'])) if FMD['m'][j]>=Mco[i]]

		dataTest=pd.DataFrame({'N':[data['N'][j] for j in indLow],'m':[data['m'][j] for j in indLow],'Nd':[data['Nd'][j] for j in indLow]})
		dataTmp=pd.DataFrame({'N':[data['N'][j] for j in indHigh],'m':[data['m'][j] for j in indHigh],'Nd':[data['Nd'][j] for j in indHigh]})

		#We will have to apply different process if data for all magnitudes is complete.
		if len(indLow)<=2:
			checkNo0=[]
		else:
			checkNo0=[j for j in range(len(dataTest['Nd'])) if dataTest['Nd'][j]!=0]
			if len(checkNo0)<=2:
				None
			else:
				dataTest=pd.DataFrame({'N':[data['N'][j] for j in checkNo0],'m':[data['m'][j] for j in checkNo0],'Nd':[data['Nd'][j] for j in checkNo0]})

				#Nmax=max(dataTmp['Nd'])
				Nmax=max(dataTest['Nd'])
				#Nmax=dataTmp['Nd'][len(max(dataTmp['Nd']))]
				Mmintmp=min(dataTest['m'])

				list1=[]
				for j in range(len(dataTest['Nd'])):
					list1.append(dataTest['Nd'][j]/Nmax)
				dataTest['Nd']=list1

				list2=[]; [list2.append(dataTest['m'][j]-Mmintmp) for j in range(len(dataTest['m']))]
				dataTest['m']=list2

				data4fit=pd.DataFrame({'N':dataTest['Nd'],'m':dataTest['m']})

				def cdf_norm(x,mean,sd,):
					return(ss.norm(loc=mean,scale=sd).cdf(x))

				curvefit_param,correlation=curve_fit(cdf_norm,data4fit['m'],data4fit['N'],p0=[mu,sig])
				params[i][2]=curvefit_param[0]; params[i][3]=curvefit_param[1]

				dataTest['N']=[(ss.norm(loc=curvefit_param[0],scale=curvefit_param[1]).cdf(j))*Nmax for j in dataTest['m']]

				list3=[]
				[list3.append(j+Mmintmp) for j in dataTest['m']]
				dataTest['m']=list3

				list4=[]
				[list4.append(j*Nmax) for j in dataTest['Nd']]
				dataTest['Nd']=list4

		dataPred=pd.DataFrame({'N':(dataTest['N'].tolist()+dataTmp['N'].tolist()),'m':(dataTest['m'].tolist()+dataTmp['m'].tolist()),'Nd':(dataTest['Nd'].tolist()+dataTmp['Nd'].tolist())})
		
		if plot and float(Mco[i])==float(add['Mc']):
			FMD=fmd(mag,mbin)
			plot_EMR(xprocessed_data=dataPred['m'],yprocessed_data=dataPred['N'],res=add,FMD=FMD,mag=mag,mbin=mbin)


		list5=[]
		[list5.append(round(j)) for j in dataPred['N']]
		dataPred['N']=list5

		k=0
		for j in (checkNo0+indHigh):
			savedmodel[i][j]=dataPred['N'][k]
			k=k+1


		#Logarithm to the basis of 10 of Poisson probability density.
		probtmp=[]; [probtmp.append(0) for j in range(nbm)]
		CheckNo0=[j for j in range(len(dataPred['N'])) if dataPred['N'][j]!=0]
		Pmodel=[dataPred['N'][j] for j in CheckNo0]
		Pdata=[dataPred['Nd'][j] for j in CheckNo0]

		k=0
		for j in CheckNo0:
			probtmp[j]=(1/math.log(10))*(-Pmodel[k]+Pdata[k]*math.log(Pmodel[k])-math.lgamma(Pdata[k]+1))

		prob[i]=-(sum(probtmp))

	indbestfit=[j for j in range(len(prob)) if prob[j]==min(prob)][0]

	res={'Mc':Mco[indbestfit],'a':params[indbestfit][0],'b':params[indbestfit][1],'mu':params[indbestfit][2],'sigma':params[indbestfit][3],'model':savedmodel[indbestfit],'Mco':Mco,'prob':prob}

	if plot and add['Mc']!=res['Mc']:
		emr(mag=mag,mbin=mbin,plot=True,add=res)

	return(res)


#5)Median-based analysis of the segment slope (MBASS)
def fmbass(a,delta=0.1,plot=False,alldisc=False):
	'''
	This function determine the number of change points in our data.
	Input:-
	a = binned magnitude list.
	delta = bin size.
	Output:-
	returns the value of change poins, unique magnitudes and cumulative counts 
	for those unique magnitude to mbass function defined in this module. '''

	a=list(a)
	tau=[]
	pva=[]
	minmag=min(a)

	#generating data which will be used in the function.
	copy_list_of_number=a[:]
	x=[]
	final_list=a

	for i in final_list:
	    if i not in x:
	        x.append(i)
	x.sort()

	count_list=[copy_list_of_number.count(i) for i in x]

	cumulative_count_list=[]

	final_count=len(a)

	count=0
	for i in count_list:
	    end_count=final_count-count
	    cumulative_count_list.append(end_count)
	    count=count+i

	ind0=[i for i in range(len(x)) if cumulative_count_list[i]!=0]
	x=[x[:][i] for i in ind0]
	cumul=[cumulative_count_list[i] for i in ind0]
	log_n=[math.log10(cumulative_count_list[i]) for i in ind0]


	#defining derivative at all points in x[1:len(x)].
	sl=[(log_n[i+1]-log_n[i])/(x[i+1]-x[i]) for i in range(len(x)-1)]
	xsl=x[1:]

	#defining the initial variables for while loop ahed.
	niter=4
	N=len(sl)
	j=0    #number of iterations.
	k=0    #number of discontinuities.
	SA=[]
	[SA.append('') for i in range(N)]

	#starting the while loop.
	while j<niter:
		for i in np.arange(0,N,1):
			SA[i]=abs(2*sum(ss.rankdata(sl)[0:i])-(i+1)*(N+1))

		n1=[]
		for i in range(len(SA)):
			if SA[i]==max(SA):
				n1.append(i)

		xn1=sl[0:n1[0]]
		xn2=sl[n1[0]:]

		if n1[0]>2 and n1[0]<=(N-2) and ss.ranksums(xn1,xn2)[1]<0.05:
			k=k+1
			[pva.append('') for i in range(k-len(pva))]; pva[k-1]=ss.ranksums(xn1,xn2)[1]
			[tau.append('') for i in range(k-len(tau))]; tau[k-1]=n1[0]

			if k>1:
				meds11=np.median(sl[0:n0])
				meds12=np.median(sl[n0:])
				for i in np.arange(0,n0,1):
					sl[i]=sl[i]+meds11
				for i in np.arange(n0,len(sl),1):
					sl[i]=sl[i]+meds12

			meds11=np.median(sl[0:n1[0]])
			meds12=np.median(sl[n1[0]:])
			for i in np.arange(0,n1[0],1):
				sl[i]=sl[i]-meds11
			for i in np.arange(n1[0],len(sl),1):
				sl[i]=sl[i]-meds12
			n0=n1[0]
		# if there is no discontinuity or it is present near the edges then we select min(mag) as first point of change.
		else:
			pva=[0]
			tau=[0]

		j=j+1

	#we have to define a function equivalent to function order in R.
	def order(a):
		list_1=a[:]
		order_list=[]
		for i in range(len(a)):
			index_min=list_1.index(min(list_1))
			order_list.append(index_min)
			list_1.remove(min(list_1))
			list_1.insert(index_min,100)
		return(order_list)

	ip=order(pva)

	#if there is only one discontinuity then I am taking the interval between xsl[tau[ip[0]]] and max(mag) for fitting the GR law.
	if len(ip)>=2:
		m0=[round(xsl[tau[ip[0]]],5),round(xsl[tau[ip[1]]],5)]
	else:
		m0=[round(xsl[tau[ip[0]]],5),max(xsl)]

	Number_discontinuities=len(ip)

	if alldisc:
		print({"discmag":[xsl[i] for i in tau],"p":pva,"m0":m0,"discmag_no":Number_discontinuities})
	return(m0,x,cumul)

def mbass(a,delta=0.1,plot=False,alldisc=False,bs=0):
	'''
	Input:-
	a : binned magnitude list.
	delta : bin size (rational number).
	plot : boolean to specify plotting.
	alldisc : boolean to find all discontinuities.
	bs : number of bootstrap samples to be used.

	Output:-
	returns the value of change poins, unique magnitudes and cumulative counts 
	for those unique magnitude to mbass function defined in this module. '''

	def mba(x):
		return(fmbass(x,delta,plot,alldisc))
	if bs==0:
		res=mba(a)    #actual FMD analyzed.
	else:
		res=ss.bootstrap((a,),statistic=mba,vectorized=False,n_resamples=abs(bs))

	if plot:
		FMD={'m':res[1],'cum':res[2]}

		if res[0][0]==max(res[0]):
			max_mag=max(FMD['m'])
		else:
			max_mag=float(res[0][1])

		plot_GR(min_mag=float(res[0][0]),max_mag=max_mag,FMD=FMD,mag=a,mbin=delta)
		plt.title('MBASS',fontsize=18, fontweight='bold')
		plt.savefig(country_name+'MBASS')
		#plt.show()
		plt.clf()

	return(res)

#6)Clauset test (Clauset)

#METHOD

#1)Tinti_aval function.
def Tinti_aval(mag,M0,mbin,bs=False):
	'''
	Function used to calculate parameters of best fit for G-R law given the data.
	Input:-
	mag : binned magnitude list.
	M0 : estimate for magnitude of completeness.
	mbins : bin size.
	bs : number of bootstrap samples used to calculate G-R law parameters.
	
	Output:-
	returns a dictionary of parameters and log likelyhood of getting the distribution
	from the model we obtained. '''

	M0=round(M0,1)
	mag=[i for i in mag if i>=M0]
	mean_mag=(sum([i-M0 for i in mag]))/len(mag)
	beta=(1/mbin)*math.log(1+(mbin/mean_mag))
	b=beta/math.log(10)     #estimate for value of b.
	FMD=fmd(mag,mbin)
	to_calc_a=[]
	for i in FMD['m']:
	    value=10**(-b*(i-mbin/2))
	    to_calc_a.append(value)
	aval=math.log10(len(mag)/sum(to_calc_a)/(1-10**(-b*mbin)))   #estimate for value of a.

	#Calculating log likelyhood.
	LL=0
	for i in mag:
	    LL=LL+(math.log(1-math.exp(-beta*mbin)) - beta*(i-M0))  #estimating log_likelyhood.

	res={'beta':beta,'aval':aval,'b':b,'LL':LL}     #reaults of above calculations.

	##Bootstrapping for error.
	if bs==True:
	    Nbs=1000
	    bootstrap_array=bootstraplist(mag,Nbs)
	    betabs=[]
	    avalabs=[]
	    for i in range(len(bootstrap_array)):
	        params=Tinti_aval(bootstrap_array[i])
	        betabs.append(params['beta'])
	        avalabs.append(params['aval'])
	    res['betaquants']=[np.quantile(betabs,0.025),np.quantile(betabs,0.975)]
	    res['avalquants']=[np.quantile(avalabs,0.025),np.quantile(avalabs,0.975)]
	    return(res)

	else:
	    return(res)

#2)ClausetTest function.
def ClausetTest(orgmag,mbin,nsim,pval,plot=False):
	def ClausetTestFunc(orgmag,Mc,mbin,Mmax,nsim,pval,fun=False,a=None,b=None,plot=False):

		'''
		orgmag = binned magnitude list
		a:  avalue (can be left empty)
		b: b value (can be left empty)
		Mc: magnitude of completeness
		mbin: magnitude bin size
		Mmax: maximum magnitude
		nsim: number of simulations (reasonable value : >=1000)
		pval: minimum fraction of synthetic cases with ks distance larger in the
		synthetic catalogs than in the real catalog

		outputs:
		stat:fraction of synthetic cases with ks distance larger in the
		synthetic catalogs than in the real catalog
		h :  1 or 0 depending on whether GR law can be rejected or accepted
		dobs: ks distance between the emprical and best fit GR law. '''

		Mc=round(Mc,1)
		# t1=datetime.datetime.now()
		orgmag=[i for i in orgmag if i>=Mc]
		FMD=fmd(orgmag,mbin)
		FMD_copy=FMD.copy()

		# when we dont have enough earthquakes to analize.
		if len(FMD['m'])<2 and fun==True:
			return(0)
		elif len(FMD['m'])<2 and fun==False:
			print('not enough data')
			return({'Mc':FMD['m'][0]})

		# t2=datetime.datetime.now()
		# print(t2-t1)
		orgcdf=[i/max(FMD['cum']) for i in FMD['cum']]

		if a==None or b==None:
		    params=Tinti_aval(orgmag,Mc,mbin,bs=False)
		    a=params['aval'] ; b=params['b']

		Norgcdf=[i/max(FMD['cum']) for i in FMD['cum']]

		#Calculating theoritical frequency distribution and cdf.
		theocount=[(10**(a-b*(i-mbin/2))-10**(a-b*(i+mbin/2))) for i in FMD['m']]
		theocount.reverse()
		theocdf=np.cumsum(theocount)/sum(theocount)
		theocdf=theocdf.tolist()
		theocdf.reverse()

		#change from mat.???????????????????
		# theo_cum_count=[10**(a-b*(i-mbin/2)) for i in orgmag]
		# theocdf=[i/max(theo_cum_count) for i in theo_cum_count]

		# Function for plotting the graph.
		if plot==True:
		    plt.clf()
		    index_of_nonzeros=[i for i in range(len(Norgcdf)) if Norgcdf[i]!=0]
		    mag_array=[FMD['m'][i] for i in index_of_nonzeros]
		    non_zero_norgcdf=[Norgcdf[i] for i in index_of_nonzeros]
		    plt.scatter(mag_array,[math.log10(i) for i in non_zero_norgcdf],color='red',s=2)
		    plt.plot(FMD['m'],[math.log10(i) for i in theocdf],color='black')
		    plt.text(x=max(FMD['m'])*4/5, y=-1, s=('Mc={Mc:.1f}').format(Mc=Mc),fontweight='bold',bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10})
		    plt.title('Cumulative Distribution Mc={0}'.format(round(Mc,1)),fontsize=18, fontweight='bold')
		    plt.ylabel('log normalized Cumulative Count',fontsize=12)
		    plt.xlabel('Magnitude',fontsize=12)

		if fun==True:

		    dsim=max([abs(Norgcdf[i]-theocdf[i]) for i in range(len(Norgcdf))])

		    Nevt=len(orgmag)

		    return(dsim)

		if fun==False:

		    dobs=max([abs(Norgcdf[i]-theocdf[i]) for i in range(len(Norgcdf))])

		    Nevt=len(orgmag)

		    def GR_fmd(x,a=a,b=b,mbin=mbin):
		        return(10**(a-b*(x-mbin/2))-10**(a-b*(x+mbin/2)))

		    theo_bootstrap=bootstrap_distribution(FMD['m'],len(orgmag),nsim,GR_fmd)

		    dsim=[]
		    for i in range(len(theo_bootstrap)):
		        FMD=fmd(theo_bootstrap[i],mbin)
		        dsim_value=ClausetTestFunc(theo_bootstrap[i],Mc,mbin,Mmax,nsim,pval,fun=True,a=None,b=None,plot=False)
		        dsim.append(dsim_value)

		    stat=sum([1 for i in dsim if i>dobs])/nsim
		    h=int(stat<pval)

		res={'stat':stat,'h':h,'dobs':dobs,'pval':pval}

		#Save the graph only if the alternate hypothesis is True.
		if plot==True and res['h']==0:
		    plt.savefig(country_name+'CLAUSET')
		    plt.clf()

		return({'FMD':FMD_copy,'mag':orgmag,'mbin':mbin,'stats':[res,params]})

	#Part of main function.
	FMD=fmd(orgmag,0.1)
	for i in FMD['m'][:-5]:
		test_result=ClausetTestFunc(a,i,mbin,max(a),nsim,pval,plot=plot)
		if test_result['stats'][0]['h']==0:
			plot_GR(i,max(test_result['FMD']['m']),test_result['FMD'],test_result['mag'],test_result['mbin'])
			plt.title('MBASS',fontsize=18, fontweight='bold')
			plt.savefig(country_name+'clauset_GR')
			#plt.show()
			plt.clf()

			# Returning the output.
			return({'Mc':i}.update(test_result))
			break

	return('\n\ndistribution does not follow GR_law')





##DEFINING THE PARAMETERS:-

'''
1)a=list of magnitude.
2)delta=mbin=bin size
3)plot=takes boolean value, 
	 if it is set True:plot will be shown and saved; 
	 else no plot will be drawn.

In mbass we have extra parametres.
4)alldisc= boolean value
		 if True then values for all the discontinuities will be shown;
		 else only value of Mc will be returned.
5)bs=number of bootstrap samples to be taken.
'''

#EXAMPLES.
#Remove the comment(#) from the examples to try them.

#Extracting magnitude data("which should be complete and should not contain any blank fields") in form of a list.
bin_size = float(input('Enter bin size :'))
a = binned_mag(data['magnitude'].tolist(),bin_size)
#User Input.
method = input('Method used to calculate Mc :')
plot_bool = bool(input('Save plots (True/False) :'))

if method == 'mbass':
	print('INPUT')
	all_discontinuities = bool(input('Return all discontinuities in MBASS (True/False) :'))
	x=mbass(a,delta=bin_size,plot=plot_bool,alldisc=True)
	print('mbass\n',x)

if method == 'maxc':
	x=maxc(a,bin_size,plot=plot_bool)
	print('maxc\n',x)

if method == 'gft':
	x=gft(a,bin_size,plot=plot_bool)
	print('gft\n',x)

if method == 'mbs':
	x=mbs(a,bin_size,plot=plot_bool)
	print('mbs\n',x)

if method == 'emr':
	x=emr(a,bin_size,plot=plot_bool)
	print('emr\n',x)

if method == 'clauset':
	print('INPUT')
	bs = int(input('Number of bootstrap samples for clauset test :'))
	p_value = float(input('P value for hypothesis testing in clauset test [0 to 1]:'))
	x=ClausetTest(a,bin_size,bs,p_value,plot=plot_bool)
	print('clauset\n',x)

if method == 'all':
	print('INPUT')
	all_discontinuities = bool(input('Return all discontinuities in MBASS (True/False) :'))
	bs = int(input('Number of bootstrap samples for clauset test :'))
	p_value = float(input('P value for hypothesis testing in clauset test [0 to 1]:'))

	x=mbass(a,delta=bin_size,plot=plot_bool,alldisc=True)
	print('mbass\n',x)
	x=maxc(a,bin_size,plot=plot_bool)
	print('maxc\n',x)
	x=gft(a,bin_size,plot=plot_bool)
	print('gft\n',x)
	x=mbs(a,bin_size,plot=plot_bool)
	print('mbs\n',x)
	x=emr(a,bin_size,plot=plot_bool)
	print('emr\n',x)
	x=ClausetTest(a,bin_size,bs,p_value,plot=plot_bool)
	print('clauset\n',x)


if method not in ['mbass','maxc','gft','emr','mbs','clauset','all']:
	print('\nERROR: Give correct input, refer \"Readme.txt\" file for correct inputs.')

print('----------- COMPLETED ------------')
