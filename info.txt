To run Mc analysis for a catalog RUN : "python3 Mc_final.py"

programme will ask for inputs :-

bin_size = floating number between (0 to 10)
methods = one string value from ['mbass','maxc','gft','emr','mbs','clauset','all']
plot = boolean (True/False)

Method specific parameters:-
"mbass" : display all discontinuities = boolean (True/False) 
"clauset" : number of bootstrap samples = dtype(integer) ; p-value = floating number between (0, 1)

OUTPUT:-
1) Save best fit plot for specified method if plot==True.
2) print Mc and other output parameters for each method.
