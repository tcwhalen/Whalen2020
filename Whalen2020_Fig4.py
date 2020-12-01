## Tim C Whalen, last edited July 2020
# Plots Figure 4D-G from Whalen 2020 JNeurophys, predicting behavior and TH from SNr measures
# 
# Creates a "random forest" (only random in data used, but each tree uses all paremters) to
# predict TH (% dopamine remaining) or Behavior PC1 (a metric of behavioral dysfunction) from
# properties of in vivo SNr neurons. Also computes "permutation" importance (actually uses 
# every derangement of the test data, not a random set of permutations) for each predictor and
# compares several simpler forests to the full model.
# Note that results differ slightly from the published figure due to random seed differences.
#
# Requires python3 and the following packages: pandas, numpy, matplotlib and sklearn
# as well as Whalen2020_data_fig4.csv for importing data.
#
# TO RUN: in a terminal, run "python3 Whalen2020_Fig4.py predict", replacing predict 
# with either "th" or "behavior" (no quotes). With no argument, behavior is default.

import random
import sys
import pandas as pd  
import numpy as np  
from sklearn.tree import DecisionTreeRegressor
from sklearn import metrics
from math import sqrt

# parse argument
predict = 'Behavior PC1'
if len(sys.argv)>=2:
	arg = sys.argv[1].lower()
	if arg=='th':
		predict = 'TH'
	elif arg!='behavior':
		print('No valid argument given; predicting behavior')
else:
	print('No valid argument given; predicting behavior')

max_depth = None
min_leaf = 1 # leaves cannot be smaller than this
min_split = 3 # nodes smaller than this won't split (but leaves could be as small as min_leaf)
nMC = 1000 # number of trees
choose = 5 # samples in bootstrapped test set
alpha = .05 # for confidence intervals

data = pd.read_csv("data/Whalen2020_data_fig4.csv")
X = data.drop(["Behavior PC1","TH","Animal"], axis=1)  
y = data.loc[:,predict]
ncol = len(X.columns)
nrow = len(X.index)
nkeep = nrow - choose

# basic forest and derangements
trees = [DecisionTreeRegressor(random_state=i*121,max_depth=max_depth,min_samples_leaf=min_leaf,min_samples_split=min_split) for i in range(nMC)]
mse = [0 for i in range(nMC)]
rms = [0 for i in range(nMC)]
mse_perm_diff = [[0 for i in range(nMC)] for f in range(ncol)] # perm = permutation for importances

# special trees (to compare RMSE)
# names are just for internal bookkeeping
names = ['FR', 'CV', 'Burst/sec', 'Sync Frac', 'Osc Frac', 'All But Osc']
drops = [ \
["CV","Sync Frac","Burst/sec","Frac Osc"],\
["FR","Sync Frac","Burst/sec","Frac Osc"],\
["FR","CV","Sync Frac","Frac Osc"],\
["FR","CV","Burst/sec","Frac Osc"],\
["FR","CV","Sync Frac","Burst/sec"],\
['Frac Osc']] # which columns to drop for each special model
rms_int = [0 for i in range(nMC)] # int = intercept

nspecial = len(names)
trees_special = [[DecisionTreeRegressor(random_state=i*647,max_depth=max_depth,min_samples_leaf=min_leaf,min_samples_split=min_split) for i in range(nMC)] for m in range(nspecial)]
rms_special = [[0 for i in range(nMC)] for m in range(nspecial)]

for i in range(nMC):
	samp = random.sample(range(nrow),choose)
	X_train = X.drop(samp,axis=0)
	X_test = X.iloc[samp]
	y_train = y.drop(samp,axis=0)
	y_test = y.iloc[samp]
	trees[i].fit(X_train,y_train)
	y_pred = trees[i].predict(X_test)
	mse[i] = metrics.mean_squared_error(y_test, y_pred)
	rms[i] = sqrt(mse[i])
	intercep = trees[i].tree_.value[0]
	rms_int[i] = sqrt(np.sum([(intercep-val)**2 for val in y_test]))
	
	# derangements, computed in an easier order and stacked into a single frame
	for f in range(ncol):
		col = X.columns[f]
		X_test_perm = X.loc[[0]]
		X_test_perm = X_test_perm.drop([0],axis=0) # hack to get empty data frame w/ columns
		y_test_perm = y.loc[[0]]
		y_test_perm = y_test_perm.drop([0],axis=0)
		for j in range(choose): # for each other row, substitute this row's value for permuted feature
			X_temp_perm = X_test.copy()
			val = X_temp_perm.loc[samp[j],col]
			X_temp_perm = X_temp_perm.drop([samp[j]])
			X_temp_perm.loc[:,col] = [val for q in range(choose-1)]
			X_test_perm = X_test_perm.append(X_temp_perm)
			y_drop = y_test.drop(samp[j],axis=0)
			y_test_perm = y_test_perm.append(y_drop)
		y_pred_perm = trees[i].predict(X_test_perm)
		mse_perm_diff[f][i] = metrics.mean_squared_error(y_test_perm, y_pred_perm) - mse[i]
	
	# special models
	for m in range(nspecial):
		X_train_special = X_train.drop(drops[m],axis=1)
		X_test_special = X_test.drop(drops[m],axis=1)
		trees_special[m][i].fit(X_train_special,y_train)
		y_pred_special = trees_special[m][i].predict(X_test_special)
		rms_special[m][i] = sqrt(metrics.mean_squared_error(y_test, y_pred_special))
	
# feature importances
sort_perm_diff = [sorted(i) for i in mse_perm_diff]
fimp_conf = [[i[int(round(nMC*alpha-1))], i[int(round(nMC*(1-alpha)-1))]] for i in sort_perm_diff]
fimp_conf_rms = [[sqrt(abs(i))*np.sign(i) for i in j] for j in fimp_conf]
fimp = [np.median(i) for i in sort_perm_diff]
fimp_rms = [sqrt(abs(np.median(i)))*np.sign(np.median(i)) for i in sort_perm_diff]

# conf for full model
rms_sort = sorted(rms)
rms_conf = [rms_sort[i] for i in [int(round(nMC*alpha-1)), int(round(nMC*(1-alpha)-1))]]

# conf for intercept model - ended up only plotting median
rms_int_sort = sorted(rms_int)
rms_int_conf = [rms_int_sort[i] for i in [int(round(nMC*alpha-1)), int(round(nMC*(1-alpha)-1))]]


## plotting
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

colors = ['cyan','blue','red','green','yellow','purple']

# importances
plt.figure()
ca = plt.gca()
for f in range(ncol):
	ca.add_patch(Rectangle((f+0.3,fimp_conf_rms[f][0]),0.4,fimp_conf_rms[f][1]-fimp_conf_rms[f][0],fill=False))
	plt.plot([f+0.3, f+0.7],[fimp_rms[f],fimp_rms[f]],color=colors[f])

plt.plot([0, ncol],[0, 0],'k--',alpha=0.5)
plt.ylim(1.1*min([fimp_conf_rms[i][0] for i in range(ncol)]),1.1*max([fimp_conf_rms[i][1] for i in range(ncol)])) # only works if min is negative, max is positive
plt.xlim(0,ncol)
plt.xticks([.5+i for i in range(ncol)],X.columns)
plt.ylabel("RMSE Change")
plt.title("Feature Importance Predicting " + predict)
plt.show()

# special comparisons
plt.figure()
ca = plt.gca()
ca.add_patch(Rectangle((.3,rms_conf[0]),0.4,rms_conf[1]-rms_conf[0],fill=False))
plt.plot([0.3,0.7],[np.median(rms),np.median(rms)])
plt.plot([0.5],[np.median(rms_int)],marker='s',linestyle='None',color='black')
for m in range(nspecial):
	plt.plot([0.5],[np.median(rms_special[m])],marker='s',linestyle='None',color=colors[m])

plt.text(0.3,rms_conf[1]+.2,'full model (95%-conf)',color='black')
plt.ylabel("RMSE")
plt.title(predict + " Error")
plt.xticks([])
plt.xlim(.25,.75)
plt.show()