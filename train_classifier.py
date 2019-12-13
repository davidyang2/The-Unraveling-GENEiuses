from sklearn.ensemble import RandomForestClassifier
import numpy as np
from numpy import genfromtxt
import sys
from sklearn import svm
from sklearn.ensemble import GradientBoostingClassifier

import gwas_feature_library
import matplotlib.pyplot as plt
import pandas

def remove_nan(Y_test):
   nanIndex = np.argwhere(np.isnan(Y_test))
   Y_test[nanIndex] = 0
   return Y_test

def plot_features(loc):
   file = pandas.read_csv(loc, sep='\t', lineterminator='\r')
   df = pandas.DataFrame(file)
   features_to_plot = ['CONTEXT', 'INTERGENIC', 'CHR_ID', 'CHR_POS', 'UP_DIST', 'DOWN_DIST', 'GO FUNCTION', 'TISSUE_EXPRESSION LEVEL']
   is_risk_factor = list(df["IS_RISK_FACTOR"].values)
   
   # get positive and negative division points
   split_point = -1
   for i in range(len(is_risk_factor)):
      if is_risk_factor[i] == 0:
         split_point = i
         break

   print(len(is_risk_factor))
   print(split_point)

   for i in range(len(features_to_plot)):
      feature_name = features_to_plot[i]
      feature = list(df[feature_name].values)
      print("Plotting " + str(feature_name))

      # get positive and negative features
      pos_features = feature[0:split_point]
      neg_features = feature[split_point:]

      # continuous features
      if feature_name == "UP_DIST" or feature_name == "DOWN_DIST" or feature_name == "CHR_POS":
         fig, ax = plt.subplots()
         bins = np.linspace(-10, 10, 30)
         colors = ['b', 'g']
         ax.hist([pos_features, neg_features], color=colors)
         ax.set_ylabel("Count")
         #plt.hist(neg_counts, bins=10, alpha=0.5, label='Negative')
         #plt.legend(loc='upper right')
         fig.savefig(feature_name + "_histograms.pdf")

      else:
         # categorical features
         pos_counts = gwas_feature_library.get_feature_counts_for_plotting(pos_features)
         neg_counts = gwas_feature_library.get_feature_counts_for_plotting(neg_features)
         gwas_feature_library.plot_feature(pos_counts, neg_counts, feature_name)

#loc = "/Users/kavya/JHU/comp_bio/project/The-Unraveling-GENEiuses/diabetes.tsv"
loc = "/Users/kavya/JHU/comp_bio/project/The-Unraveling-GENEiuses/breast_carcinomaTEMP.tsv"
#loc = "/Users/kavya/JHU/comp_bio/project/The-Unraveling-GENEiuses/cardiovascular_disease.tsv"
#loc = sys.argv[1]

plot_features(loc)

my_data = genfromtxt(loc, delimiter='\t')
#print(my_data[1:-1, -1])
#my_data = my_data[1:-2, 1:-1]
np.random.shuffle(my_data)
#print(my_data)

try:
   X = my_data[1:-1, 1:-2]
   Y = my_data[1:-1, -1]
   X_train, X_test = np.split(X, 2)
   Y_train, Y_test = np.split(Y, 2)
except ValueError:
   X = my_data[1:-2, 1:-2]
   Y = my_data[1:-2, -1]
   X_train, X_test = np.split(X, 2)
   Y_train, Y_test = np.split(Y, 2)

print(X_test)
# this is sketchy
X_train = remove_nan(X_train)
Y_train = remove_nan(Y_train)
X_test = remove_nan(X_test)
Y_test = remove_nan(Y_test)
print(Y_test)

classifier = RandomForestClassifier(n_estimators=10)
classifier = classifier.fit(X_train, Y_train)

print("RF ACCURACY: " + str(classifier.score(X_test, Y_test)))

classifier2 = svm.SVC()
classifier2 = classifier2.fit(X_train, Y_train)

print("SVM ACCURACY: " + str(classifier2.score(X_test, Y_test)))

clf = GradientBoostingClassifier(n_estimators=100, learning_rate=1.0, \
   max_depth=1, random_state=0).fit(X_train, Y_train)
clf.score(X_test, Y_test)
print("GBRT ACCURACY: " + str(clf.score(X_test, Y_test)))
