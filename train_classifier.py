from sklearn.ensemble import RandomForestClassifier
import numpy as np
from numpy import genfromtxt
import sys
from sklearn import svm

def remove_nan(Y_test):
   nanIndex = np.argwhere(np.isnan(Y_test))
   Y_test[nanIndex] = 0
   return Y_test

#loc = "/Users/kavya/JHU/comp_bio/project/diabetes.tsv"
loc = "/Users/kavya/JHU/comp_bio/project/breast_carcinomaTEMP_2.tsv"
#loc = "/Users/kavya/JHU/comp_bio/project/cardiovascular_disease.tsv"
#loc = sys.argv[1]

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