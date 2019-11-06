from sklearn.ensemble import RandomForestClassifier
import numpy as np
from numpy import genfromtxt

def remove_nan(Y_test):
   nanIndex = np.argwhere(np.isnan(Y_test))
   Y_test[nanIndex] = 0
   return Y_test

loc = "/Users/kavya/JHU/comp_bio/project/type2_diabetes.tsv"

my_data = genfromtxt(loc, delimiter='\t')
#print(my_data[1:-1, -1])
#my_data = my_data[1:-2, 1:-1]
np.random.shuffle(my_data)
X = my_data[1:-1, 1:-2]
#print(len(X))
X_train, X_test = np.split(X, 2)
Y = my_data[1:-1, -1]
Y_train, Y_test = np.split(Y, 2)

# this is sketchy
X_train = remove_nan(X_train)
Y_train = remove_nan(Y_train)
X_test = remove_nan(X_test)
Y_test = remove_nan(Y_test)
#print(Y_test)

classifier = RandomForestClassifier(n_estimators=10)
classifier = classifier.fit(X_train, Y_train)

print("OUR ACCURACY: " + str(classifier.score(X_test, Y_test)))