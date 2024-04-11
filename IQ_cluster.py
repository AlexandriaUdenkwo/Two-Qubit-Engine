# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 16:54:54 2022
​
SVM+kMeans to classify the iq blobs
​
@author: DVK
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm
from sklearn.cluster import KMeans

def kMeans(train, n_clusters=2, flip=False, plot=False):
   
    x=[]; default=True
    for j in range(len(train[0])):
       
        x.append([train[0][j], train[1][j]])
       
    km = KMeans(n_clusters=n_clusters, init='random', n_init=100, max_iter=700, tol=1e-10, random_state=0)
    y_km = km.fit_predict(x)
   
   
    if flip:
        default=False
        for j in range(len(y_km)):
            if y_km[j] == 0:
                y_km[j] = -1
            elif y_km[j] == 1:
                y_km[j] = +1
               
    if default:
        for j in range(len(y_km)):
            if y_km[j] == 0:
                y_km[j] = +1
            elif y_km[j] == 1:
                y_km[j] = -1
    #plot
    if plot:
        x = np.array(x)
        plt.figure(figsize=(7, 7))
        plt.scatter(x[y_km == 1, 0], x[y_km == 1, 1], s=6, c='red', marker='.', label='|g>')
        plt.scatter(x[y_km == -1, 0], x[y_km == -1, 1], s=6, c='blue', marker='.', label='|e>')
        plt.legend()
   
    return y_km
   
def svm_class(train, y, test, C = 0.05, plot=True):
   
    x = []; x_out = []
   
    for j in range(len(train[0])):
       
        x.append([train[0][j], train[1][j]])
       
    for i in range(len(test[0])):  
       
        x_out.append([test[0][i], test[1][i]])
       
    clf = svm.LinearSVC(C=C, dual=False, tol=1e-20, max_iter=1000)
    
   
    clf.fit(x, y); #print(f"Fitting score = {clf.score(x,y):.3f}")
   
    if plot:
#        plt.figure(figsize=(7, 7))
#        x_line = np.linspace(np.min(train[0]), np.max(train[0]));  plt.title("Train")
#        plt.scatter(train[0], train[1], s=6, alpha=0.7); plt.plot(x_line, (-clf.coef_[0][0]/clf.coef_[0][1])*x_line - clf.intercept_[0] / clf.coef_[0][1], 'r--', label="LSVM")
#        plt.ylim(np.min(train[1])-200, np.max(train[1])+200); plt.legend(); plt.show()
       
        plt.figure(figsize=(7, 7))
        clf.fit(x_out, clf.predict(x_out))
        x_line_out = np.linspace(np.min(test[0]), np.max(test[0])); plt.title("Test")
        plt.scatter(test[0], test[1], s=6, alpha=0.7); 
        
        plt.plot(x_line_out, (-clf.coef_[0][0]/clf.coef_[0][1])*x_line_out - clf.intercept_[0] / clf.coef_[0][1], 'r--', label="LSVM1")
        plt.plot(x_line_out, (-clf.coef_[1][0]/clf.coef_[1][1])*x_line_out - clf.intercept_[1] / clf.coef_[1][1], 'b--', label="LSVM2")
        plt.plot(x_line_out, (-clf.coef_[2][0]/clf.coef_[2][1])*x_line_out - clf.intercept_[2] / clf.coef_[2][1], 'g--', label="LSVM3")
        plt.ylim(np.min(test[1]), np.max(test[1])); plt.legend(); plt.show()
        
        
        y1 = (-clf.coef_[0][0]/clf.coef_[0][1])*x_line_out - clf.intercept_[0] / clf.coef_[0][1]
        y2 = (-clf.coef_[1][0]/clf.coef_[1][1])*x_line_out - clf.intercept_[1] / clf.coef_[1][1]
        y3 = (-clf.coef_[2][0]/clf.coef_[2][1])*x_line_out - clf.intercept_[2] / clf.coef_[2][1]
        
        m1 = (-clf.coef_[0][0]/clf.coef_[0][1])
        m2 = (-clf.coef_[1][0]/clf.coef_[1][1])
        m3 = (-clf.coef_[2][0]/clf.coef_[2][1])
        b1 = - clf.intercept_[0] / clf.coef_[0][1]
        b2 = - clf.intercept_[1] / clf.coef_[1][1]
        b3 = - clf.intercept_[2] / clf.coef_[2][1]
        
        bound1 = (x_line_out + (b1/m1 + b2/m2)*0.5)*(2/(1/m1+1/m2))   #(y1+y2)/2 #need perp line
        bound2 = ((m2+m3)/2)*x_line_out + (b2+b3)/2 #(y2+y3)/2 #slightly off
        bound3 = ((m1+m3)/2)*x_line_out + (b1+b3)/2 #(y1+y3)/2 #seems to work
        
        Ax = [2/(b1/m1 + b2/m2),(m2+m3)/(b2+b3),(m1+m3)/(b1+b3)] # m/b
        By = [1/((b1/m2+b2/m2)*1/(1/m1+1/m2)) ,2/(b2+b3) ,2/(b1+b3)] #1/b
        
        plt.figure(figsize=(7, 7))
#        plt.scatter(test[0], test[1], s=6, alpha=0.7); 
        plt.plot(x_line_out[:], bound1[:], 'r--', label="LSVM1")
        plt.plot(x_line_out[:], bound2[:], 'b--', label="LSVM2")
        plt.plot(x_line_out[:], bound3[:], 'g--', label="LSVM3")
#        plt.plot(x_line_out[:], y2[:], 'y--', label="y2")
#        plt.plot(x_line_out[:], y3[:], '--', label="y3")
#        plt.ylim(np.min(test[1]), np.max(test[1])); plt.legend(); plt.show()
       
        x_out = np.array(x_out)
#        plt.figure(figsize=(7, 7))
        plt.scatter(x_out[clf.predict(x_out) == 1, 0], x_out[clf.predict(x_out) == 1, 1], s=6, c='black', marker='.', label='|g>')
        plt.scatter(x_out[clf.predict(x_out) == -1, 0], x_out[clf.predict(x_out) == -1, 1], s=6, c='red', marker='.', label='|e>')
        plt.scatter(x_out[clf.predict(x_out) == 0, 0], x_out[clf.predict(x_out) == 0, 1], s=6, c='blue', marker='.', label='|f>')
        plt.legend(); #plt.show()
        plt.ylim(np.min(test[1]), np.max(test[1])); plt.legend(); plt.show()
       
    return clf,clf.predict(x_out),clf.coef_,clf.intercept_,Ax,By #[-clf.coef_[0][0]/clf.coef_[0][1], clf.intercept_[0]/clf.coef_[0][1]]