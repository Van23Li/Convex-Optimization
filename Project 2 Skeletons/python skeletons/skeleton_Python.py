#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
import time


# # The case for Survival Analysis

# In[ ]:


def mean_absolute_error(y_true, y_pred):
    return np.abs(y_true-y_pred).mean()


# In[ ]:


# load the data
X_train = np.load('X_train.npy')
y_train = np.load('y_train.npy')
z_train = np.load('z_train.npy')

X_test = np.load('X_test.npy')
y_test = np.load('y_test.npy')
z_test = np.load('z_test.npy')


# In[ ]:


plt.hist(y_train[z_train==0], bins = 20, alpha=0.5, label='Not Churned')
plt.hist(y_train[z_train==1], bins = 20, alpha=0.5, label='Churned')
plt.legend()
plt.title('Histogram of Survival Times');


# In[ ]:


# metrics for evaluation
def plot_results(y, prediction, z, save=False):
    plt.figure(figsize=(20,8))
    plt.subplot(131)
    plt.scatter(prediction, y, s=1)
    plt.xlabel('predicted')
    plt.ylabel('actual')
    plt.plot([0,70],[0,70], color='C1')
    plt.title(f'MAE {mean_absolute_error(y, prediction)}')

    plt.subplot(132)
    plt.scatter(prediction[z==1], y[z==1], s=1)
    plt.xlabel('predicted')
    plt.ylabel('actual')
    plt.plot([0,70],[0,70], color='C1')
    plt.title(f'MAE (Churned) {mean_absolute_error(y[z==1], prediction[z==1])}')

    plt.subplot(133)
    plt.scatter(prediction[z==0], y[z==0], s=1)
    plt.xlabel('predicted')
    plt.ylabel('actual')
    plt.plot([0,70],[0,70], color='C1')
    plt.title(f'MAE (Not Churned) {mean_absolute_error(y[z==0], prediction[z==0])}');
    if save:
        plt.tight_layout()
        timestamp = time.time()
        plt.savefig(f'comparison_{timestamp}.pdf')
        print(f'saved as comparison_{timestamp}.pdf')
    plt.show()
    plt.figure(figsize=(16,6))
    plt.subplot(121)
    plt.hist((y-prediction)[z==0]);
    plt.xlabel('y-f(x)')
    plt.title('Error if not Churned')
    plt.gca().axvline(x=0, color='grey', linewidth=1)
    plt.subplot(122)
    plt.hist((y-prediction)[z==1]);
    plt.xlabel('y-f(x)')
    plt.title('Error if Churned')
    plt.gca().axvline(x=0, color='grey', linewidth=1)
    if save:
        plt.tight_layout()
        timestamp = time.time()
        plt.savefig(f'error_{timestamp}.pdf')
        print(f'saved as error_{timestamp}.pdf')
    plt.show()
    
def compute_metrics(y,prediction, z):
    died_p = np.array(prediction[z==1])
    died_y = np.array(y[z==1])
    cindex = ((np.expand_dims(died_y,0)<=np.expand_dims(y,1))&(np.expand_dims(died_p,0)<=np.expand_dims(prediction,1))).sum()/(np.expand_dims(died_y,0)<=np.expand_dims(y,1)).sum()
    print(f'C-Index: {cindex}')
    underestimated = np.mean((prediction-y)[(prediction<y)&(z==1)])
    print(f'Average Underestimated Survival (Churned): {-underestimated}')
    underestimated = np.mean((prediction-y)[(prediction<y)&(z==0)])
    print(f'Average Underestimated Survival (Not Churned): {-underestimated}')


# # Linear Survival SVM

# In[ ]:


# get data, which is comparable: for a given y_i, return X_j(i), y_j(i)
def get_comparable_data(y):
    comparable_y = np.max(observed_y[observed_y<y])
    idx = np.random.randint((observed_y==comparable_y).sum())
    return observed_X[observed_y==comparable_y][idx], observed_y[observed_y==comparable_y][idx]


# In[ ]:
# define y_comp, y_comp_bar: 
# y_comp = y_i for (i,j)in E
# y_comp_bar = y_j for (i,j)in E
# define X_comp, X_comp_bar: 
# X_comp = x_i for (i,j) in E
# X_comp_bar = x_j for (i,j) in E

observed_X = X_train[z_train==1]
observed_y = y_train[z_train==1]

y_comp = y_train[y_train>1]
X_comp = X_train[y_train>1]
data_comp_bar = list(map(get_comparable_data, y_comp))

X_comp_bar = np.vstack([dcb[0] for dcb in data_comp_bar])
y_comp_bar = np.hstack([dcb[1] for dcb in data_comp_bar])


# In[ ]:


### YOUR CODE HERE

### FORMULATE AND SOLVE DM


# In[ ]:


### YOUR CODE HERE
b = ...
prediction = ...
### FIND THE PREDICTION ON THE TRAINING SET


# In[ ]:


print('Training Performance')
# plot the training performance
plot_results(y_train, prediction, z_train, save=False)
compute_metrics(y_train,prediction, z_train)


# In[ ]:


### YOUR CODE HERE
prediction = ...
### FIND THE PREDICTION ON THE TEST SET


# In[ ]:


print('Testing Performance')
# plot the test performance
plot_results(y_test, prediction, z_test, save=False)
compute_metrics(y_test, prediction, z_test)


# # Simple MAE Regression

# In[ ]:


### YOUR CODE HERE

### SOLVE THE MAE REGRESSION PROBLEM


# In[ ]:


### YOUR CODE HERE
prediction_reg = ...
### FIND THE PREDICTION ON THE TRAINING SET


# In[ ]:


print('Training Performance')
plot_results(y_train, prediction_reg, z_train, save=False)
compute_metrics(y_train, prediction_reg, z_train)


# In[ ]:


### YOUR CODE HERE
prediction_reg = ...
### FIND THE PREDICTION ON THE TRAINING SET


# In[ ]:


print('Testing Performance')
plot_results(y_test, prediction_reg, z_test, save=False)
compute_metrics(y_test, prediction_reg, z_test)


# # Kernelized Version

# In[ ]:


# sigma^2 to use for RBF Kernel
sigma2 = 1.5


# In[ ]:


# function to compute RBF Kernel Matrix between X1 and X2
def compute_rbf(X1, X2, sigma2=100):
    return np.exp(-((np.expand_dims(X2, 0)- np.expand_dims(X1, 1))**2).sum(axis=2)/(2*sigma2))


# In[ ]:


### YOUR CODE HERE

### FORMULATE AND SOLVE DM-KERNELIZED


# In[ ]:


### YOUR CODE HERE
b = ...
prediction_rbf = ...
### FIND THE PREDICTION ON THE TRAINING SET


# In[ ]:


print('Training Performance')
# training results
plot_results(y_train, prediction_rbf, z_train, save=False)
compute_metrics(y_train, prediction_rbf, z_train)


# In[ ]:


### YOUR CODE HERE
prediction_rbf = ...
### FIND THE PREDICTION ON THE TEST SET


# In[ ]:


print('Testing Performance')
# testing results
plot_results(y_test, prediction_rbf, z_test, save=False)
compute_metrics(y_test, prediction_rbf, z_test)


# In[ ]:




