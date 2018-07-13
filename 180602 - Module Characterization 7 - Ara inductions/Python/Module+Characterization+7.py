
# coding: utf-8

# In[1]:


get_ipython().magic('matplotlib inline')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import least_squares


# In[2]:


data = pd.read_excel('data.xlsx')
# This data is GFP (GeoMean MEFL)
data
######################## Some reminders ####################################
# MyData.x                                  # Data of column called "x"
# MyData.iloc[0,:]                          # First row
# MyData.iloc[:,0]                          # First column
# MyData.columns[0]                         # Name of first column
# MyData['Seconds'] = 3600 * MyData.Hours   # Create a new column named as such and caluclated as such


# In[3]:


# Let's take a look at the rawish data


# In[4]:


plt.plot(data.arabinose_uM, data.pL2f1EV, label='_nolegend_')
plt.plot(data.arabinose_uM, data.pL2f1510, label='_nolegend_')
plt.plot(data.arabinose_uM, data.pL2f1513, label='_nolegend_')
plt.plot(data.arabinose_uM, data.pL2f1542, label='_nolegend_')
plt.plot(data.arabinose_uM, data.pL2f1532, label='_nolegend_')
plt.plot(data.arabinose_uM, data.pL2f1537, label='_nolegend_')

plt.scatter(data.arabinose_uM, data.pL2f1EV)
plt.scatter(data.arabinose_uM, data.pL2f1510)
plt.scatter(data.arabinose_uM, data.pL2f1513)
plt.scatter(data.arabinose_uM, data.pL2f1542)
plt.scatter(data.arabinose_uM, data.pL2f1532)
plt.scatter(data.arabinose_uM, data.pL2f1537)

plt.legend(bbox_to_anchor=(1,0.75), fontsize=12)

plt.xlabel('Arabinose (uM)', fontsize=12)
plt.xscale('log')
plt.xlim(1,max(data.arabinose_uM))

plt.ylabel('GFP (GeoMean MEFL)', fontsize=12)
plt.yscale('log')
plt.ylim(5,3*10**5)

plt.grid()
plt.savefig('GFP vs Arabinose all.png', transparent=True, bbox_inches='tight', dpi = 400, size=(8,6))
plt.show()


# In[5]:


# The thing is that pL2f1542, pL2f1532, pL2f1537 aren't Arabinose --> GFP
# They are Arabinose --> Transcription Factor --> GFP
# So we need to use the pL2f1513 transfer curve to convert Arabinose --> Transcription Factor
# Let's build a model for pL2f1513


# In[6]:


def act_hill_function(x, basal, maximal, Kd, n):
    return basal + maximal * (x**n / (Kd**n + x**n))

def log_act_hill_function(x, basal, maximal, Kd, n):
    return np.log10(basal + maximal * (x**n / (Kd**n + x**n)))

def rep_hill_function(x, basal, maximal, Kd, n):
    return basal + maximal * (1 / (x / Kd)**n)

def log_rep_hill_function(x, basal, maximal, Kd, n):
    return np.log10(basal + maximal * (1 / (x / Kd)**n))

list_of_params_to_fit = ['Basal', 'Max', 'Kd', 'n'] 
def report_paramaters(fit_param_names, fit_param_values):
    for each in range(len(fit_param_names)):
        print(fit_param_names[each], 'is', np.round(fit_param_values[each], 1))


# In[7]:


def pL2f1513_error_function(current_parameter_guess):
    current_parameter_guess = tuple(current_parameter_guess)
    y_guess = log_act_hill_function(data.arabinose_uM, *current_parameter_guess) 
    error = y_guess - np.log10(data.pL2f1513)
    return error

initial_guess = (100, 100, 100, 2)
low_bounds = [0, 0, 0, 0]
up_bounds = [10000000, 10000000, 10000000, 10]
pL2f1513_parameters = least_squares(pL2f1513_error_function,
                                    initial_guess, 
                                    bounds=(low_bounds, 
                                            up_bounds)
                                    ).x
report_paramaters(list_of_params_to_fit, pL2f1513_parameters)


# In[8]:


plt.scatter(data.arabinose_uM, data.pL2f1513, c='b', label='Data')
plt.plot(data.arabinose_uM, act_hill_function(data.arabinose_uM, *pL2f1513_parameters), c='r', label='least squares fit')

plt.legend(loc = 'best')

plt.xlabel('Arabinose (uM)', fontsize=12)
plt.xscale('log')
plt.xlim(1,max(data.arabinose_uM))

plt.ylabel('GFP (GeoMean MEFL)', fontsize=12)
plt.yscale('log')
plt.ylim(5,3*10**5)

plt.grid()
plt.savefig('GFP vs Arabinose pL2f1513.png', transparent=True, bbox_inches='tight', dpi = 400, size=(8,6))
plt.show()


# In[9]:


# Now we have a model of how much transcription factor pBAD-RBS30, the thing driving GFP in PL2f1513, makes
# Notice that this isn't a particularly useful **unit**
# The units are something like "MEFL of GFP equivalent transcripts"


# In[10]:


# So now let's build a model for the other plasmids pL2f1542, pL2f1532, and pL2f1537
# So we need to convert Arabinose into transcription factor (TF) using the model from pL2f1513
# (The reason we don't just use the pL2f1513 data set is to be generic)
# (In the future we might not use the same arabinose concentrations for both pL2f1513 and the plasmid of interest)

data['TF_MEFL'] = act_hill_function(data.arabinose_uM, *pL2f1513_parameters)
xsmooth = np.linspace(min(data.TF_MEFL), max(data.TF_MEFL), 1000000)
data


# In[11]:


# Let's replot the total data with this new x-axis

plt.plot(data.TF_MEFL, data.pL2f1542, label='_nolegend_')
plt.plot(data.TF_MEFL, data.pL2f1532, label='_nolegend_')
plt.plot(data.TF_MEFL, data.pL2f1537, label='_nolegend_')

plt.scatter(data.TF_MEFL, data.pL2f1542)
plt.scatter(data.TF_MEFL, data.pL2f1532)
plt.scatter(data.TF_MEFL, data.pL2f1537)

plt.suptitle('GFP vs TF')
plt.legend(loc='best', fontsize=12)

plt.xlabel('TF MEFL', fontsize=12)
plt.xscale('log')

plt.ylabel('GFP (GeoMean MEFL)', fontsize=12)
plt.yscale('log')
plt.ylim(5,3*10**5)

plt.grid()
plt.savefig('GFP vs TF_MEFL all.png', transparent=True, bbox_inches='tight', dpi = 400, size=(8,6))
plt.show()


# In[12]:


# Let's build a model for pL2f1532

def pL2f1532_error_function(current_parameter_guess):
    current_parameter_guess = tuple(current_parameter_guess)
    y_guess = log_rep_hill_function(data.TF_MEFL, *current_parameter_guess) 
    error = y_guess - np.log10(data.pL2f1532)
    return error

initial_guess = (100, 100, 100, 2)
low_bounds = [0, 0, 0, 0]
up_bounds = [10000000, 10000000, 10000000, 10]
pL2f1532_parameters = least_squares(pL2f1532_error_function,
                                    initial_guess, 
                                    bounds=(low_bounds, 
                                            up_bounds)
                                    ).x
report_paramaters(list_of_params_to_fit, pL2f1532_parameters)

plt.scatter(data.TF_MEFL, data.pL2f1532, c='b', label='Data')
plt.plot(xsmooth, rep_hill_function(xsmooth, *pL2f1532_parameters), c='r', label='least squares fit')

plt.legend(loc = 'best')
plt.suptitle('pL2f1532')

plt.xlabel('TF_MEFL', fontsize=12)
plt.xscale('log')

plt.ylabel('GFP (GeoMean MEFL)', fontsize=12)
plt.yscale('log')
plt.ylim(5,3*10**5)

plt.grid()
plt.savefig('GFP vs TF_MEFL pL2f1532.png', transparent=True, bbox_inches='tight', dpi = 400, size=(8,6))
plt.show()


# In[13]:


# Let's build a model for pL2f1542

def pL2f1532_error_function(current_parameter_guess):
    current_parameter_guess = tuple(current_parameter_guess)
    y_guess = log_act_hill_function(data.TF_MEFL, *current_parameter_guess) 
    error = y_guess - np.log10(data.pL2f1542)
    return error

initial_guess = (100, 100, 100, 2)
low_bounds = [0, 0, 0, 0]
up_bounds = [10000000, 10000000, 10000000, 10]
pL2f1542_parameters = least_squares(pL2f1532_error_function,
                                    initial_guess, 
                                    bounds=(low_bounds, 
                                            up_bounds)
                                    ).x
report_paramaters(list_of_params_to_fit, pL2f1542_parameters)

plt.scatter(data.TF_MEFL, data.pL2f1542, c='b', label='Data')
plt.plot(xsmooth, act_hill_function(xsmooth, *pL2f1542_parameters), c='r', label='least squares fit')

plt.legend(loc = 'best')
plt.suptitle('pL2f1542')

plt.xlabel('TF_MEFL', fontsize=12)
plt.xscale('log')

plt.ylabel('GFP (GeoMean MEFL)', fontsize=12)
plt.yscale('log')
plt.ylim(5,3*10**5)

plt.grid()
plt.savefig('GFP vs TF_MEFL pL2f1542.png', transparent=True, bbox_inches='tight', dpi = 400, size=(8,6))
plt.show()


# In[14]:


# Let's build a model for pL2f1537

def pL2f1537_error_function(current_parameter_guess):
    current_parameter_guess = tuple(current_parameter_guess)
    y_guess = log_rep_hill_function(data.TF_MEFL, *current_parameter_guess) 
    error = y_guess - np.log10(data.pL2f1537)
    return error

initial_guess = (100, 100, 100, 2)
low_bounds = [0, 0, 0, 0]
up_bounds = [10000000, 10000000, 10000000, 10]
pL2f1537_parameters = least_squares(pL2f1537_error_function,
                                    initial_guess, 
                                    bounds=(low_bounds, 
                                            up_bounds)
                                    ).x
report_paramaters(list_of_params_to_fit, pL2f1537_parameters)


plt.scatter(data.TF_MEFL, data.pL2f1537, c='b', label='Data')
plt.plot(xsmooth, rep_hill_function(xsmooth, *pL2f1537_parameters), c='r', label='least squares fit')

plt.legend(loc = 'best')
plt.suptitle('pL2f1537')

plt.xlabel('TF_MEFL', fontsize=12)
plt.xscale('log')

plt.ylabel('GFP (GeoMean MEFL)', fontsize=12)
plt.yscale('log')
plt.ylim(5,3*10**5)

plt.grid()
plt.savefig('GFP vs TF_MEFL pL2f1537.png', transparent=True, bbox_inches='tight', dpi = 400, size=(8,6))
plt.show()

