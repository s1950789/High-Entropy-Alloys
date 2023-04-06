#!/usr/bin/env python
# coding: utf-8

# In[2]:


# Required repositories
import numpy as np
from scipy import constants
import itertools
import matplotlib.pyplot as plt


# In[3]:


# Provide Elemental Data for Radii, Atomic volume, Shear Modulus and Poisson Ratio
Elemental_data = {
    "Goldschmidt Radii": {'Ni':1.24, 'Co':1.25, 'Fe':1.26, 'Mn':1.27, 'Cr':1.28},
    "Atomic Volume": {'Ni':10.94, 'Co':11.12, 'Fe':12.09, 'Mn':12.60, 'Cr':12.27},
    "Shear Modulus": {'Ni':76, 'Co':76, 'Fe':82, 'Mn':44, 'Cr':115},
    "Poisson Ratio": { 'Ni':0.31, 'Co':0.31, 'Fe':0.28, 'Mn':0.26, 'Cr':0.21}
}


# Enthalpy Matrix for CANTOR consituents
enthalpy_matrix = ([['Ni', 'Co', 'Fe', 'Mn', 'Cr'],
                    [0,    -21,  -97, -115,  -31],
                    [-21,    0,  -60,  -19,    5],
                    [-97,  -60,    0,    9,   -8],
                    [-115, -19,    9,    0, -110],
                    [-30,    5,   -8, -110,    0]])


# In[4]:


# Elememt class to store elemental data (symbol, ratio, radii, volume, shear modulus, poisson ratio)
class Element:
    def __init__(self, symbol, composition, radii, volume, shear_modulus, poisson_ratio):
        self.symbol = symbol
        self.composition = composition 
        self.radii = radii
        self.volume = volume
        self.shear_modulus = shear_modulus 
        self.poisson_ratio = poisson_ratio
        
    def __str__(self):
        return f"Symbol: {self.symbol}, Composition: {self.composition}, Radii: {self.radii}, Volume: {self.volume}, Shear Modulus: {self.shear_modulus}, Poisson Ratio: {self.poisson_ratio}"


# In[5]:


# Effective medium class to store effective medium data (averaged: radii, shear modulus, poissons ratio, volume 
# misfit: delta, volume & burgers vector)

class effective_medium:
    def __init__(self, average_radii, average_volume, shear_modulus, poissons_ratio, 
                 misfit_delta, misfit_volume, burgers_vector):
        self.radii = r_ave
        self.average_volume = eff_v
        self.shear_modulus = mu
        self.poissons_ratio = v
        self.misfit_delta = delta 
        self.misfit_volume = misfit_vol
        self.burgers_vector = b
        
        


# In[6]:


# Experimental conditions/constants (temperature, strain_rate, boltzamn, alpha, ref_strain, taylors_factor, f_1, f_2)
class cond:
    def __init__(temperature, strain_rate, boltzamn, alpha, ref_strain, taylors_factor, f_1, f_2):
        self.temperature = T
        self.strain_rate = e
        self.boltzmann = k
        self.alpha = a
        self.ref_strain = e_0
        self.taylor_factor = taylor_factor
        self.f_1 = f_1
        self.f_2 = f_2
        
# Set constant conditions
cond.e = 0.001
cond.k = constants.k
cond.a = 0.123
cond.e_0 = 10**4
cond.taylor_factor = 3.06
cond.f_1 = 0.35
cond.f_2 = 5.70


# In[7]:


# Check an input composition is complete
def check_comp(input_comp):
    comp = sum(input_comp.values())
    if comp == 1:
        pass
    else:
        print("Incorrect Composition: total composition does not equal 1")
        print(comp)


# In[8]:


# Check data available for an input element
def check_elements(input_comp, database):
    for element in input_comp:
        if element not in database["Goldschmidt Radii"]:
            print("Incorrect Composition: elemental data not available for, " + element)
    pass


# In[9]:


# From input compostion initiate element class 
def initiate_elements (input_comp, database): 
    
    check_comp(input_comp)
    check_elements(input_comp, database)
    
    elements = {}
    
    for element in input_comp:
        comp = input_comp[element]
        radii = Elemental_data['Goldschmidt Radii'][element] * 10**(-10)
        volume = Elemental_data['Atomic Volume'][element] * 10**(-30)
        shear = Elemental_data['Shear Modulus'][element] * 10**9
        poissons = Elemental_data['Poisson Ratio'][element]

        elements[element] = Element(element, comp, radii, volume, shear, poissons)

    return elements


# In[10]:


# Function to caluclate weighted parameters based on composition ratios
def weighted_average(input_comp, param, elements):
    weighted_values = []
    
    for element in input_comp:
        comp = elements[element].composition
        parameter = getattr(elements[element], param)
        
        weighted_val = comp*parameter
        weighted_values.append(weighted_val)
        
    weighted_ave = (sum(weighted_values))
    
    return weighted_ave


# In[11]:


# From input compostion calculate effective properties & initiate class 
def effective_properties(input_comp, elements):
    # Shear Modulus GPa
    effective_medium.mu = weighted_average(input_comp, "shear_modulus", elements)

    # Poissons Ratio
    effective_medium.v = weighted_average(input_comp, "poisson_ratio", elements)

    # Effective volume
    effective_medium.eff_vol = weighted_average(input_comp, "volume", elements)

    # Average metallic radii
    effective_medium.r_ave = weighted_average(input_comp, "radii", elements)
    
    # Burgers vector
    effective_medium.b = (4*effective_medium.eff_vol)**(1/3)/(np.sqrt(2))


# In[12]:


# Calulate Delta misfit value for a compositon
def delta_misfit(input_comp, elements):
    values = []
    for element in input_comp:
        r_i = elements[element].radii
        conc = elements[element].composition
        val = conc*((1-(r_i/effective_medium.r_ave))**2)
        values.append(val)

    effective_medium.delta = np.sqrt(sum(values))


# In[13]:


# Calulate volume misfit value for a compositon
def volume_misfit(input_comp, elements):
    values = []
    for element in input_comp:
        v_i = elements[element].volume
        conc = elements[element].composition
        val = conc*((v_i - effective_medium.eff_vol)**2)
        values.append(val)
        
    effective_medium.vol_misfit = sum(values)


# In[14]:


# Initialise yield stress calulation: calculate tau_y0 & deltaE_b. 
# Specify chosen misfit parameter: volume misfit or delta misfit
def initialise_yield(misfit_param):
    a = cond.a
    mu = effective_medium.mu
    v = effective_medium.v
    f_1 = cond.f_1
    eff_vol = effective_medium.eff_vol
    delta = effective_medium.delta
    vol_misfit = effective_medium.vol_misfit
    b = effective_medium.b
    f_2 = cond.f_2
    
    if str(misfit_param) == 'delta':
        # Zero-K yield stress
        tau_y0 = 0.051*a**(-1/3)*mu*((1+v)/(1-v))**(4/3)*f_1*((9*(eff_vol**2)*(delta**2))/(b**6))**(2/3)

        # Energy barrier of edge dislocation
        deltaE_b = 0.274*a**(1/3)*mu*(b**3)*(((1+v)/(1-v))**(2/3))*f_2*((9*(eff_vol**2)*(delta**2))/(b**6))**(1/3)
        
    elif str(misfit_param) == 'volume':
        # Zero-K yield stress
        tau_y0 = 0.051*a**(-1/3)*mu*((1+v)/(1-v))**(4/3)*f_1*((vol_misfit)/(b**6))**(2/3)

        # Energy barrier of edge dislocation
        deltaE_b = 0.274*a**(1/3)*mu*(b**3)*(((1+v)/(1-v))**(2/3))*f_2*((vol_misfit)/(b**6))**(1/3)
    
    else: 
        print("Select either 'delta' or 'volume'")
        
    return tau_y0, deltaE_b


# In[15]:


# Calulate tau_y
def calc_tau_y(tau_y0, deltaE_b):
    k = cond.k
    T = cond.T
    e_0 = cond.e_0
    e = cond.e
    
    # Low Temperature Model:
    tau_y = tau_y0*(1 - ( ((k*T)/deltaE_b)*np.log(e_0/e))**(2/3) )
    
    if tau_y/tau_y0 > 0.5:
        tau_y = tau_y
        
    else: 
        # High Temperature Model:     
        tau_y = tau_y0*np.exp(-(1/0.51)*( ((k*T)/deltaE_b)*np.log(e_0/e))**(2/3) )
        
    return tau_y


# In[16]:


# Calulate yield stress
def calc_yield(tau_y):
    taylor_factor = cond.taylor_factor
    yield_stress = tau_y * taylor_factor

    return yield_stress


# In[17]:


# Combining calculation functions. Returns yield stress from: composition & temperature
# specify misfit param
def calculate_yield_stress(input_comp, temperature, misfit_param):
    cond.T = temperature

    # Get elemental data for specied elements
    elements = initiate_elements (input_comp, Elemental_data)

    # Calculate composition weighted properties
    effective_properties(input_comp, elements)

    # Calculate misfit parameters
    delta_misfit(input_comp, elements)
    volume_misfit(input_comp, elements)

    # Calculate tau_y0 and deltaE_b using delta or vol_misfit param
    tau_y0, deltaE_b = initialise_yield(misfit_param)
    
    tau_y = calc_tau_y(tau_y0, deltaE_b)
    
    yield_stress = calc_yield(tau_y)*10**(-6)
    
    return yield_stress


# In[18]:


# Function to return yield stresses of multiple compositions
def vary_composition(input_comps, temperature, misfit_param):
    yield_stresses = []
    Labels = []
        
    for composition in input_comps:
        yield_stress = calculate_yield_stress(composition, temperature, misfit_param)
        yield_stresses.append(yield_stress)
        
        name = result_string = ''.join(list(composition.keys()))
        Labels.append(name)
        
    return yield_stresses, Labels


# In[19]:


# Function to return yield stress of varying temperature
def vary_temperature(input_comps, min_temp, max_temp, step, misfit_param):

    yield_stresses = []
    
    temp_range = range(min_temp, max_temp, step)
    
    for t in temp_range:
        yield_stress = calculate_yield_stress(input_comps, t, misfit_param)
        yield_stresses.append(yield_stress)
        
    return temp_range, yield_stresses


# In[20]:


# Generate all possible unique equiatomic elements from CANTOR elements
def all_unique_combinations(database):
    
    available_elements = list(Elemental_data["Goldschmidt Radii"].keys())
    
    # Generate all combinations with a minimum of 3 elements
    all_compositions = []
    
    
    # Compositions of 3 or more elements  
    for element in range(3, len(available_elements)+1):
        combinations = itertools.combinations(available_elements, element)
        
        # Caluculate ratio of elements         
        for combo in combinations:
            # Set ratio to 0.33, 0.33, 0.34 for combinations with 3 elements
            if len(combo) == 3:
                comp = {combo[0]: 0.33, combo[1]: 0.33, combo[2]: 0.34}
            else:
                ratio = 1.0 / len(combo)  # equal ratio for all elements
                comp = {e: ratio for e in combo}
                
            all_compositions.append(comp)
    
    return all_compositions


# In[21]:


# Calulate delta values for list of compositions
def get_deltas(input_comps, database):
    
    deltas = []
    
    for composition in input_comps:
        
        # Get elemental data for specied elements
        elements = initiate_elements (composition, database)

        # Calculate composition weighted properties
        effective_properties(composition, elements)

        # Calculate misfit parameters
        delta_misfit(composition, elements)

        delta = effective_medium.delta
        deltas.append(delta)
        
        
    return deltas


# In[22]:


# Create a labeled plot for different compositions
def labelled_plot(x, y, x_label, y_label, element_labels):
    
    # create lists for the legend
    legend_colors = []
    legend_labels = []

    # iterate over the data points
    for i in range(len(x)):
        # create a scatter plot for each point with a unique color
        color = plt.cm.jet(i/len(x))
        plt.scatter(x[i], y[i], color=color)

        # append the color and label to the legend lists
        legend_colors.append(color)
        legend_labels.append(element_labels[i])
        

    # add legend
    plt.legend(legend_labels, facecolor='white', framealpha=1, ncol =1,  bbox_to_anchor=(1.35, 0.5), loc='center right')
    

    # set the labels for the axes
    plt.xlabel(str(x_label))
    plt.ylabel(str(y_label))

    # show the plot
    plt.show()



# In[23]:


# Plot & label max y point
def plot_ymax(x, y):
    
    y_max = max(y)
    
    index = y.index(y_max)
    x_coord = x[index]
    
    plt.scatter(x_coord, y_max, color = "r")
    
    label = "(" + str(round(x_coord, 1)) + ", " + str(round(y_max,1)) + ")"
    plt.annotate(label, (x_coord, y_max*0.9))


# In[24]:


# Calculate ideal entropy of mixing for a composition
def get_entropy(input_comp):
    
    R = constants.R
    
    s_values = []
    
    for element_frac in input_comp.values():
        s_i = element_frac*np.log(element_frac)
        s_values.append(s_i)
        
    s = -R*sum(s_values)
    return s


# In[25]:


# Calculate enthalpy of mixing for a composition
def get_enthlapy(enthalpy_matrix, input_comp):

    elements = list(input_comp.keys())

    pairs = []
    for i in range(len(elements)):
        for j in range(i+1, len(elements)):
            pairs.append((elements[i], elements[j]))

    enthalpys = []
    for pair in pairs:
        i_index = enthalpy_matrix[0].index(pair[0])
        j_index = enthalpy_matrix[0].index(pair[1])
        enthalpy_ij = enthalpy_matrix[j_index+1][i_index]
        c_i = input_comp[pair[0]]
        c_j = input_comp[pair[1]]
        enthalpys.append(c_i * c_j * enthalpy_ij)

    total_enthalpy = 4 * sum(enthalpys)
    
    return float(total_enthalpy)
    


# In[ ]:




