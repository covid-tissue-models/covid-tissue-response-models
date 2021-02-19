
# coding: utf-8

# In[9]:


import tellurium as te
import numpy as np
import matplotlib.pyplot as plt


# In[1]:


dosingmodel_str = '''
    model dosingmodel()
    //Time is in days!
    
    //infusion    
    //J0: -> Dpls; switch * infusion_amount / one_our // switch = (0,1) to turn on or off, infusion happens over 1h    
    //flow from plasma     
    J1: Dpls -> ; kE0 * Dpls // elimination    
    J2: Dpls -> Dperi ; kp * Dpls // to periphery     
    J3: Dpls -> Dlung ; k0 * Dpls
        
    // flow from periphery    
    J4: Dperi -> Dpls ; kpp * Dpls // to plasma     
    
    // Drug reactions / flow in lung    
    J5: Dlung -> Dpls ; k0 * Dlung    
    J6: Dlung -> Mala ; k12 * Dlung    
    J7: Dlung -> ; kE1 * Dlung
    
    // Mala reactions    
    J8: Mala -> Mnmp ; k23 * Mala    
    J9: Mala -> ; kE2 * Mala
    
    //Mnmp reactions    
    J10: Mnmp -> Mntp ; k34 * Mnmp
    J11: Mnmp ->  ; kE3 * Mnmp
    
    // Mntp reaction    
    J12: Mntp -> ; kE4 * Mntp
    
    //parameters
    // initial conditions     
    Dperi = 0    
    Dlung = 0
    Mala = 0    
    Mnmp = 0    
    Mntp = 0    
    
    
    // rates    
    k0 = 6.3335
    kE0 = 20.253
    kp = 0.41195
    kpp = 0.36502
    k12 = 1.2248
    kE1 = 7.81
    k23 = 372.61
    kE2 = 6.0801
    k34 = 181.64
    kE3 = 0.97259
    kE4 = 0.83115
    value_drug = 1
    
    const Dpls := value_drug;
    
    // events    
    
    end
'''


# In[3]:


r = te.loada(dosingmodel_str)


# In[4]:


max_Mntp = []
for i in range(0,100):
    r.resetAll()
    r.value_drug = i
    s = r.simulate(0,20,2000)
#     r.plot()
    max_Mntp.append(np.max(r['Mntp']))


# In[20]:


max_Mntp = np.array(max_Mntp)
z = np.polyfit(range(0,100),max_Mntp,1)
f = np.poly1d(z)


# In[21]:


plt.figure()
plt.plot(range(0,100,2), max_Mntp[::2], 'o')
plt.plot(range(0,100), f(range(0,100)))
plt.show()

