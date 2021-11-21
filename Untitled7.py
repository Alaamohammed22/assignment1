#!/usr/bin/env python
# coding: utf-8

# In[6]:


pip install pyopenms


# In[8]:


import pyopenms as ms
sum_=0
seq = ms.AASequence.fromString("VAKA")
for aa in seq:
    sum_+=aa.getMonoWeight()
sum_


# In[10]:


seq = ms.AASequence.fromString("VAKA")
total=seq.getMonoWeight()
total


# In[11]:


total==sum_

