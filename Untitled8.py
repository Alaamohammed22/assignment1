#!/usr/bin/env python
# coding: utf-8

# In[2]:


from pyopenms import *
dig = ProteaseDigestion()
dig.getEnzymeName()
bsa = "".join([l.strip() for l in open("uniprot-yourlist_M202111216320BA52A5CE8FCD097CB85A53697A35301FDED.fasta").readlines()[1:]])
bsa = AASequence.fromString(bsa)
result = []
dig.digest(bsa, result)
print(result[4].toString())
len(result)

