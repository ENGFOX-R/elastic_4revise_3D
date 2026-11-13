import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import os, sys
import re
import os.path
import shutil
import numpy as np

filename = "sout_data.csv"

#not use offset and scientific notion in figures.
formatter = ScalarFormatter(useOffset=True)
formatter.set_scientific(False)

# Open the file and read the first line to get the heaNoneders
with open(filename, 'r') as f:
    headers = f.readline().strip().split(',')
#print (f"headers: {headers}")
xlabel='timestep (50 in QE value)'
color=['red','green','blue','black','cyan','magenta','white']

data = np.genfromtxt(filename, names=headers, delimiter=',')
#print (data.dtype.names)
fig, ax = plt.subplots(3)
count = 0
for name in data.dtype.names[1:]:
    data[name] = ["%.6f" % item for item in data[name]]
    ax[count].set_xlabel(xlabel,color='black',fontweight='bold')
    ax[count].set_ylabel(headers[count+1],color='black',fontweight='bold')
    # Apply the formatter to the y-axis
    ax[count].yaxis.set_major_formatter(formatter) 
    ax[count].plot(data['step'], data[name],color=color[count],marker='s',markevery=5,label=name)
    ax[count].grid()# Create a ScalarFormatter object
    count = count + 1
fig.tight_layout()
if os.path.isfile("./QEout_plot.png"):
    os.remove("./QEout_plot.png")
plt.savefig("./QEout_plot.png")

