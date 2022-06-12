import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib.pyplot import MultipleLocator

if(len(sys.argv) != 2):
	print("please input file name")
else:
	file_name = (sys.argv[1])

df_news = pd.read_table(file_name,header = None)

x = df_news[0]
y = df_news[1]

plt.plot(x,y,'r',label="path")
plt.axis('equal')
plt.legend()
plt.grid('on')
plt.show()
