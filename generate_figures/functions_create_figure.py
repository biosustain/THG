import cobra
import pandas as pd 
from collections import defaultdict
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def pie(counts, ylabel):
    explode = (0, 0, 0)
    colors = ['#191970', '#001CF0', '#0038E2', '#0055D4', '#0071C6', '#008DB8', '#00AAAA',
          '#00C69C', '#00E28E', '#00FF80']

    size=15
    counts.plot(kind='pie', fontsize=size, colors=colors, explode=explode)
    plt.axis('equal')
    plt.ylabel(ylabel,fontsize=size*1.2)
    #plt.legend(labels=counts.index, loc="best",fontsize=12)
    plt.show()


def pie2(counts, ylabel):
    explode = (0, 0)
    colors = ['#191970', '#001CF0', '#0038E2', '#0055D4', '#0071C6', '#008DB8', '#00AAAA',
          '#00C69C', '#00E28E', '#00FF80']

    size=15
    counts.plot(kind='pie', fontsize=size, colors=colors, explode=explode)
    plt.axis('equal')
    plt.ylabel(ylabel,fontsize=size*1.2)
    #plt.legend(labels=counts.index, loc="best",fontsize=12)
    plt.show()


def pie5(counts, ylabel):
    explode = (0, 0, 0, 0, 0)
    colors = ['#191970', '#001CF0', '#0038E2', '#0055D4', '#0071C6', '#008DB8', '#00AAAA',
          '#00C69C', '#00E28E', '#00FF80']

    size=15
    counts.plot(kind='pie', fontsize=size, colors=colors, explode=explode)
    plt.axis('equal')
    plt.ylabel(ylabel,fontsize=size*1.2)
    #plt.legend(labels=counts.index, loc="best",fontsize=12)
    plt.show()