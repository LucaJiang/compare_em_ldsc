# Visualize the results of EM and ldsc
# Usage: python3 $code_path/visualize.py -d $output_path -o $output_path
# File: em_results.csv and ldsc_results.csv should be in the output_path
# Output: boxplot
# Input data structure:
# data,h,sigma_beta,causal_rate,hest
# 0228214204,0.1,0.1,0.001,0.077292
# 0229030114,0.1,0.1,0.001,0.068547

import argparse, os
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_theme(style="whitegrid")

# Parse arguments
parser = argparse.ArgumentParser(description='Visualize the results of EM and ldsc')
parser.add_argument('-d', '--data_path', default='.', type=str, help='path to the data')
parser.add_argument('-o', '--output', default='.', type=str, help='path to the output')
args = parser.parse_args()

# Set path
data_path = args.data_path
output_path = args.output
img_name = 'compare_' + datetime.now().strftime("%m%d%H%M%S") + '.png'

def visualize(df, algo):
    '''
    Input: df: pd.DataFrame, algo: str
    Output: list of figures
    '''
    h_list = np.unique(df['h'])
    bys = ['sigma_beta','causal_rate']
    figs = []
    for by in bys:
        fig, ax = plt.subplots()
        sns.boxplot(x='h', y='hest', data=df, hue=by, palette="Set3")
        # draw the h line for each x
        for i in range(len(h_list)):
            plt.plot([i-0.4, i+0.4], [h_list[i], h_list[i]], 'r-')
        plt.plot(0,0,'r-',label='h real')
        plt.title(algo + ' with different ' + by)
        plt.xlabel('h real')
        plt.ylabel('h estimated')
        plt.legend(title=by)
        figs.append(fig)
    return figs

# Read data
em_data_path = os.path.join(data_path, 'em_results.csv')
ldsc_data_path = os.path.join(data_path, 'ldsc_results.csv')
if not os.path.exists(em_data_path) or not os.path.exists(ldsc_data_path):
    print('No data file found!')
    exit(1)
em_data = pd.read_csv(em_data_path)
ldsc_data = pd.read_csv(ldsc_data_path)

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# Draw figures
em_figs = visualize(em_data, 'EM')
ldsc_figs = visualize(ldsc_data, 'LDSC')

fig, axs = plt.subplots(2, 2, figsize=(10, 10))

# Convert figures to images and display them
for i, fig in enumerate(em_figs + ldsc_figs):
    canvas = FigureCanvas(fig)
    canvas.draw()
    img = np.frombuffer(canvas.tostring_rgb(), dtype='uint8').reshape(fig.canvas.get_width_height()[::-1] + (3,))
    axs[i // 2, i % 2].imshow(img)
    axs[i // 2, i % 2].axis('off')

plt.savefig(os.path.join(output_path, img_name))

print('Done at ' + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + '!')
