# QGADMM
Quantized GADMM

# Linear regression
To run the main code that calls all algorithms for linear regression task, run main_QGADMM.m
To generate the energy CDF plot for the linear regression task, you need to run Q_GADMM_energy_CDFPlot.m
  Not that you should have run main_QGADMM.m first and saved all results in a mat file that will be called inside Q_GADMM_energy_CDFPlot.m.
  
# Image classification using DNN 
To run GADMM, Q-GADMM, SGD, and QSGD algorithms for image classification problem using DNN, you need to run these codes:
  GADMM_deepNN.ipynb, QGADMM_deepNN.ipynb, SGD.ipynb, and QSGD.ipynb
The main codes generate the results and save them, to plot the results for the main figure, the energy CDF, and the sesitivity analysis, you need to run the following codes after you generate the required input files.
  plot_main_result_dnn.m, plotCDF_energy_dnn.m, plot_results_dnn.m
  
# Requirements
Matlab, python 3.7.4, and tenserflow 2.0, Jupyter notebook.

# Dataset

# Citation
@misc{elgabli2019qgadmm,
    title={Q-GADMM: Quantized Group ADMM for Communication Efficient Decentralized Machine Learning},
    author={Anis Elgabli and Jihong Park and Amrit S. Bedi and Chaouki Ben Issaid and Mehdi Bennis and Vaneet Aggarwal},
    year={2019},
    eprint={1910.10453},
    archivePrefix={arXiv},
    primaryClass={cs.LG}
}

  
 
