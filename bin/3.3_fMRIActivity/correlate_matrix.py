import numpy as np
import scipy.io as io
from scipy.stats import pearsonr

def calculate_p_corr_matrix(data, lines, output_paths):
    (rows, cols) = np.shape(data)
    correlation_matrix = np.zeros((cols,cols))
    p_value_matrix = np.zeros((cols,cols))
    for i in range(cols):
        for j in range(i+1, cols):
            corr_coef, p_value = pearsonr(data[:,i], data[:,j])
            correlation_matrix[i, j] = corr_coef
            correlation_matrix[j, i] = corr_coef
            p_value_matrix[i, j] = p_value
            p_value_matrix[j, i] = p_value

    # calculate fisher-transformation
    matrix_PcorrZ = np.arctanh(correlation_matrix)
    
    io.savemat(output_paths[0], dict([('matrix', correlation_matrix),('label',lines)]))
    io.savemat(output_paths[1], dict([('matrix', p_value_matrix),('label',lines)]))
    io.savemat(output_paths[2], dict([('matrix', matrix_PcorrZ),('label',lines)]))
