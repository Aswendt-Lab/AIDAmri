import numpy as np
import scipy.io as io
from scipy.stats import pearsonr

def calculate_p_corr_matrix(data, lines, output_paths):
    correlation_matrix = np.zeros((98,98))
    p_value_matrix = np.zeros((98,98))
    for i in range(98):
        for j in range(i+1, 98):
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
