import numpy as np

def Diversity1(A, threshold):
    O = np.sum(A > threshold)
    return O

def Diversity2(A, threshold):
    A_ = A + 1e-6
    O = -np.sum(A_ * np.log(A_))
    return O

def Diversity3(A, threshold):
    if np.sum(A) < 0.9:
        O = 0
    else:
        O = 1 / np.sum(A**2)
    return O

def SimilarityBC(A, B, threshold):
    A = A * (A > threshold)
    B = B * (B > threshold)
    O = np.sum(np.minimum(A, B))
    return O

def SimilarityJS(A, B, threshold):
    A = A * (A > threshold)
    B = B * (B > threshold)
    A_ = A + 1e-10
    B_ = B + 1e-10
    x = (A_ + B_) / 2
    O = 1 - np.sqrt(0.5 * np.sum(A_ * (np.log(A_) - np.log(x))) + 0.5 * np.sum(B_ * (np.log(B_) - np.log(x))))
    return O
    #O=1-sqrt(0.5*sum(A_.*(log(A_)-log(x)))+0.5*sum(B_.*(log(B_)-log(x))));

def SimilarityJ(A, B, threshold):
    A = A * (A > threshold)
    B = B * (B > threshold)
    O = np.sum((A > threshold) & (B > threshold)) / (np.sum((A > threshold) | (B > threshold) )+ 1e-10)
    return O

def SimilarityTo1(A, A1, A2, threshold, f_similarity):
    O = f_similarity(A, A1, threshold) / ((f_similarity(A, A1, threshold) + f_similarity(A, A2, threshold))+ 1e-10)
    return O

def SimilarityTo2(A, A1, A2, threshold, f_similarity):
    O = f_similarity(A, A2, threshold) / ((f_similarity(A, A1, threshold) + f_similarity(A, A2, threshold))+ 1e-10)
    return O

def Additivity1(A, A1, A2, threshold):
    O = np.sum(A > threshold) / (np.sum(A1 > threshold) + np.sum(A2 > threshold) - np.sum((A1 > threshold) & (A2 > threshold)))
    return O

def Additivity2(A, A1, A2, threshold):
    O = np.sum((A > threshold) & ((A1 > threshold) | (A2 > threshold))) / np.sum((A1 > threshold) | (A2 > threshold))
    return O

def Additivity3(A, A1, A2, threshold):
    O = np.sum((A > threshold) & ((A1 > threshold) | (A2 > threshold))) / (np.sum(A1 > threshold) + np.sum(A2 > threshold))
    return O