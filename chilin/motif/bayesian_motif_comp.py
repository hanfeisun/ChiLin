# Author: Cliff Meyer, Dana-Farber Cancer Institute, 2008
# Clustering of position weight matrices based on
# "A Novel Bayesian DNA Motif Comparison Method for Clustering and Retrieval"

# Habib et al, 2008
# PLOS Computational Biology, Volume 4, Issue 2

import numpy


SPLINTERCUTOFF = 0.6
SENSE = 0
ANTISENSE = 1

def log_E_p_standard(n):
    """Log of the standard dirichlet prior."""
    E_p = n + 1
    E_p /= numpy.tile(E_p.sum(1), (4,1)).T
    return numpy.log(E_p)



def BLiC_score(M1,M2):
    """
    if flag_r is ANTISENSE the reverse complement is the better match
    returns alignment score shift and orientation of M1 in the match
    matrix A is shifted before reversing
    """
    n1 = M1.shape[0]
    n2 = M2.shape[0]
    n = min(n1,n2)

    if n1 > n2:
        A,B = M1,M2 # A is longer one
        i_max = 1
    else:
        A,B = M2,M1
        n1,n2 = n2,n1
        i_max = -1

    max_score = -999
    #for i in range( abs(n1-n2) + 1 ):
    #    score   = BLiC_score_aligned( A[i:i+n, :], B )
    #    #score_r = BLiC_score_aligned( A[i+n:i:-1, ::-1], B )
    #    score_r = BLiC_score_aligned( A[i:i+n, :][::-1,::-1], B )
    #    print score, score_r
    for i in range(1-n2, n1):
        if i<0:
            Bsub = B[-i:n2-i, :]
            Brev_sub = B[::-1,::-1][-i:n2-i, :]
            Asub = A[:n2+i, :]
        elif i <= n1-n2:
            Bsub = B
            Brev_sub = B[::-1,::-1]
            Asub = A[i:i+n2, :]
        elif n1-i < n2:
            Asub = A[i:n1, :]      #B:    xCATCGCxxx
            Bsub = B[:n1-i, :]     #A: xxxxxxTCGC
            Brev_sub = B[::-1,::-1][:n1-i, :]

        score   = BLiC_score_aligned( Asub, Bsub )
        score_r = BLiC_score_aligned( Asub, Brev_sub )
        #print i, '\t', Asub.shape, '\t', score, score_r
        if score_r > score:
            flag, score = ANTISENSE, score_r
        else:
            flag = SENSE

        if score > max_score:
            max_score = score
            flag_r = flag
            ii = i

    i_max *= ii
    return max_score, i_max, flag_r

def BLiC_score_aligned_old(M1, M2):

    A1 = numpy.zeros(M1.shape)
    for i in xrange(M1.shape[0]):
        sum_M1_i = M1[i].sum()
        if sum_M1_i != 0:
            A1[i] = M1[i] / sum_M1_i
        else:
            A1[i] = M1[i]

    A2 = numpy.zeros(M1.shape)
    for i in range(M1.shape[0]):
        sum_M2_i = M2[i].sum()
        if sum_M2_i != 0:
            A2[i] = M2[i] / sum_M2_i
        else:
            A2[i] = M2[i]

    A12 = A1 + A2

    log_p1 = log_E_p_standard(A1)
    log_p2 = log_E_p_standard(A2)
    log_p12 = log_p1
    log_pBG = log_E_p_standard(numpy.ones(A1.shape))

    s = 2 * (A12 * log_p12).sum() - (A1 * log_p1 + A2 * log_p2 + A12 * log_pBG).sum()
    return s


def BLiC_score_aligned(M1, M2):
    #for i in xrange(M1.shape[0]):

    sum_M1_i = M1.sum(axis=1)
    sum_M1_i = 1.0*(sum_M1_i==0) + sum_M1_i
    A1 = M1.transpose()/sum_M1_i
    A1 = A1.transpose()

    sum_M2_i = M2[0:M1.shape[0],].sum(axis=1)
    sum_M2_i = 1.0*(sum_M2_i==0) + sum_M2_i
    A2 = M2[0:M1.shape[0],].transpose()/sum_M2_i
    A2 = A2.transpose()

    #A2 = numpy.zeros(M1.shape)
    #for i in range(M1.shape[0]):
    #    sum_M2_i = M2[i].sum()
    #    if sum_M2_i != 0:
    #        A2[i] = M2[i] / sum_M2_i
    #    else:
    #        A2[i] = M2[i]

    A12 = A1 + A2

    log_p1 = log_E_p_standard(A1)
    log_p2 = log_E_p_standard(A2)
    log_p12 = log_p1
    log_pBG = log_E_p_standard(numpy.ones(A1.shape))

    s = 2 * (A12 * log_p12).sum() - (A1 * log_p1 + A2 * log_p2 + A12 * log_pBG).sum()
    return s


def merge( M1, M2, shift, antisense, normalize=True ):
    """
    Merge PWMs M1 and M2 where shift is used to align matrics of unequal length.
    Rows are scaled to the max total count of either matrix in that row.
    """

    #print 'M1', M1
    #print 'M2', M2
    #print 'shift', shift
    #print 'anti', antisense

    ashift = abs( shift )
    M3 = numpy.zeros( ( max( M1.shape[0], M2.shape[0] ), 4 ), float )
    n = min( M1.shape[0], M2.shape[0] )

    if M1.shape[0] > M2.shape[0]:
        A,B = M1,M2
    else:
        A,B = M2,M1

    #if antisense == True:
    #    print 'A', A[::-1,::-1]
    #    print 'B', B
    #    M3[:,:] = A[::-1,::-1]
    #    M3[M3.shape[0]-abs(shift)-n : M3.shape[0]-abs(shift)] += B
    #else:
    #    M3[:,:] = A[:,:]
    #    M3[abs(shift) : n+abs(shift)] += B

    try:
        if antisense == True:
            #if ashift > 0:
            #    M3 = A[ ashift+n:ashift:-1, ::-1 ] + B
            #else:
            #    M3 = A[ n::-1, ::-1 ] + B
            M3 = A[ashift:ashift+n,:][::-1,::-1] + B
        else:
            M3 = A[ashift:ashift+n,:] + B
    except:
        print 'exception', A.shape, B.shape, A[ashift:ashift+n,:][::-1,::-1].shape

        #scale each row
    if normalize == True:
        for i in xrange(M3.shape[0]):
            if M3[i].sum() != 0:
                M3[i] /= M3[i].sum()

    return M3
