#!/usr/bin/env python
"""
Simple python wrapper for SSW library
Please put the path of libssw.so into LD_LIBRARY_PATH or pass it explicitly as a parameter
By Yongan Zhao (March 2016)
"""

import sys
import os.path as op
import ctypes as ct



lBlosum50 = [
	#  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
     	5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -1, -1, -5,	# A
       -2,  7, -1, -2, -4,  1,  0, -3,  0, -4, -3,  3, -2, -3, -3, -1, -1, -3, -1, -3, -1,  0, -1, -5,	# R
       -1, -1,  7,  2, -2,  0,  0,  0,  1, -3, -4,  0, -2, -4, -2,  1,  0, -4, -2, -3,  5,  0, -1, -5,	# N
       -2, -2,  2,  8, -4,  0,  2, -1, -1, -4, -4, -1, -4, -5, -1,  0, -1, -5, -3, -4,  6,  1, -1, -5,	# D
       -1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -3, -1, -5,	# C
       -1,  1,  0,  0, -3,  7,  2, -2,  1, -3, -2,  2,  0, -4, -1,  0, -1, -1, -1, -3,  0,  4, -1, -5,	# Q
       -1,  0,  0,  2, -3,  2,  6, -3,  0, -4, -3,  1, -2, -3, -1, -1, -1, -3, -2, -3,  1,  5, -1, -5,	# E
     	0, -3,  0, -1, -3, -2, -3,  8, -2, -4, -4, -2, -3, -4, -2,  0, -2, -3, -3, -4, -1, -2, -1, -5,	# G
       -2,  0,  1, -1, -3,  1,  0, -2, 10, -4, -3,  0, -1, -1, -2, -1, -2, -3,  2, -4,  0,  0, -1, -5,	# H
       -1, -4, -3, -4, -2, -3, -4, -4, -4,  5,  2, -3,  2,  0, -3, -3, -1, -3, -1,  4, -4, -3, -1, -5,	# I
       -2, -3, -4, -4, -2, -2, -3, -4, -3,  2,  5, -3,  3,  1, -4, -3, -1, -2, -1,  1, -4, -3, -1, -5,	# L
       -1,  3,  0, -1, -3,  2,  1, -2,  0, -3, -3,  6, -2, -4, -1,  0, -1, -3, -2, -3,  0,  1, -1, -5,	# K
       -1, -2, -2, -4, -2,  0, -2, -3, -1,  2,  3, -2,  7,  0, -3, -2, -1, -1,  0,  1, -3, -1, -1, -5,	# M
       -3, -3, -4, -5, -2, -4, -3, -4, -1,  0,  1, -4,  0,  8, -4, -3, -2,  1,  4, -1, -4, -4, -1, -5,	# F
       -1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -1, -1, -5,	# P
     	1, -1,  1,  0, -1,  0, -1,  0, -1, -3, -3,  0, -2, -3, -1,  5,  2, -4, -2, -2,  0,  0, -1, -5,	# S
    	0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  2,  5, -3, -2,  0,  0, -1, -1, -5, 	# T
       -3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1,  1, -4, -4, -3, 15,  2, -3, -5, -2, -1, -5, 	# W
       -2, -1, -2, -3, -3, -1, -2, -3,  2, -1, -1, -2,  0,  4, -3, -2, -2,  2,  8, -1, -3, -2, -1, -5, 	# Y
     	0, -3, -3, -4, -1, -3, -3, -4, -4,  4,  1, -3,  1, -1, -3, -2,  0, -3, -1,  5, -3, -3, -1, -5, 	# V
       -2, -1,  5,  6, -3,  0,  1, -1,  0, -4, -4,  0, -3, -4, -2,  0,  0, -5, -3, -3,  6,  1, -1, -5, 	# B
       -1,  0,  0,  1, -3,  4,  5, -2,  0, -3, -3,  1, -1, -4, -1,  0, -1, -2, -2, -3,  1,  5, -1, -5, 	# Z
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5, 	# X
       -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1 	# *
       ]



class CAlignRes(ct.Structure):
    """
    @typedef	structure of the alignment result
    @field	nScore	the best alignment score
    @field	nScore2	sub-optimal alignment score
    @field	nRefBeg	0-based best alignment beginning position on reference;	ref_begin1 = -1 when the best alignment beginning
                                            position is not available
    @field	nRefEnd	0-based best alignment ending position on reference
    @field	nQryBeg	0-based best alignment beginning position on read; read_begin1 = -1 when the best alignment beginning
                                            position is not available
    @field	nQryEnd	0-based best alignment ending position on read
    @field	nRefEnd2	0-based sub-optimal alignment ending position on read
    @field	sCigar	best alignment cigar; stored the same as that in BAM format, high 28 bits: length, low 4 bits: M/I/D (0/1/2);
                                    cigar = 0 when the best alignment path is not available
    @field	nCigarLen	length of the cigar string; cigarLen = 0 when the best alignment path is not available
    """
    _fields_ = [('nScore', ct.c_uint16), 
                ('nScore2', ct.c_uint16), 
                ('nRefBeg', ct.c_int32), 
                ('nRefEnd', ct.c_int32), 
                ('nQryBeg', ct.c_int32), 
                ('nQryEnd', ct.c_int32), 
                ('nRefEnd2', ct.c_int32), 
                ('sCigar', ct.POINTER(ct.c_uint32)), 
                ('nCigarLen', ct.c_int32)] 



class CProfile(ct.Structure):
    """
    @typedef	structure of the query profile
    @field	pByte	byte array for profile
    @field	pWord	word array for profile
    @field	pRead	number array for read
    @field	pMat	score matrix
    @field	nReadLen	read length
    @field	nN	edge length of score matrix
    @field	nBias	bias
    """
    _fields_ = [('pByte', ct.POINTER(ct.c_int32)),
                ('pWord', ct.POINTER(ct.c_int32)),
                ('pRead', ct.POINTER(ct.c_int8)),
                ('pMat', ct.POINTER(ct.c_int8)),
                ('nReadLen', ct.c_int32),
                ('nN', ct.c_int32),
                ('nBias', ct.c_uint8)]



class CSsw(object):
    """
    A class for libssw
    """
    def __init__(self, sLibPath):
        """
        init all para
        @para   sLibpath    argparse object
        """
# load libssw
        sLibName = 'libssw.so'
        if not sLibPath:
# user doesn't give the path explicitly
            if not op.exists(op.join(sLibPath, sLibName)):
                print >> sys.stderr, 'libssw.so does not exist in the input path'
                sys.exit()
            self.ssw = ct.cdll.LoadLibrary(op.join(sLibPath,sLibName))
        else:
# otherwise just search in PATH
            bFound = False
            for s in sys.path:
                if op.exists(op.join(s,sLibName)):
                    bFound = True
                    self.ssw = ct.cdll.LoadLibrary(op.join(s,sLibName))
            if bFound == False:
                print >> sys.stderr, 'libssw.so does not exist in PATH'
                sys.exit()

# init ssw_init
        """
	@function	Create the query profile using the query sequence.
	@param	read	pointer to the query sequence; the query sequence needs to be numbers
	@param	readLen	length of the query sequence
	@param	mat	pointer to the substitution matrix; mat needs to be corresponding to the read sequence
	@param	n	the square root of the number of elements in mat (mat has n*n elements)
	@param	score_size	estimated Smith-Waterman score; if your estimated best alignment score is surely < 255 please set 0; if
						your estimated best alignment score >= 255, please set 1; if you don't know, please set 2
	@return	pointer to the query profile structure
	@note	example for parameter read and mat:
			If the query sequence is: ACGTATC, the sequence that read points to can be: 1234142
			Then if the penalty for match is 2 and for mismatch is -2, the substitution matrix of parameter mat will be:
			//A  C  G  T
			  2 -2 -2 -2 //A
			 -2  2 -2 -2 //C
			 -2 -2  2 -2 //G
			 -2 -2 -2  2 //T
			mat is the pointer to the array {2, -2, -2, -2, -2, 2, -2, -2, -2, -2, 2, -2, -2, -2, -2, 2}
        """
        self.ssw_init = self.ssw.ssw_init
        self.ssw_init.argtypes = [ct.POINTER(ct.c_int8), ct.c_int32, ct.POINTER(ct.c_int8), ct.c_int32, ct.c_int8]
        self.ssw_init.restype = ct.POINTER(CProfile)
# init init_destroy
        """
	@function	Release the memory allocated by function ssw_init.
	@param	p	pointer to the query profile structure
        """
        self.init_destroy = self.ssw.init_destroy
        self.init_destroy.argtypes = [ct.POINTER(CProfile)]
        self.init_destroy.restype = None
# init ssw_align
        """
!	@function	Do Striped Smith-Waterman alignment.
	@param	prof	pointer to the query profile structure
	@param	ref	pointer to the target sequence; the target sequence needs to be numbers and corresponding to the mat parameter of
				function ssw_init
	@param	refLen	length of the target sequence
	@param	weight_gapO	the absolute value of gap open penalty
	@param	weight_gapE	the absolute value of gap extension penalty
	@param	flag	bitwise FLAG; (from high to low) bit 5: when setted as 1, function ssw_align will return the best alignment
					beginning position; bit 6: when setted as 1, if (ref_end1 - ref_begin1 < filterd && read_end1 - read_begin1
					< filterd), (whatever bit 5 is setted) the function will return the best alignment beginning position and
					cigar; bit 7: when setted as 1, if the best alignment score >= filters, (whatever bit 5 is setted) the function
  					will return the best alignment beginning position and cigar; bit 8: when setted as 1, (whatever bit 5, 6 or 7 is
 					setted) the function will always return the best alignment beginning position and cigar. When flag == 0, only
					the optimal and sub-optimal scores and the optimal alignment ending position will be returned.
	@param	filters	score filter: when bit 7 of flag is setted as 1 and bit 8 is setted as 0, filters will be used (Please check the
 					decription of the flag parameter for detailed usage.)
	@param	filterd	distance filter: when bit 6 of flag is setted as 1 and bit 8 is setted as 0, filterd will be used (Please check
					the decription of the flag parameter for detailed usage.)
	@param	maskLen	The distance between the optimal and suboptimal alignment ending position >= maskLen. We suggest to use
					readLen/2, if you don't have special concerns. Note: maskLen has to be >= 15, otherwise this function will NOT
					return the suboptimal alignment information. Detailed description of maskLen: After locating the optimal
					alignment ending position, the suboptimal alignment score can be heuristically found by checking the second
					largest score in the array that contains the maximal score of each column of the SW matrix. In order to avoid
					picking the scores that belong to the alignments sharing the partial best alignment, SSW C library masks the
					reference loci nearby (mask length = maskLen) the best alignment ending position and locates the second largest
					score from the unmasked elements.
	@return	pointer to the alignment result structure
	@note	Whatever the parameter flag is setted, this function will at least return the optimal and sub-optimal alignment score,
			and the optimal alignment ending positions on target and query sequences. If both bit 6 and 7 of the flag are setted
			while bit 8 is not, the function will return cigar only when both criteria are fulfilled. All returned positions are
			0-based coordinate.
        """
        self.ssw_align = self.ssw.ssw_align
        self.ssw_align.argtypes = [ct.c_void_p, ct.POINTER(ct.c_int8), ct.c_int32, ct.c_uint8, ct.c_uint8, ct.c_uint8, ct.c_uint16, ct.c_int32, ct.c_int32]
        self.ssw_align.restype = ct.POINTER(CAlignRes)
# init align_destroy
        """
	@function	Release the memory allocated by function ssw_align.
	@param	a	pointer to the alignment result structure
        """
        self.align_destroy = self.ssw.align_destroy
        self.align_destroy.argtypes = [ct.POINTER(CAlignRes)]
        self.align_destroy.restype = None



def read_matrix(sFile):
    """
    read a score matrix for either DNA or protein
    assume the format of the input score matrix is the same as that of http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt

    """
    with open(args.sMatrix, 'r') as f:
        for l in f:
            if not l.startswith('#'):
                break
        lEle = l.strip().split()
        dEle2Int = {}
        dInt2Ele = {}
        for i,ele in enumerate(lEle):
            dEle2Int[ele] = i
            dEle2Int[ele.lower()] = i
            dInt2Ele[i] = ele
        nEleNum = len(lEle)
        lScore = []
        for l in f:
            lScore.extend([int(x) for x in l.strip().split()[1:]])

        return lEle, dEle2Int, dInt2Ele, lScore


if __name__ == '__main__':
    pass
