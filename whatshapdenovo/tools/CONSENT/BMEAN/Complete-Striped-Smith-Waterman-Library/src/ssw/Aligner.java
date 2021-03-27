package ssw;


/**
 * JNI wrapper for Striped Smith-Waterman alignment
 * 
 * @author Daniel Cameron
 *
 */
public class Aligner {
	/**
	 * your estimated best alignment score is surely < 255
	 */
	public static final int MAX_SCORE_LESS_THAN_255 = 0;
	/**
	 * if your estimated best alignment score >= 255
	 */
	public static final int MAX_SCORE_GREATER_THAN_255 = 1;
	public static final int MAX_SCORE_UNSURE = 2;
	//public static final int FLAG_UNUSED = 0x80;
	//public static final int FLAG_UNUSED = 0x40;
	//public static final int FLAG_UNUSED = 0x20;
	//public static final int FLAG_UNUSED = 0x10;
	/**
	 *  bit 5: when setted as 1, the function will return the best alignment beginning position
	 */
	public static final int FLAG_INCLUDE_BEST_ALIGNMENT_POSITION = 0x08;
	/**
	 * bit 6: when setted as 1, if (ref_end1 - ref_begin1 < filterd && read_end1 - read_begin1 < filterd), (whatever bit 5 is setted) the function will return the best alignment beginning position and cigar;  
	 */
	public static final int FLAG_INCLUDE_BEST_ALIGNMENT_POSITION_AND_CIGAR_IF_MORE_THAN_FILTER_DISTANCE_POSITIONS_ALIGNED = 0x04;
	/**
	 * bit 7: when setted as 1, if the best alignment score >= filters, (whatever bit 5 is setted) the function will return the best alignment beginning position and cigar;
	 */
	public static final int FLAG_INCLUDE_BEST_ALIGNMENT_POSITION_AND_CIGAR_IF_SCORE_GREATER_THAN_FILTER_SORE = 0x02;
	/**
	 * bit 8: when setted as 1, (whatever bit 5, 6 or 7 is setted) the function will always return the best alignment beginning position and cigar. When flag == 0, only the optimal and sub-optimal scores and the optimal alignment ending position will be returned.
	 */
	public static final int FLAG_INCLUDE_BEST_ALIGNMENT_POSITION_AND_CIGAR = 0x01;
	/**
	 * Do Striped Smith-Waterman alignment. Warning: No parameter checking is performed. Incorrect arguments are likely to crash the JVM.
	 * 
	 * @param read pointer to the query sequence; the query sequence needs to be numbers
	 * @param flattenedMatrix pointer to the substitution matrix; mat needs to be corresponding to the read sequence
	 * @param n the square root of the number of elements in mat (mat has n*n elements)
	 * @param score_size estimated Smith-Waterman score; if your estimated best alignment score is surely < 255 please set 0; if
					your estimated best alignment score >= 255, please set 1; if you don't know, please set 2
	 * @param ref pointer to the target sequence; the target sequence needs to be numbers and corresponding to the mat parameter of function initprofile
	 * @param gapOpen the absolute value of gap open penalty.
	 * @param gapExtend the absolute value of gap extension penalty.
	 * @param flag bitwise FLAG; (from high to low) bit 5: when setted as 1, the function will return the best alignment
					beginning position; bit 6: when setted as 1, if (ref_end1 - ref_begin1 < filterd && read_end1 - read_begin1
					< filterd), (whatever bit 5 is setted) the function will return the best alignment beginning position and
					cigar; bit 7: when setted as 1, if the best alignment score >= filters, (whatever bit 5 is setted) the function
  					will return the best alignment beginning position and cigar; bit 8: when setted as 1, (whatever bit 5, 6 or 7 is
 					setted) the function will always return the best alignment beginning position and cigar. When flag == 0, only
					the optimal and sub-optimal scores and the optimal alignment ending position will be returned.
						Note: Whatever the parameter flag is setted, this function will at least return the optimal and sub-optimal alignment score,
						and the optimal alignment ending positions on target and query sequences. If both bit 6 and 7 of the flag are setted
						while bit 8 is not, the function will return cigar only when both criteria are fulfilled. All returned positions are
						0-based coordinate.
	 * @param filterscore score filter: when bit 7 of flag is setted as 1 and bit 8 is setted as 0, filters will be used (Please check the
 					decription of the flag parameter for detailed usage.)
	 * @param filterdistance distance filter: when bit 6 of flag is setted as 1 and bit 8 is setted as 0, filterd will be used (Please check
					the decription of the flag parameter for detailed usage.)
	 * @param maskLen The distance between the optimal and suboptimal alignment ending position >= maskLen. We suggest to use
					readLen/2, if you don't have special concerns. Note: maskLen has to be >= 15, otherwise this function will NOT
					return the suboptimal alignment information. Detailed description of maskLen: After locating the optimal
					alignment ending position, the suboptimal alignment score can be heuristically found by checking the second
					largest score in the array that contains the maximal score of each column of the SW matrix. In order to avoid
					picking the scores that belong to the alignments sharing the partial best alignment, SSW C library masks the
					reference loci nearby (mask length = maskLen) the best alignment ending position and locates the second largest
					score from the unmasked elements.
	 * @return Smith-Waterman alignment
	 */
	public static native Alignment align(byte[] read, byte[] flattenedMatrix, int n, int score_size, byte[] ref, int gapOpen, int gapExtend, int flag, short filterscore, int filterdistance, int maskLen);
	/**
	 * Performs striped Smith-Waterman alignment
	 * 
	 * @param read read sequence
	 * @param ref reference sequence
	 * @param matrix scoring matrix
	 * @param gapOpen absolute value of gap open penalty
	 * @param gapExtend absolute value of gap extension penalty
	 * @param ignoreCase upper and lowercase sequence values are treated as the same value.
	 * @return Smith-Waterman alignment
	 */
	public static Alignment align(byte[] read, byte[] ref, int[][] matrix, int gapOpen, int gapExtend, boolean ignoreCase) {
		if (gapOpen < 0 || gapExtend < 0) throw new IllegalArgumentException("Gap open and extension penalties must be positive");
		if (gapOpen >= 256 || gapExtend >= 256) throw new IllegalArgumentException("Gap open and extension penalties must fit into unsigned 8-bit integer");
		//if (flag != flag & 0xFF) throw new IllegalArgumentException("Only lowest 8 bits of flag are meaningful");
		int[] lookup = new int[257]; // lookup[256] is used as sentinal for number of unique bases/matrix size
		java.util.Arrays.fill(lookup, -1);
		lookup[256] = 0;
		byte[] readNum = convertToNumeric(lookup, read, ignoreCase);
		byte[] refNum = convertToNumeric(lookup, ref, ignoreCase);
		byte[] flattenedMatrix = flatten(lookup, matrix);
		int uniqueBases = lookup[256];
		assert(flattenedMatrix.length == uniqueBases * uniqueBases);
		assert(maxValue(readNum) < uniqueBases);
		assert(maxValue(refNum) < uniqueBases);
		Alignment alignment = align(
				readNum, flattenedMatrix, uniqueBases, MAX_SCORE_UNSURE,
				refNum, gapOpen, gapExtend, FLAG_INCLUDE_BEST_ALIGNMENT_POSITION_AND_CIGAR, (short)0, 0, Math.max(15, readNum.length / 2));
		return alignment;
	}
	private static int maxValue(byte[] array) {
		int max = Integer.MIN_VALUE;
		for (int i = 0; i < array.length; i++) {
			if (array[i] > max) {
				max = array[i];
			}
		}
		return max;
	}
	/**
	 * Converts an ASCII sequence into numeric successive 0-based values 
	 * @param lookup ASCII to numeric conversion lookup array
	 * @param sequence ASCII sequence
	 * @param ignoreCase treat upper and lowercase ASCII values as the same 
	 * @return numeric sequence
	 */
	private static byte[] convertToNumeric(int[] lookup, byte[] sequence, boolean ignoreCase) {
		byte[] numericSeq = new byte[sequence.length];
		for (int i = 0; i < sequence.length; i++) {
			int b = sequence[i];
			if (ignoreCase) {
				b = Character.toUpperCase(b);
			}
			if (lookup[b] == -1) {
				lookup[b] = lookup[256]++;
			}
			numericSeq[i] = (byte)lookup[b];
		}
		return numericSeq;
	}
	/**
	 * Generates a flattened ssw scoring matrix  
	 * @param lookup ASCII to numeric conversion lookup array
	 * @param matrix scoring matrix
	 * @return flattened ssw numeric scoring matrix
	 */
	private static byte[] flatten(int[] lookup, int[][] matrix) {
		int size = lookup[256];
		byte[] flattened = new byte[size * size];
		for (int i = 0; i < matrix.length; i++) {
			int newi = lookup[i];
			if (newi == -1) continue;
			for (int j = 0; j < matrix[i].length; j++) {
				int newj = lookup[j];
				if (newj == -1) continue;
				int score = matrix[i][j];
				if (score < Byte.MIN_VALUE || score > Byte.MAX_VALUE) {
					throw new IllegalArgumentException("Scoring matrix values must fit into signed 8-bit integer");
				}
				flattened[newi * size + newj] = (byte)score;
			}
		}
		return flattened;
	}
}
