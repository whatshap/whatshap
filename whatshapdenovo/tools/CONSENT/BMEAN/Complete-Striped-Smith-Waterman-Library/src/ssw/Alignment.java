package ssw;

/**
 * Smith-Waterman alignment
 * 
 * @author Daniel Cameron
 *
 */
public class Alignment {
	/**
	 * score1	the best alignment score
	 */
	public final short score1;
	/**
	 * sub-optimal alignment score
	 */
	public final short score2;
	/**
	 * 0-based best alignment beginning position on reference;	ref_begin1 = -1 when the best alignment beginning position is not available
	 */
	public final int ref_begin1;
	/**
	 * 0-based best alignment ending position on reference
	 */
	public final int ref_end1;
	/**
	 * 0-based best alignment beginning position on read; read_begin1 = -1 when the best alignment beginning
						position is not available
	 */
	public final int read_begin1;
	/**
	 * 0-based best alignment ending position on read
	 */
	public final int read_end1;
	/**
	 * 0-based sub-optimal alignment ending position on read
	 */
	public final int ref_end2;
	/**
	 * best alignment cigar; stored the same as that in BAM format
	 */
	public final String cigar;
	public Alignment(
			short score1,
			short score2,
			int ref_begin1,
			int ref_end1,
			int	read_begin1,
			int read_end1,
			int ref_end2,
			String cigar) {
		this.score1 = score1;
		this.score2 = score2;
		this.ref_begin1 = ref_begin1;
		this.ref_end1 = ref_end1;
		this.read_begin1 = read_begin1;
		this.read_end1 = read_end1;
		this.ref_end2 = ref_end2;
		this.cigar = cigar;
	}
	@Override
	public String toString() {
		return String.format("score1=%d score2=%d, ref=%d,%d read=%d,%d refend2=%d, cigar=%s",
				score1, score2, ref_begin1, ref_end1, read_begin1, read_end1, ref_end2, cigar);
	}
}
