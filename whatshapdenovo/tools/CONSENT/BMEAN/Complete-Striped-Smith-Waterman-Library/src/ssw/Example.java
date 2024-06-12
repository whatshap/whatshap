package ssw;

/**
 * Example code for SSW JNI wrapper
 * 
 * @author Daniel Cameron
 *
 */
public class Example {
	public static void main(String[] args) {
		try {
			System.loadLibrary("sswjni");
		} catch (java.lang.UnsatisfiedLinkError e) {
			System.out.println(String.format("Cannot find libsswjni.so. Has the library been built and LD_LIBRARY_PATH or -Djava.library.path set appropriately?\n%s", e));
			throw e;
		}
		int[][] score = new int[128][128];
		for (int i = 0; i < 128; i++) {
			for (int j = 0; j < 128; j++) {
				if (i == j) score[i][j] = 2;
				else score[i][j] = -2;
			}
		}
		System.out.println("Aligning nucleotides");
		Alignment aln = Aligner.align("CTGAGCCGGTAAATC".getBytes(), "CAGCCTTTCTGACCCGGAAATCAAAATAGGCACAACAAA".getBytes(), score, 3, 1, true);
		if (aln == null) {
			throw new RuntimeException();
		}
		System.out.print(String.format("score1=%d ", aln.score1));
		System.out.print(String.format("score2=%d ", aln.score2));
		System.out.print(String.format("ref_begin1=%d ", aln.ref_begin1));
		System.out.print(String.format("ref_end1=%d ", aln.ref_end1));
		System.out.print(String.format("read_begin1=%d ", aln.read_begin1));
		System.out.print(String.format("read_end1=%d ", aln.read_end1));
		System.out.print(String.format("ref_end2=%d ", aln.ref_end2));
		System.out.print(String.format("cigar=%s ", aln.cigar));
		System.out.println();
	}
}
