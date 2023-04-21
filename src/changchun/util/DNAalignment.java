package changchun.util;


public class DNAalignment {
	
	
	public static void main(String[] args) {
		String dna1 = "acbcdb";
		String dna2 = "cadbd";
		
		Alignment score = align(dna1, dna2);
		
		System.out.println(score);
		System.out.println("Score = " + score.score());
	}
	
	
	// Computes and returns the optimal alignment between the
	// two DNA sequences. Sequences do not have to be the same
	// length.
	public static Alignment align(String dna1, String dna2) {
		if (dna1.length() == 0 && dna2.length() == 0) {
			return new Alignment();
		} else if (dna1.length() == 0) {
			Alignment result = align(dna1, dna2.substring(1));
			result.addMatch('-', dna2.charAt(0));
			return result;
		} else if (dna2.length() == 0) {
			Alignment result = align(dna1.substring(1), dna2);
			result.addMatch(dna1.charAt(0), '-');
			return result;
		} else {
			Alignment first = align(dna1.substring(1), dna2);
			first.addMatch(dna1.charAt(0), '-');
			
			Alignment second = align(dna1, dna2.substring(1));
			second.addMatch('-', dna2.charAt(0));
			
			Alignment both = align(dna1.substring(1), dna2.substring(1));
			both.addMatch(dna1.charAt(0), dna2.charAt(0));
			
			if (first.score() >= second.score() && first.score() >= both.score()) {
				return first;
			} else if (second.score() >= first.score() && second.score() >= both.score()) {
				return second;
			} else {
				return both;
			}
		}
	}	
	
	
	//----------------HELPER CLASS------------------
	// Represents an alignment of two strands of DNA.
	public static class Alignment {
		String dna1;
		String dna2;
		
		// Create a new alignment with no characters included.
		public Alignment() {
			dna1 = "";
			dna2 = "";
		}
		
		// Adds c1 to the front of the first DNA sequence,
		// adds c2 to the front of the 2nd DNA sequence.
		public void addMatch(char c1, char c2) {
			dna1 = c1 + dna1;
			dna2 = c2 + dna2;
		}
		
		// Computes the score of this alignment.
		// Match is +2, mismatch is -1.
		public int score() {
			int score = 0;
			for (int i = 0; i < dna1.length(); i++) {
				if (dna1.charAt(i) == dna2.charAt(i)) {
					score += 2;
				} else {
					score -= 1;
				}
			}
			return score;		
		}
		
		// Returns each DNA sequence on its own line.
		public String toString() {
			return dna1 + "\n" + dna2;
		}
	}
}