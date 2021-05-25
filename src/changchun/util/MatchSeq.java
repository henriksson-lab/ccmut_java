package changchun.util;


/**
 * 
 * Match original sequence with edited sequence. Somewhat quick and dirty. Assumes perfect match at the stard and end.
 * Does not attempt to reverse complement either sequence
 * 
 * @author Johan Henriksson
 *
 */
public class MatchSeq {

	public String fasta;
	public String newseq;
	
	public int fasta_i1, fasta_i2;  //indices of the part that is to be removed
	public String fastaSub;
	
	
	public int newseq_i1;
	public int newseq_i2;
	
	public boolean match(String fasta, String newseq) {
		return match(fasta, newseq, 30);
	}
	
	public boolean match(String fasta, String newseq, int len) {
		
		this.fasta=fasta;
		this.newseq=newseq;
		
		String seq1=newseq.substring(0,len);
		String seq2=newseq.substring(newseq.length()-len);
		
		//Replace the content. Return new seq
		fasta_i1=fasta.indexOf(seq1);
		fasta_i2=fasta.indexOf(seq2);
		
		if(fasta_i1!=-1 && fasta_i2!=-1) {
		
			//Calculate substitute
			int fasta_end=fasta_i2+len;
			fastaSub=fasta.substring(0,fasta_i1)+newseq+fasta.substring(fasta_end);

			
			//Narrow in on the precise interval: left side
			newseq_i1=len;
			fasta_i1+=len;
			while(fasta.charAt(fasta_i1)==newseq.charAt(newseq_i1)) {
				fasta_i1++;
				newseq_i1++;
			}
			
			//Narrow in on the precise interval: right side
			newseq_i2=newseq.length()-len;
			while(fasta.charAt(fasta_i2)==newseq.charAt(newseq_i2)) {
				fasta_i2--;
				newseq_i2--;
			}
			
			return true;
		} else
			return false;
	}
	
	
	public static void main(String[] args) {
		
		MatchSeq m=new MatchSeq();
		
		
		
		System.out.println(m.match("123123123123123ABCabcabcabcabcDEF456456456456", "ABCeeeeeeepDEF",3));
		System.out.println(m.fastaSub);
		
		System.out.println("-- "+m.fasta_i1+"  "+m.fasta_i2);
		System.out.println(m.fasta.substring(m.fasta_i1,m.fasta_i2+1));
		
		System.out.println("-- "+m.newseq_i1+"  "+m.newseq_i2);
		System.out.println(m.newseq.substring(m.newseq_i1,m.newseq_i2+1));
		
	}
}
