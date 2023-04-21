package changchun;

import java.io.Serializable;
import java.util.LinkedList;

import htsjdk.samtools.util.Interval;

/**
 * 
 * Summary information for one gene in one strain, about how it was edited
 * 
 * @author Johan Henriksson
 *
 */
public class OneEdit implements Serializable {
	private static final long serialVersionUID = 1L;
	
	public String geneid;
	public Interval interval;
	
	public String chrom;
	public int from, to;
	
	
	public int sumdel;
	public int sumins;
	public int numread;
	public int fine;
	public LinkedList<String> reads=new LinkedList<String>();
	public String alignment;

	
	/// After better figuring out region
	public int reducedFrom;		
	public int reducedTo;
	public String alignedGenome;
	public String alignedTemplate;		
	
	
	
}