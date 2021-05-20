package changchun;

import java.io.Serializable;
import java.util.LinkedList;

import htsjdk.samtools.util.Interval;

public class OneEdit implements Serializable {
	private static final long serialVersionUID = 1L;
	public String geneid;
	public Interval interval;
	public int sumdel;
	public int sumins;
	public int numread;
	public int fine;
	public LinkedList<String> reads=new LinkedList<String>();
	public String alignment;		
}