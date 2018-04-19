

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.TreeMap;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;

public class sumgenescod {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new sumgenescod(args );
	}
	
	sumgenescod(String[] args)
	{
		file = args[0] ;
		file2 = args[1] ;
		file3 = args[1]+"merged.idx" ;
		
		readcsv();
		readidx();
		writecsv();
	} 
	
	String file  =  "";
	String file2  =  "";
	String file3  =  "";
	
	TreeMap<String , String> genes= new TreeMap<String , String> ();
	TreeMap<String , Integer[]> geneexp= new TreeMap<String , Integer[]> ();
	
	public void readcsv()
	{
		try 
		{ 
			CSVReader reader = new CSVReader(new FileReader(file) );
			String[] nextLine;
			 
			//nextLine = reader.readNext();
			 
			while ((nextLine = reader.readNext()) != null) 
			{    
				 genes.put(nextLine[1], nextLine[0]);
				 
			}
			reader.close();

		 

		} catch (IOException ee) {
			System.out.println("error " + ee);
		}
	}
	public void readidx()
	{
		try 
		{ 
			CSVReader reader = new CSVReader(new FileReader(file2),'\t' );
			String[] nextLine;
			 
			//nextLine = reader.readNext();
			 
			while ((nextLine = reader.readNext()) != null) 
			{    
				String g=nextLine[0];
				if(genes.containsKey(nextLine[0]))
				{
					g=genes.get(nextLine[0]);
				}
				
				Integer[] exp = geneexp.get(g);
				if(exp==null)
				{
					exp=new Integer[2];
					geneexp.put(g, exp);
					exp[0] = 0;
					exp[1] = 0;
					 
				}
			 
				exp[1] =exp[1]+ Integer.parseInt(nextLine[2]);
				exp[0] =exp[0]+ 1;
			}
			reader.close(); 

		} catch (IOException ee) {
			System.out.println("error " + ee);
		}
	}
	
	public void writecsv()
	{
		try 
		{ 
			FileWriter fout = new  FileWriter(file3); 
			CSVWriter writer = new CSVWriter(fout);
		 	String[] nextLine;
			 int c=0;
			for(String k:geneexp.keySet())
			{
				Integer[] data=geneexp.get(k);
				String[] data2=new String[3];
				data2[0]=k;
				data2[1]=data[0]+"";
				data2[2]=data[1]+"";
				if(c>0)
				{
					writer.writeNext(data2);
				}
				c++;
			}
			writer.close(); 
			fout.close();

		} catch (IOException ee) {
			System.out.println("error " + ee);
		}
	}
}
