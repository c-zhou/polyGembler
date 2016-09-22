package cz1.test;

import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

public class StreamTest
{
    public static class MyOutputStream extends FilterOutputStream
    {
        public MyOutputStream(OutputStream out)
        {
            super(out);
        }

        @Override
        public void write(byte[] b, int off, int len) throws IOException
        {
            byte[] text = "MyOutputStream called: ".getBytes();         
            super.write(text, 0, text.length);
            super.write(b, off, len);
        }
    }
 
    public static void main(String[] args)
    {       
        PrintStream outStream = new PrintStream(new MyOutputStream(System.out), true); 
        System.setOut(outStream); 

        System.out.println("");
        System.out.print("TEST");
        System.out.println("Another test");
    }
}