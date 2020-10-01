package cz1.util;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.BufferedWriter;
import java.io.BufferedReader;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPOutputStream;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.nio.file.Files;
import java.util.Calendar;
import java.util.UUID;
import java.text.SimpleDateFormat;

public class Utils {
	
	private static final Logger myLogger = LogManager.getLogger(Utils.class);

	public static BufferedWriter getGZIPBufferedWriter(String path) throws IOException {
		return new BufferedWriter(new OutputStreamWriter(new
				GZIPOutputStream(new FileOutputStream(path))));
	}

	public static String getSystemTime(){
		return new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").
				format(Calendar.getInstance().getTime());
	}
	
	public static boolean makeOutputDir(File file) {
		// TODO Auto-generated method stub
		if(!file.exists() || file.exists()&&!file.isDirectory()) {
			file.mkdir();
			return true;
		} else {
			return false;
		}
	}
	
	public static void deleteDirectory(File file) {
		File[] contents = file.listFiles();
	    if (contents != null) {
	        for (File f : contents) {
	            if (! Files.isSymbolicLink(f.toPath())) {
	            	deleteDirectory(f);
	            }
	        }
	    }
	    file.delete();
	}

    public static BufferedReader getBufferedReader(String inSourceName) {

        try {
            if (inSourceName.startsWith("http")) {
                if (inSourceName.endsWith(".gz")) {
                    return new BufferedReader(new InputStreamReader(new GZIPInputStream((new URL(inSourceName)).openStream())));
                } else {
                    return new BufferedReader(new InputStreamReader((new URL(inSourceName)).openStream()));
                }
            } else {
                if (inSourceName.endsWith(".gz")) {
                    return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inSourceName))));
                } else {
                    return new BufferedReader(new InputStreamReader(new FileInputStream(inSourceName)));
                }
            }
        } catch (Exception e) {
            myLogger.error("getBufferedReader: Error getting reader for: " + inSourceName);
            e.printStackTrace();
        }
        return null;
    }

    public static BufferedReader getBufferedReader(String inSourceName, int bufSize) {

        try {
            if (bufSize < 1) {
                return getBufferedReader(inSourceName);
            } else {
                if (inSourceName.startsWith("http")) {
                    if (inSourceName.endsWith(".gz")) {
                        return new BufferedReader(new InputStreamReader(new GZIPInputStream((new URL(inSourceName)).openStream())), bufSize);
                    } else {
                        return new BufferedReader(new InputStreamReader((new URL(inSourceName)).openStream()), bufSize);
                    }
                } else {
                    if (inSourceName.endsWith(".gz")) {
                        return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inSourceName))), bufSize);
                    } else {
                        return new BufferedReader(new InputStreamReader(new FileInputStream(inSourceName)), bufSize);
                    }
                }
            }
        } catch (Exception e) {
            myLogger.error("getBufferedReader: Error getting reader for: " + inSourceName);
            e.printStackTrace();
        }
        return null;
    }

    public static BufferedReader getBufferedReader(File file, int bufSize) {
        return getBufferedReader(file.getAbsolutePath(), bufSize);
    }

    public static BufferedReader getBufferedReader(InputStream is) {
		// TODO Auto-generated method stub
		return new BufferedReader(new InputStreamReader(is));
	}
    
    public static BufferedReader getBufferedReader(InputStream is, int bufSize) {
		// TODO Auto-generated method stub
		return new BufferedReader(new InputStreamReader(is), bufSize);
	}
    
    public static BufferedWriter getBufferedWriter(String filename) {
        return getBufferedWriter(filename, false);
    }

    public static BufferedWriter getBufferedWriter(String filename, boolean append) {
        return getBufferedWriter(new File(filename), append);
    }

    public static BufferedWriter getBufferedWriter(File file) {
        return getBufferedWriter(file, false);
    }

    public static BufferedWriter getBufferedWriter(File file, boolean append) {

        String filename = null;

        try {
            filename = file.getCanonicalPath();
            if (filename.endsWith(".gz")) {
                return new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(file, append))));
            } else {
                return new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file, append)));
            }
        } catch (Exception e) {
            myLogger.error("getBufferedReader: Error getting reader for: " + filename);
            e.printStackTrace();
        }
        return null;
    }

    public static DataOutputStream getDataOutputStream(String filename, int bufSize) {

        try {
            if (filename.endsWith(".gz")) {
                return new DataOutputStream(new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(filename), bufSize)));
            } else {
                return new DataOutputStream(new BufferedOutputStream(new FileOutputStream(filename), bufSize));
            }
        } catch (Exception e) {
            myLogger.error("getDataOutputStream: Error getting reader for: " + filename);
            e.printStackTrace();
        }
        return null;

    }
    
    public static InputStream getInputStream(String filename) {

        try {
            if (filename.startsWith("http")) {
                if (filename.endsWith(".gz")) {
                    return new GZIPInputStream((new URL(filename)).openStream());
                } else {
                    return (new URL(filename)).openStream();
                }
            } else {
                if (filename.endsWith(".gz")) {
                    return new GZIPInputStream(new FileInputStream(filename));
                } else {
                    return new FileInputStream(filename);
                }
            }
        } catch (Exception e) {
            myLogger.error("getInputStream: Error getting reader for: " + filename);
            e.printStackTrace();
        }
        return null;
    }
    
	public static String paste(double[] array, String collapse) {
		StringBuilder s = new StringBuilder();
		s.append(array[0]);
		for(int i=1; i<array.length; i++) {
			s.append(collapse);
			s.append(array[i]);
		}
		return s.toString();
	}
	
	public static String paste(int[] array, String collapse) {
		StringBuilder s = new StringBuilder();
		s.append(array[0]);
		for(int i=1; i<array.length; i++) {
			s.append(collapse);
			s.append(array[i]);
		}
		return s.toString();
	}
	
	public static String fixedLengthPaddingString(String str, int len) {
	    return String.format("%1$"+len+ "s", str);
	}

	public static String paste(String[] array, String collapse) {
		// TODO Auto-generated method stub
		StringBuilder s = new StringBuilder();
		s.append(array[0]);
		for(int i=1; i<array.length; i++) {
			s.append(collapse);
			s.append(array[i]);
		}
		return s.toString();
	}

	public static String paste(boolean[] array, String collapse) {
		// TODO Auto-generated method stub
		StringBuilder s = new StringBuilder();
		s.append(array[0]);
		for(int i=1; i<array.length; i++) {
			s.append(collapse);
			s.append(array[i]);
		}
		return s.toString();
	}

	public static String makeTempDir() {
		// TODO Auto-generated method stub
		String uuid = "."+UUID.randomUUID().toString();
		while(new File(uuid).exists())
			uuid = UUID.randomUUID().toString();
		return uuid;
	}

	public static void renameFile(String file1, String file2) {
		// TODO Auto-generated method stub
		new File(file1).renameTo(new File(file2));
	}

	public static void deleteFile(String file) {
		// TODO Auto-generated method stub
		new File(file).delete();
	}
	
	
}
