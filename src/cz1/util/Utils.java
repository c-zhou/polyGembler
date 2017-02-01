/*
 * Utils.java
 *
 * Created on May 27, 2003, 2:03 AM
 */
package cz1.util;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.net.URL;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.apache.log4j.Logger;

/**
 *
 * @author terryc
 */
public final class Utils {

    private static final Logger myLogger = Logger.getLogger(Utils.class);
    private static Collection<String> myJavaPackages = null;

    private Utils() {
        // Utility Class
    }

    /**
     * Returns the base name of a string delimited with periods (i.e. Java
     * Class).
     *
     * @param str string to parse
     *
     * @return base name
     */
    public static String getBasename(String str) {
        int index = str.lastIndexOf('.');
        index++;
        return str.substring(index);
    }

    /**
     * This returns the filename only. Preceding directories are removed and
     * everything after last . is removed.
     *
     * @param str original filename
     * @param suffix suffix
     *
     * @return trimmed filename
     */
    public static String getFilename(String str) {

        int indexForwardSlash = str.lastIndexOf('/');
        int indexBackwardSlash = str.lastIndexOf('\\');

        int index = 0;
        if ((indexForwardSlash == -1) && (indexBackwardSlash == -1)) {
            index = 0;
        } else if (indexForwardSlash > indexBackwardSlash) {
            index = indexForwardSlash + 1;
        } else {
            index = indexBackwardSlash + 1;
        }

        String result = str.substring(index);
        if (result.indexOf('.') > 0) {
            result = result.substring(0, result.indexOf('.'));
        }

        return result;

    }

    /**
     * This returns the filename only. Preceding directories are removed and
     * suffix. If suffix not found, then everything after last . is removed.
     *
     * @param str original filename
     * @param suffix suffix
     *
     * @return trimmed filename
     */
    public static String getFilename(String str, String suffix) {

        int indexForwardSlash = str.lastIndexOf('/');
        int indexBackwardSlash = str.lastIndexOf('\\');

        int index = 0;
        if ((indexForwardSlash == -1) && (indexBackwardSlash == -1)) {
            index = 0;
        } else if (indexForwardSlash > indexBackwardSlash) {
            index = indexForwardSlash + 1;
        } else {
            index = indexBackwardSlash + 1;
        }

        String result = str.substring(index);
        if ((suffix != null) && (result.lastIndexOf(suffix) > 0)) {
            result = result.substring(0, result.lastIndexOf(suffix));
        } else if (result.lastIndexOf('.') > 0) {
            result = result.substring(0, result.lastIndexOf('.'));
        }

        return result;

    }

    public static List<String> getFullyQualifiedClassNames(String simpleName) {

        if (myJavaPackages == null) {
            myJavaPackages = getJavaPackages();
        }

        List<String> fqns = new ArrayList<String>();
        for (String aPackage : myJavaPackages) {
            try {
                String fqn = aPackage + "." + simpleName;
                Class.forName(fqn);
                fqns.add(fqn);
            } catch (Exception e) {
                // Do Nothing
            }
        }
        return fqns;

    }

    public static Collection<String> getJavaPackages() {
        String classpath = System.getProperty("java.class.path");
        return getPackagesFromClassPath(classpath);
    }

    public static Set<String> getPackagesFromClassPath(String classpath) {
        Set<String> packages = new HashSet<String>();
        String[] paths = classpath.split(File.pathSeparator);
        for (String path : paths) {
            if (path.trim().length() == 0) {
                continue;
            } else {
                File file = new File(path);
                if (file.exists()) {
                    String childPath = file.getAbsolutePath();
                    if (childPath.endsWith(".jar")) {
                        packages.addAll(readZipFile(childPath));
                    } else {
                        packages.addAll(readDirectory(childPath));
                    }
                }
            }

        }
        return packages;
    }

    public static Set<String> readDirectory(String path) {
        Set<String> packages = new HashSet<String>();
        File file = new File(path);
        int startIndex = path.length() + 1;
        for (File child : file.listFiles()) {
            recursiveRead(child, startIndex, packages);
        }
        return packages;
    }

    public static void recursiveRead(File file, int startIndex, Set<String> packages) {
        if (!file.isFile()) {
            packages.add(file.getAbsolutePath().substring(startIndex).replace(File.separator, "."));
            for (File child : file.listFiles()) {
                recursiveRead(child, startIndex, packages);
            }
        }
    }

    public static Set<String> readZipFile(String path) {
        Set<String> packages = new HashSet<String>();
        try {
            ZipFile zFile = new ZipFile(path);
            Enumeration<? extends ZipEntry> entries = zFile.entries();
            while (entries.hasMoreElements()) {
                ZipEntry entry = entries.nextElement();
                if (!entry.isDirectory()) {
                    String dirName = new File(entry.getName()).getParent();
                    if (dirName != null) {
                        String name = dirName.replace(File.separator, ".");
                        if (name.endsWith(".")) {
                            name = name.substring(0, name.length() - 1);
                        }
                        packages.add(name);
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return packages;
    }

    public static String shortenStrLineLen(String str, int preferredLen) {
        return shortenStrLineLen(str, preferredLen, -1);
    }

    /**
     */
    public static String shortenStrLineLen(String str, int preferredLen, int preferredLines) {

        StringBuilder buffer = new StringBuilder();

        int startIndex = 0;
        int endIndex = preferredLen;
        int strLen = str.length();
        int numLines = 0;

        while (startIndex < strLen - 1) {

            int place = str.indexOf(' ', endIndex);
            int newLine = str.indexOf('\n', startIndex);

            if ((newLine != -1) && (newLine < place)) {
                place = newLine;
            }

            String part = null;
            if (place == -1) {
                part = str.substring(startIndex);
                buffer.append(part);
                buffer.append("\n");
                break;
            } else {
                place++;
                part = str.substring(startIndex, place);
                buffer.append(part);
                buffer.append("\n");
                startIndex = place;
                endIndex = place + preferredLen;
            }

            numLines++;

            if ((preferredLines > 0) && (numLines >= preferredLines)) {
                return buffer.toString();
            }

        }

        return buffer.toString();

    }

    /**
     * Adds suffix (i.e. .txt) to end of filename if it's not already there.
     *
     * @param filename filename
     * @param suffix suffix
     *
     * @return filename with suffix
     */
    public static String addSuffixIfNeeded(String filename, String suffix) {

        String temp = filename.toLowerCase();

        if (suffix.charAt(0) != '.') {
            suffix = '.' + suffix;
        }

        if (temp.endsWith(suffix)) {
            return filename;
        } else {
            return filename + suffix;
        }

    }

    /**
     * Adds default suffix if not already one of the possible suffixes.
     *
     * @param filename filename
     * @param defaultSuffix default suffix
     * @param possible possible suffixes
     *
     * @return filename with suffix
     */
    public static String addSuffixIfNeeded(String filename, String defaultSuffix, String[] possible) {

        String temp = filename.toLowerCase();

        for (int i = 0; i < possible.length; i++) {
            String current = possible[i].toLowerCase();
            if (current.charAt(0) != '.') {
                current = '.' + current;
            }

            if (temp.endsWith(current)) {
                return filename;
            }

        }

        if (defaultSuffix.charAt(0) != '.') {
            defaultSuffix = '.' + defaultSuffix;
        }
        return filename + defaultSuffix;

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

    /**
     * Finds index of Nth occurrence of character in string.
     *
     * @param str string
     * @param match character to match
     * @param n Nth occurrence
     *
     * @return index
     */
    public static int findNthOccurrenceInString(String str, char match, int n) {
        int result = str.indexOf(match);
        while (--n > 0 && result != -1) {
            result = str.indexOf(match, result + 1);
        }
        return result;
    }

    /**
     * Returns max heap size in MB.
     *
     * @return max heap size
     */
    public static long getMaxHeapSizeMB() {
        return Runtime.getRuntime().maxMemory() / 1048576l;
    }

    /**
     * Gets input stream for given file.
     *
     * @param filename file name
     *
     * @return input stream
     */
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

    /**
     * Return number of lines in given file.
     *
     * @param filename file name
     *
     * @return number of lines
     */
    public static int getNumberLines(String filename) {

        InputStream input = getInputStream(filename);
        try {

            byte[] buffer = new byte[1024];
            int result = 0;
            int numChrsRead = 0;
            while ((numChrsRead = input.read(buffer)) != -1) {
                for (int i = 0; i < numChrsRead; ++i) {
                    if (buffer[i] == '\n') {
                        ++result;
                    }
                }
            }
            return result;

        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalStateException("Utils: getNumberLines: Problem getting number lines: " + filename);
        } finally {
            try {
                input.close();
            } catch (Exception e) {
                // do nothing
            }
        }
    }

    public static String readLineSkipComments(BufferedReader br) throws IOException {
        String s = br.readLine();
        while ((s.startsWith("#"))) {
            s = br.readLine();
        }
        return s;
    }
    
    public static String getSystemTime(){
		return new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").
				format(Calendar.getInstance().getTime());
	}
}
