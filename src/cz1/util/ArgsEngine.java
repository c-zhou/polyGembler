/*
 * Args - A Reusable Solution for Command Line Arguments Parsing in Java
 *
 * Copyright 2008 Adarsh Ramamurthy
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Web site: http://www.adarshr.com/papers/args
 */
package cz1.util;

import java.util.HashMap;
import java.util.Map;

/**
 * A reusable solution for Command Line Arguments Parsing in Java.
 * <p>
 * The following shows how to use this API.
 *
 * <pre>
 * // Initiate the arguments engine.
 * ArgsEngine engine = new ArgsEngine();
 *
 * // Configure the switches/options. Use true for valued options
 * engine.add("-q", "--quiet");
 * engine.add("-o", "--redirect-output", true);
 * engine.add("-h", "--help");
 *
 * // Perform the parsing.
 * engine.parse(args);
 *
 * // Start fetching states of switches
 * boolean quiet = engine.getBoolean("-q");
 *
 * if(engine.getBoolean("-o"))
 * {
 *     // For valued options, use getString.
 *     String redir = engine.getString("-o");
 * }
 * </pre>
 *
 * This has been modified from the original.
 *
 * @author Adarsh Ramamurthy
 *
 * @version 1.0, 12th May 2008
 */
public class ArgsEngine {

    /**
     * Holds all options configured.
     */
    private Map<String, Option> options = new HashMap<String, Option>();
    /**
     * Indicates if the parse method has been called.
     */
    private boolean parseCalled = false;

    /**
     * Configures the engine by adding options.
     *
     * @param shortForm the short form of an option. Example: "-h".
     *
     * @param longForm the long form of an option. Example: "--help".
     */
    public void add(String shortForm, String longForm) {
        this.add(shortForm, longForm, false);
    }

    /**
     * Configures the engine by adding options.
     *
     * @param shortForm the short form of an option. Example: "-h".
     *
     * @param longForm the long form of an option. Example: "--help".
     *
     * @param valued indicates if the option expects a value. The next argument
     * will be considered as the value for this option.
     */
    public void add(String shortForm, String longForm, boolean valued) {
        Option option = new Option(shortForm, longForm, valued);

        this.options.put(shortForm, option);
        this.options.put(longForm, option);
    }

    /**
     * Parses the input command line arguments. The result of parsing will
     * be stored in the current instance of <tt>ArgsEngine</tt>.
     * <p>Operations {@link #getBoolean(String)}, {@link #getNonOptions()}
     * and {@link #getString(String)} can be then used to extract the values.
     *
     * @param args the command line arguments.
     */
    public void parse(String[] args) {
        this.parseCalled = true;

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];

            // An option.
            if ((arg.startsWith("--") || arg.startsWith("-"))
                    && this.options.containsKey(arg)) {
                Option option = this.options.get(arg);

                if (option.isValued()) {
                    // This is the last arg.
                    if (i + 1 >= args.length) {
                        throw new RuntimeException("Value required for option "
                                + arg);
                    }

                    option.setValue(args[++i]);
                } else {
                    option.setValue("not-null");
                }
            } else {
                throw new RuntimeException("Unrecognized option " + arg);
            }
        }
    }

    /**
     * Gets the value for a valued option.
     *
     * @param key the valued option key. Either short form or long form
     * can be input.
     *
     * @return the value for the option if it's valued, <tt>null</tt>
     * otherwise.
     */
    public String getString(String key) {
        if (!this.parseCalled) {
            throw new IllegalStateException("Method parse not invoked");
        }

        Option option = this.options.get(key);

        // If option is non null and valued, return the value, null otherwise.
        return option != null
                ? (option.isValued() ? option.getValue() : null) : null;
    }

    /**
     * Gets the option.
     *
     * @param key the option's short form or long form name.
     *
     * @return <tt>true</tt> if the option is found in the args parsed,
     * <tt>false</tt> otherwise.
     */
    public boolean getBoolean(String key) {
        if (!this.parseCalled) {
            throw new IllegalStateException("Method parse not invoked");
        }

        Option option = this.options.get(key);

        return option != null ? option.getValue() != null : false;
    }

    /**
     * An object for representing and manipulating with the options.
     * <p>
     * This class is internal to the <tt>ArgsEngine</tt>.
     *
     * @author Adarsh Ramamurthy
     *
     * @version 1.0, 12th April 2008
     */
    private static class Option {

        /**
         * Short form of the option.
         */
        private String shortForm;
        /**
         * Long form of the option.
         */
        private String longForm;
        /**
         * Indicates if the option is valued.
         */
        private boolean valued;
        /**
         * The value of a valued option.
         */
        private String value;

        /**
         * Constructs an instance of <tt>Option</tt> taking the short and
         * long forms provided.
         *
         * @param shortForm the short form.
         *
         * @param longForm the long form.
         */
        public Option(String shortForm, String longForm) {
            this.shortForm = shortForm;
            this.longForm = longForm;
        }

        /**
         * Constructs an instance of <tt>Option</tt> taking the short and
         * long forms provided.
         *
         * @param shortForm the short form.
         *
         * @param longForm the long form.
         */
        public Option(String shortForm, String longForm, boolean valued) {
            this.shortForm = shortForm;
            this.longForm = longForm;
            this.valued = valued;
        }

        /**
         * Gets the short form name.
         *
         * @return the short form name.
         */
        public String getShortForm() {
            return this.shortForm;
        }

        /**
         * Gets the long form name.
         *
         * @return the long form name.
         */
        public String getLongForm() {
            return this.longForm;
        }

        /**
         * Tells if the option is valued.
         *
         * @return <tt>true</tt> if valued, <tt>false</tt> otherwise.
         */
        public boolean isValued() {
            return this.valued;
        }

        /**
         * Gets the value.
         *
         * @return the value.
         */
        public String getValue() {
            return this.value;
        }

        /**
         * Sets the value.
         *
         * @param value the value to set.
         */
        public void setValue(String value) {
            this.value = value;
        }
    }
}