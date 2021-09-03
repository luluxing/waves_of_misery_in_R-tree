package xxl.core.indexStructures.rtrees;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class PropertyReader {

    private Properties m_prop;

    public String m_treeName;
    public int m_insertionSeed;
    public int m_primaryIndexSize;
    public String m_outputDir;
    public boolean m_bulkLoading;
    public String m_bulkLoadingFileName;
    public int m_numBatches;
    public boolean m_performRangeQuery;
    public boolean m_performknnQuery;
    public String m_compareMethod;
    public double m_bulkLoadingInitialCap;
    public boolean m_isUniform;
    public boolean m_unequalRandSplit;
    public double m_queryArea;
    public String m_queryFile;
    public int m_queryNumber;
    public boolean m_queryTailLatency;
    public boolean m_practicalRemedyLinear;
    public boolean m_practicalRemedyRandom;
    public boolean m_regularSplit;
    public int m_splitFreq = 0;
    public double m_minCap = 0;
    public int m_overflowPageNum = 1;
    public boolean m_fastVersion = false;
    public int m_sibling;
    public boolean m_reinsertion;
    public int m_batchSize;
    public boolean m_realData;
    public double m_xMin;
    public double m_xMax;
    public double m_yMin;
    public double m_yMax;
    public String m_batchInsertFile;
    public boolean m_isCSV;
    public double m_minMaxFactor;
    public int m_dimension;
    public int m_blockSize;
    public int m_bufferSize;
    public boolean m_soundRemedy;
    public int m_remedySeed;
    public boolean m_inMemory;

    static public <T> void wsm_print(T... ts) {
        System.out.print("[WSM_LOG] ");
        for (T t : ts)
            System.out.print(t);
        System.out.println("");
    }

    private <T> void print_prop(String prop, T val) {
        System.out.println(prop + " has value: " + val);
    }

    public void print_properties() {
        System.out.println("========= Core parameters of the tree =========");
        print_prop("Tree_name", m_treeName);
        print_prop("Min_max_factor", m_minMaxFactor);
        print_prop("Dimension", m_dimension);
        print_prop("Block_size", m_blockSize);
        print_prop("Buffer_size", m_bufferSize);

        System.out.println("========= Experiment details =========");
        print_prop("Primary_index_size", m_primaryIndexSize);
        print_prop("Output_dir", m_outputDir);
        print_prop("Load data file", m_bulkLoadingFileName);
        print_prop("Bulk_loading", m_bulkLoading);

        print_prop("Node_initial_capacity", m_bulkLoadingInitialCap);
        print_prop("Batch_size", m_batchSize);
        print_prop("Number_of_batches", m_numBatches);

        System.out.println("========= Querying =========");
        if (m_performRangeQuery) {
            print_prop("Query_number", m_queryNumber);
            print_prop("Query_area", m_queryArea);
            print_prop("Query_file", m_queryFile);
            if (m_queryTailLatency)
                print_prop("Query_tail_latency", m_queryTailLatency);
        }
        if (m_performknnQuery) {
            print_prop("Perform knn query", m_performknnQuery);
        }

        System.out.println("========= Optional experiment parameters =========");
        print_prop("Comparison method (hilbert, peano or x)", m_compareMethod);
        print_prop("Batch inserted data is uniform (otherwise normal)", m_isUniform);
        print_prop("Batch inserted seed", m_insertionSeed);

        System.out.println("========= Using real dataset =========");
        if (m_realData) {
            print_prop("Using real-world data", m_realData);
            print_prop("Real-world data file", m_batchInsertFile);
            print_prop("Batch-insertion file is a csv tile", m_isCSV);
        }

        System.out.println("========= Using in memory index =========");
        if (m_inMemory) {
            print_prop("In-memory index", m_inMemory);
        }

        System.out.println("========= Remedies =========");
        if (m_soundRemedy)
            print_prop("Sound remedy", m_soundRemedy);
        if (m_practicalRemedyLinear)
            print_prop("Practical remedy (Linear)", m_practicalRemedyLinear);
        if (m_practicalRemedyRandom)
            print_prop("Practical remedy (Random)", m_practicalRemedyRandom);
        if (m_unequalRandSplit)
            print_prop("Unequal random split", m_unequalRandSplit);
        if (m_regularSplit) {
            print_prop("Regular elective split", m_regularSplit);
            print_prop("split_freq", m_splitFreq);
            print_prop("min_capacity", m_minCap);
            print_prop("Overflow_page_num", m_overflowPageNum);
            print_prop("fast_version", m_fastVersion);
        }

        System.out.println("========= Tree specific modifications =========");
        if (m_treeName.startsWith("hilbert") && m_sibling != 1)
            print_prop("Hilbert R-tree deferred split cooperative sibling number", m_sibling);
        if (m_treeName.startsWith("rstar") && m_reinsertion)
            print_prop("R*-tree reinsertion", m_reinsertion);

        System.out.println("====================================");
    }

    public PropertyReader(String fileName) {
        try (InputStream input = new FileInputStream(fileName)) {
            m_prop = new Properties();
            m_prop.load(input);

            // Core parameters of the tree
            m_treeName = m_prop.getProperty("Tree_name");
            m_minMaxFactor = Double.parseDouble(m_prop.getProperty("Min_max_factor", "0.5"));
            m_dimension = Integer.parseInt(m_prop.getProperty("Dimension", "2"));
            m_blockSize = Integer.parseInt(m_prop.getProperty("Block_size", "15371"));
            m_bufferSize = Integer.parseInt(m_prop.getProperty("Buffer_size", "20"));

            // Mandatory parameters of the experiment
            m_primaryIndexSize = Integer.parseInt(m_prop.getProperty("Primary_index_size"));
            m_outputDir = m_prop.getProperty("Output_dir") + "/";
            m_bulkLoading = Boolean.parseBoolean(m_prop.getProperty("Bulk_loading"));
            m_bulkLoadingFileName = m_prop.getProperty("Bulk_loading_file");

            m_bulkLoadingInitialCap = Double.parseDouble(m_prop.getProperty("initial_capacity"));
            m_batchSize = Integer.parseInt(m_prop.getProperty("Batch_size", "10000"));
            m_numBatches = Integer.parseInt(m_prop.getProperty("Number_of_batches"));
            m_performRangeQuery = Boolean.parseBoolean(m_prop.getProperty("Perform_query"));
            if (m_performRangeQuery) {
                m_queryNumber = Integer.parseInt(m_prop.getProperty("Query_number"));
                m_queryArea = Double.parseDouble(m_prop.getProperty("Query_area"));
                m_queryFile = m_prop.getProperty("Query_file");
                m_queryTailLatency = Boolean.parseBoolean(m_prop.getProperty("Tail_latency", "false"));
            }
            m_performknnQuery = Boolean.parseBoolean(m_prop.getProperty("kNN_query", "false"));

            // Optional parameters
            m_compareMethod = m_prop.getProperty("compare", "hilbert");
            m_isUniform = Boolean.parseBoolean(m_prop.getProperty("Uniform_distribution", "true"));
            m_insertionSeed = Integer.parseInt(m_prop.getProperty("Insertion_seed", "43"));

            m_realData = Boolean.parseBoolean(m_prop.getProperty("Use_real_data", "false"));
            if (m_realData) {
                m_batchInsertFile = m_prop.getProperty("Batch_insertion_file");
                m_isCSV = Boolean.parseBoolean(m_prop.getProperty("is_csv", "false"));
                m_xMin = Double.parseDouble(m_prop.getProperty("x_min"));
                m_xMax = Double.parseDouble(m_prop.getProperty("x_max"));
                m_yMin = Double.parseDouble(m_prop.getProperty("y_min"));
                m_yMax = Double.parseDouble(m_prop.getProperty("y_max"));
            }
            m_inMemory = Boolean.parseBoolean(m_prop.getProperty("In_memory", "false"));

            // Sound remedy
            m_soundRemedy = Boolean.parseBoolean(m_prop.getProperty("SoundRemedy", "false"));
            m_remedySeed = Integer.parseInt(m_prop.getProperty("RemedySeed", "1"));

            // Practical remedy
            m_practicalRemedyLinear = Boolean.parseBoolean(m_prop.getProperty("Practical_remedy_linear", "false"));
            m_practicalRemedyRandom = Boolean.parseBoolean(m_prop.getProperty("Practical_remedy_random", "false"));

            // Number of cooperative siblings for Hilbert R-tree deferred split. (only 2 is supported)
            m_sibling = Integer.parseInt(m_prop.getProperty("defer", "1"));
            if (m_treeName.equalsIgnoreCase("hilbert") && m_sibling != 1) {
                m_treeName += "_defer";
            }

            // To enable re-insertion in R*-tree or not
            m_reinsertion = Boolean.parseBoolean(m_prop.getProperty("reinsertion", "false"));

            // Unequal random split (URS)
            m_unequalRandSplit = Boolean.parseBoolean(m_prop.getProperty("Unequal_random_split", "false"));

            // Regular elective split (RES)
            m_regularSplit = Boolean.parseBoolean(m_prop.getProperty("Regular_split", "false"));
            if (m_regularSplit) {
                m_splitFreq = Integer.parseInt(m_prop.getProperty("split_freq", "400"));
                m_overflowPageNum = Integer.parseInt(m_prop.getProperty("Overflow_page_num", "1"));
                m_fastVersion = Boolean.parseBoolean(m_prop.getProperty("fast_version", "true"));
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
