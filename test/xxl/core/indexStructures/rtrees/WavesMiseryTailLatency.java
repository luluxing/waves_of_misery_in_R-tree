/* XXL: The eXtensible and fleXible Library for data processing

Copyright (C) 2000-2011 Prof. Dr. Bernhard Seeger
                        Head of the Database Research Group
                        Department of Mathematics and Computer Science
                        University of Marburg
                        Germany

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library;  If not, see <http://www.gnu.org/licenses/>.

    http://code.google.com/p/xxl/

*/

package xxl.core.indexStructures.rtrees;

import xxl.core.collections.containers.Container;
import xxl.core.collections.containers.CounterContainer;
import xxl.core.collections.containers.io.BlockFileContainer;
import xxl.core.collections.containers.io.BufferedContainer;
import xxl.core.collections.containers.io.ConverterContainer;
import xxl.core.collections.containers.io.RawAccessContainer;
import xxl.core.collections.queues.Queue;
import xxl.core.collections.queues.io.BlockBasedQueue;
import xxl.core.cursors.Cursor;
import xxl.core.cursors.Cursors;
import xxl.core.cursors.sorters.MergeSorter;
import xxl.core.cursors.sources.io.FileInputCursor;
import xxl.core.functions.AbstractFunction;
import xxl.core.functions.Function;
import xxl.core.indexStructures.*;
import xxl.core.indexStructures.Tree.Query.Candidate;
import xxl.core.io.FilesystemOperations;
import xxl.core.io.LRUBuffer;
import xxl.core.io.LogFilesystemOperations;
import xxl.core.io.XXLFilesystem;
import xxl.core.io.converters.*;
import xxl.core.io.raw.RAMRawAccess;
import xxl.core.io.raw.RawAccess;
import xxl.core.spatial.SpaceFillingCurves;
import xxl.core.spatial.points.DoublePoint;
import xxl.core.spatial.rectangles.DoublePointRectangle;
import xxl.core.spatial.rectangles.Rectangle;
import xxl.core.spatial.rectangles.Rectangles;

import java.io.*;
import java.lang.management.ManagementFactory;
import java.util.*;

/**
 * Simple Test of the HilbertRTree with random DoublePoint ({@link DoublePoint})
 */
public class WavesMiseryTailLatency {

    static public int dimension = 2;

    /**Precision of the Hilbert Space filling curve*/
    static public final int FILLING_CURVE_PRECISION = 1<<20;

    /** Universe of the data space
     * ( DoublePoints uniformly distributed in rectangle with xmin = 0 , ymin = 0, xmax = 1 , ymax = 1)
     */
    static public DoublePointRectangle universe = new DoublePointRectangle(new double[]{0,0},
            new double[]{1.0,1.0});

    /**
     * Converter for the keys of the data. Keys are the java.long values(hilbert values of the MBRs middle points).
     */
    static public MeasuredConverter keyConverter = new MeasuredFixedSizeConverter<Long>(LongConverter.DEFAULT_INSTANCE);

    /**
     * Measured converter for DoublePoints
     */
    static public MeasuredConverter dataConverter = new MeasuredConverter(){
        public int getMaxObjectSize() {
            return dimension * 8;
        }
        public Object read(DataInput dataInput, Object object) throws IOException {
            return  ConvertableConverter.DEFAULT_INSTANCE.read(dataInput, new DoublePoint(dimension));
        }
        public void write(DataOutput dataOutput, Object object) throws IOException {
            ConvertableConverter.DEFAULT_INSTANCE.write(dataOutput, (DoublePoint)object);
        }
    };

    /**
     * Function creating a descriptor for a point.
     */
//    static Function<Object,Object> GET_DESCRIPTOR = new AbstractFunction<Object,Object>() {
//        @Override
//        public Object invoke (Object o) {
//            DoublePoint p = (DoublePoint)o;
//            return new DoublePointRectangle(p, p);
//        }
//    };

    /**
     * Function creating a MBR ({@link DoublePointRectangle})  for a given doublePoint object .
     */
    static Function getEntryMBR  = new AbstractFunction() {
        public Object invoke (Object o) {
            DoublePoint p = (DoublePoint)o;
            return new DoublePointRectangle(p, p);
        }
    };

    /**
     * A function that sorts the points according to their hilbert value.
     * This is needed for sort-based bulk insertion.
     */
    static Comparator COMPARE_HILBERT = new Comparator() {
        protected double uni[];
        protected double uniDeltas[];

        public int compare(Object o1, Object o2) {
            if (dimension!=2) throw new IllegalArgumentException();
            if (uni == null) {
                uni = (double[])universe.getCorner(false).getPoint();
                uniDeltas = universe.deltas();
            }
//			double leftBorders1[] = (double[])((Rectangle)((KPE)o1).getData()).getCorner(false).getPoint();
//			double leftBorders2[] = (double[])((Rectangle)((KPE)o2).getData()).getCorner(false).getPoint();
            double leftBorders1[] = (double[]) ((DoublePoint)o1).getPoint();
            double leftBorders2[] = (double[]) ((DoublePoint)o2).getPoint();

            double x1 = (leftBorders1[0]-uni[0])/uniDeltas[0];
            double y1 = (leftBorders1[1]-uni[1])/uniDeltas[1];
            double x2 = (leftBorders2[0]-uni[0])/uniDeltas[0];
            double y2 = (leftBorders2[1]-uni[1])/uniDeltas[1];

            long h1 = SpaceFillingCurves.hilbert2d((int) (x1*FILLING_CURVE_PRECISION),(int) (y1*FILLING_CURVE_PRECISION));
            long h2 = SpaceFillingCurves.hilbert2d((int) (x2*FILLING_CURVE_PRECISION),(int) (y2*FILLING_CURVE_PRECISION));

            return (h1<h2)?-1: ((h1==h2)?0:+1);
        }
    };

    /**
     *
     * A function that sorts the points according to their peano value.
     * This is needed for sort-based bulk insertion.
     */
    static Comparator COMPARE_PEANO = new Comparator() {
        protected double uni[];
        protected double uniDeltas[];

        public int compare(Object o1, Object o2) {
            if (dimension!=2) throw new IllegalArgumentException();
            if (uni == null) {
                uni = (double[])universe.getCorner(false).getPoint();
                uniDeltas = universe.deltas();
            }
            double leftBorders1[] = (double[]) ((DoublePoint)o1).getPoint();
            double leftBorders2[] = (double[]) ((DoublePoint)o2).getPoint();

            double x1 = (leftBorders1[0]-uni[0])/uniDeltas[0];
            double y1 = (leftBorders1[1]-uni[1])/uniDeltas[1];
            double x2 = (leftBorders2[0]-uni[0])/uniDeltas[0];
            double y2 = (leftBorders2[1]-uni[1])/uniDeltas[1];

            long h1 = SpaceFillingCurves.peano2d((int) (x1*FILLING_CURVE_PRECISION),(int) (y1*FILLING_CURVE_PRECISION));
            long h2 = SpaceFillingCurves.peano2d((int) (x2*FILLING_CURVE_PRECISION),(int) (y2*FILLING_CURVE_PRECISION));
            return (h1<h2) ? -1 : ( (h1==h2) ? 0 : 1 );
        }
    };

    /**
     *  A function that sorts the points according to x coordinate
     *  This is needed for sort-based bulk insertion.
     */
    static Comparator COMPARE = new Comparator() {
        public int compare(Object o1, Object o2) {
            double p1[] = (double[]) ((DoublePoint)o1).getPoint();
            double p2[] = (double[]) ((DoublePoint)o2).getPoint();
            return (p1[0]<p2[0]) ? -1 : ( (p1[0]==p2[0]) ? 0 : 1);
//			return ((Rectangle)(((KPE)o1).getData())).compareTo(((KPE)o2).getData());
        }
    };

    static void writeIntToFile(String fileName, ArrayList arr) throws FileNotFoundException {
        if (arr.size() == 0)
            return;
        PrintWriter pw = new PrintWriter(fileName);
        for (int i = 0; i < arr.size(); i++) {
            if (arr.get(i) instanceof Integer) {
                pw.printf("%d\n", arr.get(i));
            }
            else {
                ArrayList<Integer> l = (ArrayList<Integer>) arr.get(i);
                for (int j = 0; j < l.size(); j++)
                    pw.printf("%d,", l.get(j));
                pw.printf("\n");
            }
        }
        pw.close();
    }

    static void writeDoubleToFile(String fileName, ArrayList arr) throws FileNotFoundException {
        if (arr.size() == 0)
            return;
        PrintWriter pw = new PrintWriter(fileName);
        for (int i = 0; i < arr.size(); i++) {
            if (arr.get(i) instanceof Double) {
                pw.printf("%f\n", arr.get(i));
            }
            else {
                ArrayList<Double> l = (ArrayList<Double>) arr.get(i);
                for (int j = 0; j < l.size(); j++)
                    pw.printf("%f,", l.get(j));
                pw.printf("\n");
            }
        }
        pw.close();
    }

    /**
     * Returns a comparator which evaluates the distance of two candidate objects
     * to the specified <tt>queryObject</tt>. This comparator
     * is used for nearest neighbor queries and defines an order on candidate-
     * descriptors. With the help of a priority queue (Min-heap) and this
     * comparator the nearest neighbor query can be performed.
     *
     * @param queryObject a KPE to which the nearest neighbors should be determined
     * @return a comparator defining an order on candidate objects
     */
    public static Comparator getDistanceBasedComparator (DoublePoint queryObject) {
        final Rectangle query = new DoublePointRectangle(queryObject, queryObject);
        return new Comparator () {
            public int compare (Object candidate1, Object candidate2) {
                Rectangle r1 =( (HilbertRTree.ORSeparator) ((Candidate) candidate1).descriptor() ).getIndexEntryMBR();
                Rectangle r2 = ( (HilbertRTree.ORSeparator) ((Candidate) candidate2).descriptor() ).getIndexEntryMBR();
                double d1 = query.distance(r1, 2);
                double d2 = query.distance(r2, 2);
                return (d1<d2) ? -1 : ( (d1==d2) ? 0 : 1 );
            }
        };
    }

    /**
     * This function computes hilbert curve value of the MBRs middle point
     * There is no need to do normalization as Hilbert R-tree does this internally
     */
    static public Function getHilbertValue = new AbstractFunction(){

        public Object invoke(Object point){
            if (point instanceof DoublePoint){

                DoublePoint middlePoint = (DoublePoint)point;
                double x = middlePoint.getValue(0);// - uni[0]) / uniDeltas[0];
                double y = middlePoint.getValue(1);// - uni[1]) / uniDeltas[1];
                return new Long(SpaceFillingCurves.hilbert2d((int)
                        (x*FILLING_CURVE_PRECISION),(int) (y*FILLING_CURVE_PRECISION)));
            }
            throw new IllegalArgumentException();
        }
    };

    /**
     * Factory function that creates ORSeparator of the HilbertRTree
     */
    public static final Function createORSeparator =
            new AbstractFunction() {
                public Object invoke(Object key, Object mbr) {
                    return new DoublePointRectangleSep((Long)key, (DoublePointRectangle)mbr);
                }
            };
    /**
     *  Factory function that creates ORKeyRange of the HilbertRTree
     */
    public static final Function createORKeyRange =
            new AbstractFunction() {
                public Object invoke(List arguments) {
                    if(arguments.size() !=3 ) throw new IllegalArgumentException();
                    Long min = (Long)arguments.get(0);
                    Long max = (Long)arguments.get(1);
                    DoublePointRectangle entryMBR =  (DoublePointRectangle)arguments.get(2);
                    return new LongRange(min, max, entryMBR);
                }
            };

    /**
     * This class represents the ORSeparator of the HilbertRTree
     */
    public static class DoublePointRectangleSep extends HilbertRTree.ORSeparator{
        public DoublePointRectangleSep(Long separatorValue, DoublePointRectangle entryMBR) {
            super(separatorValue, entryMBR);
        }
        @Override
        public Object clone() {
            return new DoublePointRectangleSep(((Long)this.sepValue()).longValue(),
                    (DoublePointRectangle) ((Descriptor) this.entryMBR).clone()) ;
        }
    }

    /**
     * This class represents the ORKeyRange of the HilbertRTree
     */
    public static class LongRange extends HilbertRTree.ORKeyRange{
        public LongRange(Long min, Long max, DoublePointRectangle entryMBR) {
            super(min, max, entryMBR);
        }

        @Override
        public Object clone() {
            return new LongRange(((Long)this.sepValue).longValue(), ((Long)this.maxBound).longValue(),
                    (DoublePointRectangle) ((Descriptor) this.entryMBR).clone());
        }
    }

    /**
     * This method saves the parameters of the HilbertRTree, this allows to make the tree persistent
     * @param tree HilbertRTree
     * @param pfad
     * @throws IOException
     */
    public static void saveParams(HilbertRTree tree, String pfad, Converter keyKonverter) throws IOException{
        Tree.IndexEntry root = tree.rootEntry();
        HilbertRTree.ORKeyRange desc = (HilbertRTree.ORKeyRange)tree.rootDescriptor();
        int height =  tree.height();
        Long id = (Long)root.id();
        Rectangle rec = (Rectangle)desc.getIndexEntryMBR();
        Object minBound = desc.minBound();
        Object maxBound = desc.maxBound();
        DataOutputStream outPut = new DataOutputStream(new FileOutputStream(new File(pfad)));
        LongConverter.DEFAULT_INSTANCE.write(outPut, id);
        IntegerConverter.DEFAULT_INSTANCE.writeInt(outPut, height);
        keyKonverter.write(outPut, minBound);
        keyKonverter.write(outPut, maxBound);
        ConvertableConverter.DEFAULT_INSTANCE.write(outPut, (DoublePointRectangle)rec);
        outPut.close();
        ((Container)tree.getContainer.invoke()).flush();
        ((Container)tree.getContainer.invoke()).close();
    }

    private static double computeAvg(ArrayList<Long> list) {
        double s = 0.0;
        double base = 1.0 * list.get(0);
        for (int i = 0; i < list.size(); i++) {
            s += list.get(i) - base;
        }
        double avg = 1.0 * s / list.size() + base;
        return avg;
    }

    /**
     * Starts simpleHilbertRTree Test
     * @param args path To HilbertRTree
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

        if (args.length != 1) {
            System.out.println("Usage: Experiment.properties");
            return;
        }

        PropertyReader reader = new PropertyReader(args[0]);
        reader.print_properties();

        if (reader.m_realData) {
            universe = new DoublePointRectangle(new double[]{reader.m_xMin, reader.m_yMin},
                    new double[]{reader.m_xMax, reader.m_yMax});
        }
        PropertyReader.wsm_print("Universe is " + universe);

        Comparator<DoublePoint> sFCComparator;
        if (reader.m_compareMethod.equalsIgnoreCase("x"))
            sFCComparator = COMPARE;
        else if (reader.m_compareMethod.equalsIgnoreCase("hilbert"))
            sFCComparator = COMPARE_HILBERT;
        else if (reader.m_compareMethod.equalsIgnoreCase("peano"))
            sFCComparator = COMPARE_PEANO;
        else {
            PropertyReader.wsm_print("Compare method not found. Accept only: x, hilbert, peano.");
            return;
        }

        PropertyReader.wsm_print("Init new tree");
        String treePath = reader.m_outputDir + "/" + reader.m_treeName + reader.m_primaryIndexSize;

        /** Create a in-memory index or disk-based index */
        XXLFilesystem fs = new XXLFilesystem(true);

        FilesystemOperations fso = fs.getFilesystemOperations();
        fso = new LogFilesystemOperations(fso, System.out, true);

        Container rawContainer;
        if (reader.m_inMemory) {
            // Number of blocks is hard-coded
            RawAccess ra = new RAMRawAccess(1048576, reader.m_blockSize); // 4GB
//            RawAccess ra = new RAMRawAccess(10576, reader.m_blockSize);
            rawContainer = new RawAccessContainer(ra, 100);
        } else {
            rawContainer = new BlockFileContainer(treePath, reader.m_blockSize, fso);
        }

        Container ioContainer = new CounterContainer(rawContainer);

        // size of a data entry = size of DoublePoint (1 double per dimension)
        int dataSize = reader.m_dimension * 8;
        // size of a descriptor = size of DoublePointRectangle (2 doubles per dimension)
        int descriptorSize = reader.m_dimension * 2 * 8;
        Tree tree = null;
//        CounterContainer lowerCounterContainer;
        Container treeContainer;
        FilteredCounterContainer container;
//        FilteredCounterContainer lowerLeafCounter, upperLeafCounter;
        if (reader.m_treeName.startsWith("hilbert")) {
            tree = new HilbertRTree(reader.m_blockSize, universe, reader.m_minMaxFactor, reader.m_sibling);
            Container converterContainer = new ConverterContainer(ioContainer, ((HilbertRTree)tree).nodeConverter());

            treeContainer = new BufferedContainer(converterContainer,
                                                        new LRUBuffer(reader.m_bufferSize), true);

            // unbuffered container that counts the access to the RTree
            container = new FilteredCounterContainer(treeContainer,
                                                            FilteredCounterContainer.TREE_LEAF_NODE_COUNTER_FUNCTION);
            ((HilbertRTree)tree).initialize(getHilbertValue , getEntryMBR, container,
                        keyConverter,
                        dataConverter,
                        createORSeparator,
                        createORKeyRange);

        } else {
            if (reader.m_treeName.equalsIgnoreCase("rstar"))
                tree = new RTree(reader.m_reinsertion);
            else if (reader.m_treeName.equalsIgnoreCase("linear"))
                tree = new LinearRTree();
            else if (reader.m_treeName.equalsIgnoreCase("quadratic"))
                tree = new QuadraticRTree();
            else if (reader.m_treeName.equalsIgnoreCase("rrstar"))
                tree = new RevisedRTree(dimension);
            else {
                PropertyReader.wsm_print("Tree name not found. Accept only: hilbert, rstar, linear, quadratic, rrstar.");
                return;
            }

            Converter converter = ((RTree)tree).nodeConverter(new ConvertableConverter(
                                                        DoublePoint.factoryFunctionDoublePoint(reader.m_dimension)),
                                                        reader.m_dimension);
            ConverterContainer converterContainer = new ConverterContainer(rawContainer, converter);
            treeContainer = new BufferedContainer(converterContainer,
                                                        new LRUBuffer(reader.m_bufferSize), true);

            // unbuffered container that counts the access to the RTree
            container = new FilteredCounterContainer(treeContainer,
                                                        FilteredCounterContainer.TREE_LEAF_NODE_COUNTER_FUNCTION);

            ((RTree)tree).initialize(getEntryMBR, container, reader.m_blockSize, dataSize,
                        descriptorSize, reader.m_minMaxFactor);

        }

        tree.outDir = reader.m_outputDir;
        tree.unequalRandomSplit = reader.m_unequalRandSplit;
        tree.regularElectiveSplit = reader.m_regularSplit;
        if (tree.regularElectiveSplit) {
            tree.splitFrequency = reader.m_splitFreq;
//            tree.minElectiveSplit = reader.m_minCap;
            tree.fastElectiveSplit = reader.m_fastVersion;
            tree.overflowPageNum = reader.m_overflowPageNum;
        }

        long begin = System.nanoTime();

        /** Bulk-loading */
        PropertyReader.wsm_print("Construct the tree by bulk loading");

        Converter<DoublePoint> objectConverter =
                new ConvertableConverter<DoublePoint>(DoublePoint.factoryFunctionDoublePoint(reader.m_dimension));
        Cursor unsortedCursor = new FileInputCursor<DoublePoint>(objectConverter, new File(reader.m_bulkLoadingFileName));
        Container sortedRunsContainer = new BlockFileContainer(
                reader.m_outputDir + "sortedRuns_" + reader.m_treeName, reader.m_blockSize);
        Function<Function<?, Integer>, Queue<DoublePoint>> queueFactoryFunction =
                BlockBasedQueue.createBlockBasedQueueFunctionForMergeSorter(sortedRunsContainer,
                        reader.m_blockSize, objectConverter);
        // with this value we provide an available memory for initial run generation;
        int memorySizeForRuns = 1024 * 1024 * 10;
        // here we provide how much memory is needed for last merge
        int memoryLastRuns = memorySizeForRuns;
        Iterator<DoublePoint> sortedPointIterator = new MergeSorter<DoublePoint>(unsortedCursor, sFCComparator, dataSize, memorySizeForRuns,
                memoryLastRuns, queueFactoryFunction, false);

        if (reader.m_treeName.startsWith("hilbert")) {
            // Assuming ascending order only (for deferred splitting)
            if (reader.m_soundRemedy)
                new HilbertRTreeBulkLoading((HilbertRTree) tree, sortedPointIterator,
                        reader.m_remedySeed, reader.m_primaryIndexSize);
            else if (reader.m_practicalRemedyLinear)
                new HilbertRTreeBulkLoading((HilbertRTree) tree, sortedPointIterator,
                        reader.m_bulkLoadingInitialCap, reader.m_remedySeed, true);
            else if (reader.m_practicalRemedyRandom)
                new HilbertRTreeBulkLoading((HilbertRTree) tree, sortedPointIterator,
                        reader.m_bulkLoadingInitialCap, reader.m_remedySeed, false);
            else
                new HilbertRTreeBulkLoading((HilbertRTree) tree, sortedPointIterator,
                        reader.m_bulkLoadingInitialCap);
        } else {
            if (reader.m_soundRemedy)
                new SortBasedBulkLoading((RTree) tree, sortedPointIterator, tree.determineContainer,
                        reader.m_remedySeed, reader.m_primaryIndexSize, reader.m_outputDir);
            else if (reader.m_practicalRemedyLinear)
                new SortBasedBulkLoading((RTree) tree, sortedPointIterator, tree.determineContainer,
                        reader.m_remedySeed, reader.m_bulkLoadingInitialCap, reader.m_outputDir, true);
            else if (reader.m_practicalRemedyRandom)
                new SortBasedBulkLoading((RTree) tree, sortedPointIterator, tree.determineContainer,
                        reader.m_remedySeed, reader.m_bulkLoadingInitialCap, reader.m_outputDir, false);
            else
                new SortBasedBulkLoading((RTree) tree, sortedPointIterator, tree.determineContainer,
                        reader.m_bulkLoadingInitialCap);
        }

        long done = System.nanoTime();

        PropertyReader.wsm_print("Loading finishes in ", (done - begin) / (1000.0*1000*1000*60), " min.");
        PropertyReader.wsm_print("Total data in leaf nodes is ", tree.numOfDataInLeafNodes);
        PropertyReader.wsm_print("Total leaf node capacity is ", tree.totalLeafNodeCapacity);
        PropertyReader.wsm_print("Initial leaf node utilization is ", tree.numOfDataInLeafNodes/tree.totalLeafNodeCapacity);
        PropertyReader.wsm_print("Single leaf node capacity is ", tree.singleLeafNodeCapacity);
        PropertyReader.wsm_print("Single index node capacity is ", tree.singleIndexNodeCapacity);
        PropertyReader.wsm_print("Total leaf node number is ", tree.leafNodeNum);

        container.reset();

        if (reader.m_bufferSize > 0) {
            container.flush();
        }

        /** Load Queries */
        List<DoublePointRectangle> queries = new ArrayList<>();
        if (reader.m_performRangeQuery) {

            if (reader.m_queryTailLatency) {
                Converter<DoublePointRectangle> rectangleConverter = new ConvertableConverter<DoublePointRectangle>(
                        Rectangles.factoryFunctionDoublePointRectangle(reader.m_dimension));
                String[] queryFileNames = {"query_0.000001_num_1k", "query_0.00001_num_1k", "query_0.0001_num_1k",
                                        "query_0.001_num_1k", "query_0.01_num_1k", "query_0.05_num_1k"};
                for (int i = 0; i < queryFileNames.length; i++) {
                    String queryPath = reader.m_queryFile + "/" + queryFileNames[i];
                    PropertyReader.wsm_print("Load queries from " + queryPath);
                    int number = reader.m_queryNumber;
                    Cursor queryCursor = new FileInputCursor<DoublePointRectangle>(rectangleConverter, new File(queryPath));
                    while (queryCursor.hasNext() && number > 0) {
                        Object obj = queryCursor.next();
                        queries.add((DoublePointRectangle) obj);
                        number--;
                    }
                    queryCursor.close();
                }
            }
        }


        /*********************************************************************/
        /*        			     INSERTING RANDOM DATA                       */
        /*********************************************************************/
        PropertyReader.wsm_print("Begin inserting data and counting splits");

        begin = System.nanoTime();

        if (!reader.m_isUniform)
            PropertyReader.wsm_print("Using normal distribtuion");

        BufferedReader br = null;
        if (reader.m_realData) {
            br = new BufferedReader(new FileReader(reader.m_batchInsertFile));
        }
        String delimit = reader.m_isCSV ? "," : " ";

        Random random = new Random(reader.m_insertionSeed);

        ArrayList<Integer> bufferEvictionList = new ArrayList<>();
        ArrayList<Double> bufferUtilList = new ArrayList<>();

        double nanoToMilliSeconds = 1000 * 1000;

        double[] point;
        long start;
        long end;

        int leafSplitCount, nonLeafSplitCount;
        int insertionLeafAccess, insertionNonLeafAccess;
        double insertionTime = 0, queryCPUTime = 0;
        double avgLeafUtilRatio, avgLeafDataOverLeafNum, avgTotalDataOverLeafNum, avgTotalDataOverAllPageNum;
        int rootLevel;
        ArrayList<Long> leafNodeNum = new ArrayList<>();

        int numBatches = reader.m_numBatches;
        int batchSize = reader.m_batchSize;
        int dimension = reader.m_dimension;
        boolean isUniform = reader.m_isUniform;
        String line;

        double totalQt = 0, bufferUtil = 0, totalIOtime = 0;
        int bufferDisplace = 0, preInsertBufferDisplace = 0, preQueryBufferDisplace = 0;
        boolean end_early = false, first_query = false;

        String[] outputFileNames = {
                "leafSplits", "nonLeafSplits", "leafSplits2", "nonLeafSplits2", // 0 - 3
                "avgLeafRatio", "avgLeafDataOverLeafNumList", // 4 - 5
                "avgTotalDataOverLeafNumList", "avgTotalDataOverAllPageNumList", // 6 - 7
                "insertTime", "insertionLeafAccess", "insertionNonLeafAccess", // 8 - 10
                "insertionLeafAccess2", "insertionNonLeafAccess2", // 11 - 12
                "treeHeight", "leafNodeNum", // 13 - 14
                "queryResultSize", "queryLeafAccessList", "queryNonLeafAccessList", // 15 - 17
                "queryTime", "queryAvgIO", "queryAvgLeaf", // 18 - 20
                    "bufferEviction", "bufferUtil",  // 21 - 22
                "queryBufferEviction", // 23
                "queryCPUTime", "queryIOTime"}; // 24, 25
        ArrayList<PrintWriter> pwArray = new ArrayList<>();
        for (int i = 0; i < outputFileNames.length; i++) {
            PrintWriter tempPW = new PrintWriter(reader.m_outputDir + reader.m_treeName + "_" +
                                                                                    outputFileNames[i] + ".txt");
            pwArray.add(tempPW);
        }


        ((BufferedContainer) treeContainer).flush();
        ((BufferedContainer) treeContainer).clearBuffers();


        PropertyReader.wsm_print("Performing warm-up queries");

        Cursor preResults;
        if (reader.m_treeName.startsWith("hilbert"))
            preResults = ((HilbertRTree)tree).queryOR();
        else {
            int[] tmp = {0,0};
            preResults = ((RTree) tree).query(((RTree) tree).rootDescriptor(), 0, tmp);
        }
        Cursors.count(preResults);


        for (int j = 0; j < numBatches; j++) {
            if (end_early)
                break;
            if (j % 50 == 0)
                PropertyReader.wsm_print("Currently running " + j + "th batch.");

            bufferUtil = 0.0;
            preInsertBufferDisplace = ((BufferedContainer) treeContainer).getBufferDisplaceNum();

            tree.curBatchNum++;

            leafSplitCount = 0;
            nonLeafSplitCount = 0;

            insertionLeafAccess = 0;
            insertionNonLeafAccess = 0;
            insertionTime = 0;
            avgLeafUtilRatio = 0;
            avgLeafDataOverLeafNum = 0;
            avgTotalDataOverLeafNum = 0;
            avgTotalDataOverAllPageNum = 0;
            rootLevel = 0;
//            leafNodeNum = 0;
            leafNodeNum.clear();

            container.reset();
//            container.flush();
            for (int k = 0; k < batchSize; k++) {
                // Generate a new point
                point = new double[dimension];
                if (reader.m_realData) {
                    if ((line = br.readLine()) != null) {
                        String[] arrOfStr = line.split(delimit, 2);
                        point = new double[]{Double.parseDouble(arrOfStr[0]), Double.parseDouble(arrOfStr[1])};
                    } else {
                        PropertyReader.wsm_print("Real world data not enough");
                        end_early = true;
                        break;
                    }
                } else {
                    for (int i = 0; i < dimension; i++) {
                        if (isUniform)
                            point[i] = random.nextDouble();
                        else
                            point[i] = random.nextGaussian() * 0.5 / 3 + 0.5;
                    }
                }

                // Insert this point
                start = System.nanoTime();
                tree.insert(new DoublePoint(point));
                end = System.nanoTime();

                insertionTime += (end - start);
                leafSplitCount += tree.leafSplitCount;
                nonLeafSplitCount += tree.nonLeafSplitCount;

                ((BufferedContainer) treeContainer).getBufferUtilization();
                double temp = ((BufferedContainer) treeContainer).getBufferStoredLeafEntry() * dataSize;
                temp += ((BufferedContainer) treeContainer).getBufferStoredIndexEntry() * descriptorSize;
                bufferUtil += temp / (reader.m_bufferSize * reader.m_blockSize);


                insertionLeafAccess += tree.insertionLeafAccess;
                insertionNonLeafAccess += tree.insertionNonLeafAccess;

                rootLevel += tree.rootEntry().level();
                leafNodeNum.add(tree.leafNodeNum);

                avgLeafUtilRatio += tree.numOfDataInLeafNodes / tree.totalLeafNodeCapacity;
                avgLeafDataOverLeafNum += (tree.numOfDataInLeafNodes - tree.totalOverflowEntryNum)
                        / (1.0 * tree.singleLeafNodeCapacity * tree.leafNodeNum);
                avgTotalDataOverLeafNum += tree.numOfDataInLeafNodes / (1.0 * tree.singleLeafNodeCapacity * tree.leafNodeNum);
                avgTotalDataOverAllPageNum += 1.0 * tree.numOfDataInLeafNodes /
                        (tree.singleLeafNodeCapacity *(tree.leafNodeNum + tree.totalOverflowPageNum));
            }


            pwArray.get(0).printf("%d\n", leafSplitCount);
            pwArray.get(1).printf("%d\n", nonLeafSplitCount);
            pwArray.get(2).printf("%d\n", container.insertPredicates);
            pwArray.get(3).printf("%d\n", container.inserts - container.insertPredicates);
            avgLeafUtilRatio /= batchSize;
            pwArray.get(4).printf("%f\n", avgLeafUtilRatio);
            avgLeafDataOverLeafNum /= batchSize;
            pwArray.get(5).printf("%f\n", avgLeafDataOverLeafNum);
            avgTotalDataOverLeafNum /= batchSize;
            pwArray.get(6).printf("%f\n", avgTotalDataOverLeafNum);
            avgTotalDataOverAllPageNum /=batchSize;
            pwArray.get(7).printf("%f\n", avgTotalDataOverAllPageNum);
            pwArray.get(8).printf("%f\n", insertionTime / (batchSize * nanoToMilliSeconds));
            pwArray.get(9).printf("%d\n", insertionLeafAccess);
            pwArray.get(10).printf("%d\n", insertionNonLeafAccess);
            pwArray.get(11).printf("%d\n", container.getsPredicates);
            pwArray.get(12).printf("%d\n", container.gets - container.getsPredicates);
            rootLevel /= batchSize;
            pwArray.get(13).printf("%d\n", rootLevel);
            pwArray.get(14).printf("%f\n", computeAvg(leafNodeNum));

            bufferDisplace = ((BufferedContainer) treeContainer).getBufferDisplaceNum() - preInsertBufferDisplace;
            preQueryBufferDisplace = ((BufferedContainer) treeContainer).getBufferDisplaceNum();
            pwArray.get(21).printf("%d\n", bufferDisplace);
            pwArray.get(22).printf("%f\n", bufferUtil / batchSize);

            for (int i = 0; i < outputFileNames.length; i++) {
                pwArray.get(i).flush();
            }

            if (!reader.m_performRangeQuery)
                continue;

            int resultSize = 0;
            double qt = 0;
            int writeToDisk = 0;
            container.reset();
            container.flush();

            if (reader.m_queryTailLatency) {
                int number = reader.m_queryNumber;
                for (int q = 0; q < queries.size(); q++) {
                    ((BlockFileContainer) rawContainer).ioTime = 0;
                    int[] nodeAccess = {0, 0};
                    preQueryBufferDisplace = ((BufferedContainer) treeContainer).getBufferDisplaceNum();

                    Cursor results;
                    if (reader.m_treeName.startsWith("hilbert"))
                        results = ((HilbertRTree)tree).queryOR(queries.get(q), 0, nodeAccess);
                    else
                        results = ((RTree)tree).query(queries.get(q), 0, nodeAccess);
                    start = ManagementFactory.getThreadMXBean().getCurrentThreadCpuTime();
//                start = System.nanoTime();
                    resultSize += Cursors.count(results);
//                end = System.nanoTime();
                    end = ManagementFactory.getThreadMXBean().getCurrentThreadCpuTime();
                    results.close();
                    qt += end - start;
                    qt += ((BlockFileContainer) rawContainer).ioTime;

                    writeToDisk = ((BufferedContainer) treeContainer).getBufferDisplaceNum() - preQueryBufferDisplace;
                    number--;
                    String delim = number == 0 ? ";" : ",";
                    pwArray.get(23).printf("%d%s", writeToDisk, delim);
                    pwArray.get(24).printf("%f%s", (end - start) / nanoToMilliSeconds, delim);
                    pwArray.get(25).printf("%f%s", ((BlockFileContainer) rawContainer).ioTime / nanoToMilliSeconds, delim);
                    pwArray.get(19).printf("%d%s", container.gets, delim);
                    pwArray.get(20).printf("%d%s", container.getsPredicates, delim);
                    if (number == 0)
                        number = reader.m_queryNumber;
                }
            }

            pwArray.get(15).printf("%d\n", resultSize);
//            pwArray.get(16).printf("\n");
//            pwArray.get(17).printf("\n");
            pwArray.get(18).printf("%f\n", qt / nanoToMilliSeconds);
            pwArray.get(23).printf("\n");
            pwArray.get(24).printf("\n");
            pwArray.get(25).printf("\n");
            pwArray.get(19).printf("\n");
            pwArray.get(20).printf("\n");

            for (int i = 0; i < outputFileNames.length; i++) {
                pwArray.get(i).flush();
            }

        }


        done = System.nanoTime();
        PropertyReader.wsm_print("Done in ", (done - begin) / (1000.0*1000*1000*60), " min.");

    }

}
