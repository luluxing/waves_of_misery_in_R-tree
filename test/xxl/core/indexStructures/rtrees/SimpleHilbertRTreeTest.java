//
// Source code recreated from a .class file by IntelliJ IDEA
// (powered by FernFlower decompiler)
//

package xxl.core.indexStructures.rtrees;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Properties;
import java.util.Random;
import java.util.Scanner;
import xxl.core.collections.containers.Container;
import xxl.core.collections.containers.CounterContainer;
import xxl.core.collections.containers.io.BlockFileContainer;
import xxl.core.collections.containers.io.BufferedContainer;
import xxl.core.collections.containers.io.ConverterContainer;
import xxl.core.cursors.Cursor;
import xxl.core.cursors.sorters.MergeSorter;
import xxl.core.cursors.sources.io.FileInputCursor;
import xxl.core.functions.AbstractFunction;
import xxl.core.functions.Function;
import xxl.core.indexStructures.HilbertRTree;
import xxl.core.indexStructures.HilbertRTreeBulkLoading;
import xxl.core.indexStructures.HilbertRTree.ORKeyRange;
import xxl.core.indexStructures.HilbertRTree.ORSeparator;
import xxl.core.indexStructures.Tree.IndexEntry;
import xxl.core.indexStructures.Tree.Query.Candidate;
import xxl.core.io.LRUBuffer;
import xxl.core.io.converters.ConvertableConverter;
import xxl.core.io.converters.Converter;
import xxl.core.io.converters.IntegerConverter;
import xxl.core.io.converters.LongConverter;
import xxl.core.io.converters.MeasuredConverter;
import xxl.core.io.converters.MeasuredFixedSizeConverter;
import xxl.core.spatial.SpaceFillingCurves;
import xxl.core.spatial.points.DoublePoint;
import xxl.core.spatial.rectangles.DoublePointRectangle;
import xxl.core.spatial.rectangles.Rectangle;

public class SimpleHilbertRTreeTest {
    public static String rectFile = "rect_file.rec";
    public static String pointsFile = "point_file.pt";
    public static int batchSize = 10000;
    public static boolean soundRemedy = false;
    public static int remedySeed = 1;
    public static double minMaxFactor = 0.5D;
    public static int iters = 7000;
    public static int dimension = 2;
    public static int blockSize = 15371;
    public static int bufferSize = 0;
    public static final int FILLING_CURVE_PRECISION = 1048576;
    public static DoublePointRectangle universe = new DoublePointRectangle(new double[]{0.0D, 0.0D}, new double[]{1.0D, 1.0D});
    public static MeasuredConverter keyConverter;
    public static MeasuredConverter dataConverter;
    static Function getEntryMBR;
    public static Function getHilbertValue;
    public static final Function createORSeparator;
    public static final Function createORKeyRange;
    static Comparator COMPARE_HILBERT;
    static Comparator COMPARE_PEANO;
    static Comparator COMPARE;

    public SimpleHilbertRTreeTest() {
    }

    public static Comparator getDistanceBasedComparator(DoublePoint queryObject) {
        final Rectangle query = new DoublePointRectangle(queryObject, queryObject);
        return new Comparator() {
            public int compare(Object candidate1, Object candidate2) {
                Rectangle r1 = ((ORSeparator)((Candidate)candidate1).descriptor()).getIndexEntryMBR();
                Rectangle r2 = ((ORSeparator)((Candidate)candidate2).descriptor()).getIndexEntryMBR();
                double d1 = query.distance(r1, 2);
                double d2 = query.distance(r2, 2);
                return d1 < d2 ? -1 : (d1 == d2 ? 0 : 1);
            }
        };
    }

    public static void saveParams(HilbertRTree tree, String pfad, Converter keyKonverter) throws IOException {
        IndexEntry root = tree.rootEntry();
        ORKeyRange desc = (ORKeyRange)tree.rootDescriptor();
        int height = tree.height();
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

    static void writeIntToFile(String fileName, ArrayList arr) throws FileNotFoundException {
        if (arr.size() != 0) {
            PrintWriter pw = new PrintWriter(fileName);

            for(int i = 0; i < arr.size(); ++i) {
                if (arr.get(i) instanceof Integer) {
                    pw.printf("%d\n", arr.get(i));
                } else {
                    ArrayList<Integer> l = (ArrayList)arr.get(i);

                    for(int j = 0; j < l.size(); ++j) {
                        pw.printf("%d,", l.get(j));
                    }

                    pw.printf("\n");
                }
            }

            pw.close();
        }
    }

    static void writeDoubleToFile(String fileName, ArrayList arr) throws FileNotFoundException {
        if (arr.size() != 0) {
            PrintWriter pw = new PrintWriter(fileName);

            for(int i = 0; i < arr.size(); ++i) {
                if (arr.get(i) instanceof Double) {
                    pw.printf("%f\n", arr.get(i));
                } else {
                    ArrayList<Double> l = (ArrayList)arr.get(i);

                    for(int j = 0; j < l.size(); ++j) {
                        pw.printf("%f,", l.get(j));
                    }

                    pw.printf("\n");
                }
            }

            pw.close();
        }
    }

    static List generateFixedRangeQueries(double queryArea, int queryNum, Random random) {
        List<DoublePointRectangle> res = new ArrayList();
        if (queryArea < 0.0D) {
            System.out.println("Invalid query area: Negative!");
            return res;
        } else {
            double width = Math.sqrt(queryArea);

            for(int i = 0; i < queryNum; ++i) {
                double[] leftCorner = new double[2];
                double[] rightCorner = new double[2];
                double x = random.nextDouble() * (1.0D - width) + width / 2.0D;
                double y = random.nextDouble() * (1.0D - width) + width / 2.0D;
                leftCorner[0] = x - width / 2.0D;
                leftCorner[1] = y - width / 2.0D;
                rightCorner[0] = x + width / 2.0D;
                rightCorner[1] = y + width / 2.0D;
                res.add(new DoublePointRectangle(leftCorner, rightCorner));
            }

            return res;
        }
    }

    public static void main(String[] args) throws IOException {
        if (args.length != 1) {
            System.out.println("Need an input file");
        } else {
            int queryNumber = 0;
            double[] left = new double[2];
            double[] right = new double[2];
            String osmFile = null;

            int constructSeed;
            int insertionSeed;
            Integer primaryIndexSize;
            String outputDir;
            String outputTreeDir;
            Boolean bulk;
            String bulkLoadingFileName;
            int numBatches;
            int constructIterNum;
            boolean performRangeQuery;
            String compareMethod;
            double bulkLoadingInitialCap;
            Boolean isUniform;
            Boolean useNewAlgo;
            Boolean isMixedRangeQueries;
            Double queryArea;
            int querySeed;
            boolean practicalRemedyLinear;
            boolean practicalRemedyRandom;
            boolean spreadSplit;
            boolean use_osm;
            try {
                FileInputStream input = new FileInputStream(args[0]);

                try {
                    Properties prop = new Properties();
                    prop.load(input);
                    primaryIndexSize = Integer.parseInt(prop.getProperty("Primary_index_size"));
                    outputDir = prop.getProperty("Output_dir");
                    outputTreeDir = prop.getProperty("Output_tree_dir");
                    bulk = Boolean.parseBoolean(prop.getProperty("Bulk_loading"));
                    bulkLoadingFileName = prop.getProperty("Bulk_loading_file");
                    numBatches = Integer.parseInt(prop.getProperty("Number_of_batches"));
                    performRangeQuery = Boolean.parseBoolean(prop.getProperty("Perform_query"));
                    compareMethod = prop.getProperty("compare");
                    isMixedRangeQueries = Boolean.parseBoolean(prop.getProperty("Mixed_query"));
                    if (performRangeQuery) {
                        queryNumber = Integer.parseInt(prop.getProperty("Query_number"));
                    }

                    constructIterNum = Integer.parseInt(prop.getProperty("Construct_tree_iter", "0"));
                    constructSeed = Integer.parseInt(prop.getProperty("Construction_seed", "42"));
                    insertionSeed = Integer.parseInt(prop.getProperty("Insertion_seed", "43"));
                    batchSize = Integer.parseInt(prop.getProperty("Batch_size", "10000"));
                    minMaxFactor = Double.parseDouble(prop.getProperty("Min_max_factor", "0.5"));
                    dimension = Integer.parseInt(prop.getProperty("Dimension", "2"));
                    blockSize = Integer.parseInt(prop.getProperty("Block_size", "15371"));
                    bufferSize = Integer.parseInt(prop.getProperty("Buffer_size", "20"));
                    soundRemedy = Boolean.parseBoolean(prop.getProperty("SoundRemedy", "false"));
                    remedySeed = Integer.parseInt(prop.getProperty("RemedySeed", "1"));
                    bulkLoadingInitialCap = Double.parseDouble(prop.getProperty("initial_capacity", "1"));
                    isUniform = Boolean.parseBoolean(prop.getProperty("Uniform_distribution", "true"));
                    useNewAlgo = Boolean.parseBoolean(prop.getProperty("Use_new_algo", "false"));
                    queryArea = Double.parseDouble(prop.getProperty("Query_area", "0.0001"));
                    querySeed = Integer.parseInt(prop.getProperty("Query_seed", "1"));
                    practicalRemedyLinear = Boolean.parseBoolean(prop.getProperty("Practical_remedy_linear", "false"));
                    practicalRemedyRandom = Boolean.parseBoolean(prop.getProperty("Practical_remedy_random", "false"));
                    spreadSplit = Boolean.parseBoolean(prop.getProperty("spread_split", "false"));
                    left[0] = Double.parseDouble(prop.getProperty("left1", "0.0"));
                    left[1] = Double.parseDouble(prop.getProperty("left2", "0.0"));
                    right[0] = Double.parseDouble(prop.getProperty("right1", "1.0"));
                    right[1] = Double.parseDouble(prop.getProperty("right2", "1.0"));
                    use_osm = Boolean.parseBoolean(prop.getProperty("osm", "false"));
                    if (use_osm) {
                        osmFile = prop.getProperty("osm_file");
                    }
                } catch (Throwable var94) {
                    try {
                        input.close();
                    } catch (Throwable var93) {
                        var94.addSuppressed(var93);
                    }

                    throw var94;
                }

                input.close();
            } catch (IOException var95) {
                var95.printStackTrace();
                System.out.println("Cannot read the properties file!");
                return;
            }

            Comparator sFCComparator;
            if (compareMethod.equalsIgnoreCase("x")) {
                sFCComparator = COMPARE;
            } else if (compareMethod.equalsIgnoreCase("hilbert")) {
                sFCComparator = COMPARE_HILBERT;
            } else {
                if (!compareMethod.equalsIgnoreCase("peano")) {
                    System.out.println("Incorrect compare method");
                    return;
                }

                sFCComparator = COMPARE_PEANO;
            }

            universe = new DoublePointRectangle(left, right);
            HilbertRTree tree = new HilbertRTree(blockSize, universe, minMaxFactor);
            tree.outDir = outputDir;
            tree.regularElectiveSplit = spreadSplit;
            System.out.println("minMaxFactor:" + minMaxFactor + ", SoundRemedy: " + Boolean.toString(soundRemedy) + ", perform querying: " + Boolean.toString(performRangeQuery));
            System.out.println("Output dir is:" + outputDir + ", initial index size: " + primaryIndexSize + ", use compare: " + compareMethod + ", spread split: " + Boolean.toString(spreadSplit));
            System.out.println("Bulk loading capacity is: " + bulkLoadingInitialCap + ", practical remedy (linear): " + Boolean.toString(practicalRemedyLinear) + ", practical remedy (random): " + Boolean.toString(practicalRemedyRandom));
            System.out.printf("universe is: (%f, %f) (%f, %f)\n", left[0], left[1], right[0], right[1]);
            if (useNewAlgo) {
                System.out.println("Use the new splitting algorithm");
                tree.unequalRandomSplit = true;
            }

            Container fileContainer = null;
            BufferedContainer treeContainer = null;
            System.out.println("Init new HilbertRTree");
            String treePath = outputTreeDir + "/Hilbert_" + primaryIndexSize;
            fileContainer = new CounterContainer(new BlockFileContainer(treePath, blockSize));
            treeContainer = new BufferedContainer(new ConverterContainer(fileContainer, tree.nodeConverter()), new LRUBuffer(bufferSize), false);
            tree.initialize(getHilbertValue, getEntryMBR, treeContainer, keyConverter, dataConverter, createORSeparator, createORKeyRange);
            long start;
            long begin;
            if (!bulk) {
                System.out.println("Construct the tree by insertion");
                Random random = new Random((long)constructSeed);
//                int splitCount = false;
                double nanoToSeconds = 1.0E9D;
                List<double[]> points = new ArrayList();
                double[] pointQuery = new double[dimension];
                PrintWriter fw = new PrintWriter(outputDir + "/total_splits_in_construction.txt");
                double[] point = new double[dimension];
                start = System.nanoTime();
                long end = System.nanoTime();
                double duration = 0.0D;

                for(int j = 0; j < constructIterNum; ++j) {
//                    int splitCount = 0;
                    duration = 0.0D;

                    for(int k = 0; k < batchSize; ++k) {
                        point = new double[dimension];

                        for(int i = 0; i < dimension; ++i) {
                            point[i] = random.nextDouble();
                        }

                        start = System.nanoTime();
                        tree.insert(new DoublePoint(point));
                        end = System.nanoTime();
                        duration += (double)(end - start);
//                        splitCount += tree.splitCount;
                        if (k % 100 == 0) {
                            points.add(point);
                        }
                    }

//                    fw.printf("%d\n", splitCount);
                    fw.flush();
                }
            } else {
                System.out.println("Construct the tree by bulk loading");
                begin = System.nanoTime();
                int dataSize = dimension * 8;
                Converter<DoublePoint> objectConverter = new ConvertableConverter(DoublePoint.factoryFunctionDoublePoint(dimension));
                Cursor cursor = new FileInputCursor(objectConverter, new File(bulkLoadingFileName));
                int memorySizeForRuns = 10485760;
                Iterator<DoublePoint> sortedPointIterator = new MergeSorter(cursor, sFCComparator, dataSize, memorySizeForRuns, memorySizeForRuns);

                if (soundRemedy) {
                    new HilbertRTreeBulkLoading(tree, sortedPointIterator, remedySeed, primaryIndexSize);
                } else if (practicalRemedyLinear) {
                    new HilbertRTreeBulkLoading(tree, sortedPointIterator, bulkLoadingInitialCap, remedySeed, true);
                } else if (practicalRemedyRandom) {
                    new HilbertRTreeBulkLoading(tree, sortedPointIterator, bulkLoadingInitialCap, remedySeed, false);
                } else {
                    new HilbertRTreeBulkLoading(tree, sortedPointIterator, bulkLoadingInitialCap);
                }

                start = System.nanoTime();
                System.out.printf("Done bulk loading in %f min.\n", (double)(start - begin) / 6.0E10D);
                System.out.printf("Leaf node utilization: %f/%f=%f\n", tree.numOfDataInLeafNodes, tree.totalLeafNodeCapacity, tree.numOfDataInLeafNodes / tree.totalLeafNodeCapacity);
//                System.out.printf("Leaf node num = %d, leaf capacity = %d\n", tree.leafNodeNum, tree.singleNodeCapacity);
            }

//            if (isMixedRangeQueries) {
//                System.out.println("Generate range queries with mixed ranges");
//            } else {
//                System.out.println("Generate range queries with one fixed area:" + queryArea);
//            }
//
//            begin = System.nanoTime();
//            Random queryRandom = new Random((long)querySeed);
//            List<DoublePointRectangle> queries = new ArrayList();
//            List warmupQueries;
//            if (isMixedRangeQueries) {
//                warmupQueries = generateFixedRangeQueries(1.0E-4D, queryNumber, queryRandom);
//                List<DoublePointRectangle> q2 = generateFixedRangeQueries(0.001D, queryNumber, queryRandom);
//                List<DoublePointRectangle> q3 = generateFixedRangeQueries(0.01D, queryNumber, queryRandom);
//                List<DoublePointRectangle> q4 = generateFixedRangeQueries(0.05D, queryNumber, queryRandom);
//                List<DoublePointRectangle> q5 = generateFixedRangeQueries(0.2D, queryNumber, queryRandom);
//                ((List)queries).addAll(warmupQueries);
//                ((List)queries).addAll(q2);
//                ((List)queries).addAll(q3);
//                ((List)queries).addAll(q4);
//                ((List)queries).addAll(q5);
//            } else {
//                queries = generateFixedRangeQueries(queryArea, queryNumber, queryRandom);
//            }
//
//            warmupQueries = generateFixedRangeQueries(1.0E-4D, 5, queryRandom);
//            long done = System.nanoTime();
//            System.out.printf("Done query generating in %f min.\n", (double)(done - begin) / 6.0E10D);
//            System.out.println("Begin inserting data and counting splits");
//            if (!isUniform) {
//                System.out.println("Using normal distribtuion");
//            }

            begin = System.nanoTime();
            Random random = new Random((long)insertionSeed);
            List<double[]> osmPoints = new ArrayList();
            if (use_osm) {
                System.out.println("Read in OpenStreetMap data");
                File myFile = new File(osmFile);
                Scanner myReader = new Scanner(myFile);
                int counter = 0;
                while(myReader.hasNextLine()) {
                    String line = myReader.nextLine();
                    String[] arrOfStr = line.split(" ", 2);
                    double[] point = new double[]{Double.parseDouble(arrOfStr[0]), Double.parseDouble(arrOfStr[1])};
                    if (!(point[0] < left[0]) && !(point[0] > right[0]) && !(point[0] < left[1]) && !(point[1] > right[1])) {
                        osmPoints.add(point);
                    } else {
                        System.out.println(point[0] + "," + point[1]);
                    }
                    counter++;
                    if (counter >= batchSize * numBatches) {
                        break;
                    }
                }

                myReader.close();
            }

            ArrayList<Integer> totalSplits = new ArrayList();
            ArrayList<Integer> leafSplits = new ArrayList();
            ArrayList<Integer> nonLeafSplits = new ArrayList();
            ArrayList<Integer> indexSizeList = new ArrayList();
            ArrayList<Integer> leafNodeAccess = new ArrayList();
            ArrayList<Integer> nonLeafNodeAccess = new ArrayList();
            ArrayList<Integer> treeHeight = new ArrayList();
            ArrayList<Integer> leafNumber = new ArrayList();
            ArrayList<Integer> queryResultSizeList = new ArrayList();
            ArrayList<ArrayList<Integer>> queryLeafAccessList = new ArrayList();
            ArrayList<ArrayList<Integer>> queryNonLeafAccessList = new ArrayList();
            ArrayList<ArrayList<Double>> queryTimeList = new ArrayList();
            ArrayList<Double> insertTimeList = new ArrayList();
            ArrayList<Double> avgLeafUtil = new ArrayList();
            ArrayList<Double> overlapArea = new ArrayList();
            ArrayList<Double> areaDifference = new ArrayList();
            double nanoToSeconds = 1.0E9D;
            int indexSize = primaryIndexSize;
            boolean fstTimeQuery = true;
            int osmIndex = 0;

            for(int j = 0; j < numBatches && (!use_osm || osmIndex < osmPoints.size()); ++j) {
                ++tree.curBatchNum;
                System.out.println(j);
                int splitCount = 0;
                int leafSplitCount = 0;
                int nonLeafSplitCount = 0;
                int leafAccess = 0;
                int nonLeafAccess = 0;
                double insertionTime = 0.0D;
                double avgLeafUtilRatio = 0.0D;
                int rootLevel = 0;
                int leafNodeNum = 0;
                double overlapValue = 0.0D;
                double areaDiff = 0.0D;

//                long start;
                long end;
                int q;
                for(q = 0; q < batchSize; ++q) {
                    double[] point = new double[dimension];
                    if (!use_osm) {
                        for(int i = 0; i < dimension; ++i) {
                            if (isUniform) {
                                point[i] = random.nextDouble();
                            } else {
                                point[i] = random.nextGaussian() * 0.5D / 3.0D + 0.5D;
                            }
                        }
                    } else {
                        point = (double[])osmPoints.get(osmIndex);
                        ++osmIndex;
                        if (osmIndex >= osmPoints.size()) {
                            break;
                        }
                    }

                    start = System.nanoTime();
                    tree.insert(new DoublePoint(point));
                    end = System.nanoTime();
                    insertionTime += (double)(end - start);
//                    splitCount += tree.splitCount;
                    leafSplitCount += tree.leafSplitCount;
                    System.out.println("Split is " + tree.leafSplitCount);
                    nonLeafSplitCount += tree.nonLeafSplitCount;
                    overlapValue += tree.overlapArea;
                    areaDiff += tree.areaDifference;
//                    leafAccess += tree.leafAccess;
//                    nonLeafAccess += tree.nonLeafAccess;
                    avgLeafUtilRatio += tree.numOfDataInLeafNodes / tree.totalLeafNodeCapacity;
                    rootLevel += tree.rootEntry().level();
                    leafNodeNum += tree.leafNodeNum;
                }

                rootLevel /= batchSize;
                leafNodeNum /= batchSize;
                avgLeafUtilRatio /= (double)batchSize;
                indexSize += batchSize;
                indexSizeList.add(indexSize);
                totalSplits.add(splitCount);
                leafSplits.add(leafSplitCount);
                nonLeafSplits.add(nonLeafSplitCount);
                leafNodeAccess.add(leafAccess);
                nonLeafNodeAccess.add(nonLeafAccess);
                insertTimeList.add(insertionTime / nanoToSeconds);
                avgLeafUtil.add(avgLeafUtilRatio);
                overlapArea.add(overlapValue);
                areaDifference.add(areaDiff);
                treeHeight.add(rootLevel);
                leafNumber.add(leafNodeNum);
//                if (performRangeQuery) {
//                    if (fstTimeQuery) {
//                        q = 0;
//
//                        while(true) {
//                            if (q >= warmupQueries.size()) {
//                                fstTimeQuery = false;
//                                break;
//                            }
//
//                            Cursor results = tree.queryOR((Rectangle)warmupQueries.get(q), 0);
//
//                            while(results.hasNext()) {
//                                results.next();
//                            }
//
//                            results.close();
//                            ++q;
//                        }
//                    }

//                    ArrayList<Double> queryTime = new ArrayList();
//                    ArrayList<Integer> queryLeafAccess = new ArrayList();
//                    ArrayList<Integer> queryNonLeafAccess = new ArrayList();
//
//                    for(int q = 0; q < ((List)queries).size(); ++q) {
//                        int[] nodeAccess = new int[]{0, 0};
//                        start = System.nanoTime();
//                        Cursor results = tree.queryOR((Rectangle)((List)queries).get(q), 0, nodeAccess, outputDir);
//
//                        while(results.hasNext()) {
//                            results.next();
//                        }
//
//                        end = System.nanoTime();
//                        double qt = (double)(end - start);
//                        queryTime.add(qt / 1000.0D);
//                        queryLeafAccess.add(nodeAccess[0]);
//                        queryNonLeafAccess.add(nodeAccess[1]);
//                        results.close();
//                    }
//
//                    queryTimeList.add(queryTime);
//                    queryLeafAccessList.add(queryLeafAccess);
//                    queryNonLeafAccessList.add(queryNonLeafAccess);
//                }
            }

            writeIntToFile(outputDir + "/totalSplits.txt", totalSplits);
            writeIntToFile(outputDir + "/leafSplits.txt", leafSplits);
            writeIntToFile(outputDir + "/nonLeafSplits.txt", nonLeafSplits);
            writeIntToFile(outputDir + "/indexSize.txt", indexSizeList);
            writeIntToFile(outputDir + "/leafNodeAccess.txt", leafNodeAccess);
            writeIntToFile(outputDir + "/nonLeafNodeAccess.txt", nonLeafNodeAccess);
            writeIntToFile(outputDir + "/treeHeight.txt", treeHeight);
            writeIntToFile(outputDir + "/leafNodeNum.txt", leafNumber);
            writeIntToFile(outputDir + "/queryResultSize.txt", queryResultSizeList);
            writeIntToFile(outputDir + "/queryLeafNodeAccess.txt", queryLeafAccessList);
            writeIntToFile(outputDir + "/queryNonLeafNodeAccess.txt", queryNonLeafAccessList);
            writeDoubleToFile(outputDir + "/insertTime.txt", insertTimeList);
            writeDoubleToFile(outputDir + "/avgLeafRatio.txt", avgLeafUtil);
            writeDoubleToFile(outputDir + "/queryTime.txt", queryTimeList);
            writeDoubleToFile(outputDir + "/overlappedArea.txt", overlapArea);
            writeDoubleToFile(outputDir + "/areaDiff.txt", areaDifference);
//            writeDoubleToFile(outputDir + "/splitOverlapRatio.txt", tree.splitOverlapRatioList);
//            writeDoubleToFile(outputDir + "/splitOverlapRatio2.txt", tree.splitOverlapRatioList2);
//            writeDoubleToFile(outputDir + "/splitRatio.txt", tree.splitRatioList);
//            writeIntToFile(outputDir + "/nonLeafSplitLevel.txt", tree.nonLeafSplitLevel);
            long done = System.nanoTime();
            System.out.printf("Done in %f min.\n", (double)(done - begin) / 6.0E10D);
        }
    }

    static {
        keyConverter = new MeasuredFixedSizeConverter(LongConverter.DEFAULT_INSTANCE);
        dataConverter = new MeasuredConverter() {
            public int getMaxObjectSize() {
                return SimpleHilbertRTreeTest.dimension * 8;
            }

            public Object read(DataInput dataInput, Object object) throws IOException {
                return ConvertableConverter.DEFAULT_INSTANCE.read(dataInput, new DoublePoint(SimpleHilbertRTreeTest.dimension));
            }

            public void write(DataOutput dataOutput, Object object) throws IOException {
                ConvertableConverter.DEFAULT_INSTANCE.write(dataOutput, (DoublePoint)object);
            }
        };
        getEntryMBR = new AbstractFunction() {
            public Object invoke(Object o) {
                DoublePoint p = (DoublePoint)o;
                return new DoublePointRectangle(p, p);
            }
        };
        getHilbertValue = new AbstractFunction() {
            public Object invoke(Object point) {
                if (point instanceof DoublePoint) {
                    DoublePoint middlePoint = (DoublePoint)point;
                    double x = middlePoint.getValue(0);
                    double y = middlePoint.getValue(1);
                    return new Long(SpaceFillingCurves.hilbert2d((int)(x * 1048576.0D), (int)(y * 1048576.0D)));
                } else {
                    throw new IllegalArgumentException();
                }
            }
        };
        createORSeparator = new AbstractFunction() {
            public Object invoke(Object key, Object mbr) {
                return new SimpleHilbertRTreeTest.DoublePointRectangleSep((Long)key, (DoublePointRectangle)mbr);
            }
        };
        createORKeyRange = new AbstractFunction() {
            public Object invoke(List arguments) {
                if (arguments.size() != 3) {
                    throw new IllegalArgumentException();
                } else {
                    Long min = (Long)arguments.get(0);
                    Long max = (Long)arguments.get(1);
                    DoublePointRectangle entryMBR = (DoublePointRectangle)arguments.get(2);
                    return new SimpleHilbertRTreeTest.LongRange(min, max, entryMBR);
                }
            }
        };
        COMPARE_HILBERT = new Comparator() {
            protected double[] uni;
            protected double[] uniDeltas;

            public int compare(Object o1, Object o2) {
                if (SimpleHilbertRTreeTest.dimension != 2) {
                    throw new IllegalArgumentException();
                } else {
                    if (this.uni == null) {
                        this.uni = (double[])SimpleHilbertRTreeTest.universe.getCorner(false).getPoint();
                        this.uniDeltas = SimpleHilbertRTreeTest.universe.deltas();
                    }

                    double[] leftBorders1 = (double[])((DoublePoint)o1).getPoint();
                    double[] leftBorders2 = (double[])((DoublePoint)o2).getPoint();
                    double x1 = (leftBorders1[0] - this.uni[0]) / this.uniDeltas[0];
                    double y1 = (leftBorders1[1] - this.uni[1]) / this.uniDeltas[1];
                    double x2 = (leftBorders2[0] - this.uni[0]) / this.uniDeltas[0];
                    double y2 = (leftBorders2[1] - this.uni[1]) / this.uniDeltas[1];
                    long h1 = SpaceFillingCurves.hilbert2d((int)(x1 * 1048576.0D), (int)(y1 * 1048576.0D));
                    long h2 = SpaceFillingCurves.hilbert2d((int)(x2 * 1048576.0D), (int)(y2 * 1048576.0D));
//                    System.out.printf("(%f,%f) and (%f,%f) => %d and %d\n", x1, y1,x2,y2, h1, h2);
                    return h1 < h2 ? -1 : (h1 == h2 ? 0 : 1);
                }
            }
        };
        COMPARE_PEANO = new Comparator() {
            protected double[] uni;
            protected double[] uniDeltas;

            public int compare(Object o1, Object o2) {
                if (SimpleHilbertRTreeTest.dimension != 2) {
                    throw new IllegalArgumentException();
                } else {
                    if (this.uni == null) {
                        this.uni = (double[])SimpleHilbertRTreeTest.universe.getCorner(false).getPoint();
                        this.uniDeltas = SimpleHilbertRTreeTest.universe.deltas();
                    }

                    double[] leftBorders1 = (double[])((DoublePoint)o1).getPoint();
                    double[] leftBorders2 = (double[])((DoublePoint)o2).getPoint();
                    double x1 = (leftBorders1[0] - this.uni[0]) / this.uniDeltas[0];
                    double y1 = (leftBorders1[1] - this.uni[1]) / this.uniDeltas[1];
                    double x2 = (leftBorders2[0] - this.uni[0]) / this.uniDeltas[0];
                    double y2 = (leftBorders2[1] - this.uni[1]) / this.uniDeltas[1];
                    long h1 = SpaceFillingCurves.peano2d((int)(x1 * 1048576.0D), (int)(y1 * 1048576.0D));
                    long h2 = SpaceFillingCurves.peano2d((int)(x2 * 1048576.0D), (int)(y2 * 1048576.0D));
                    return h1 < h2 ? -1 : (h1 == h2 ? 0 : 1);
                }
            }
        };
        COMPARE = new Comparator() {
            public int compare(Object o1, Object o2) {
                double[] p1 = (double[])((DoublePoint)o1).getPoint();
                double[] p2 = (double[])((DoublePoint)o2).getPoint();
                return p1[0] < p2[0] ? -1 : (p1[0] == p2[0] ? 0 : 1);
            }
        };
    }

    public static class LongRange extends ORKeyRange {
        public LongRange(Long min, Long max, DoublePointRectangle entryMBR) {
            super(min, max, entryMBR);
        }

        public Object clone() {
            return new SimpleHilbertRTreeTest.LongRange((Long)this.sepValue, (Long)this.maxBound, (DoublePointRectangle)this.entryMBR.clone());
        }
    }

    public static class DoublePointRectangleSep extends ORSeparator {
        public DoublePointRectangleSep(Long separatorValue, DoublePointRectangle entryMBR) {
            super(separatorValue, entryMBR);
        }

        public Object clone() {
            return new SimpleHilbertRTreeTest.DoublePointRectangleSep((Long)this.sepValue(), (DoublePointRectangle)this.entryMBR.clone());
        }
    }
}
