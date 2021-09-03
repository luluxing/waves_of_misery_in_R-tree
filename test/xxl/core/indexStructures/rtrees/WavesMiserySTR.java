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
import xxl.core.cursors.Cursor;
import xxl.core.cursors.Cursors;
import xxl.core.cursors.sources.io.FileInputCursor;
import xxl.core.functions.AbstractFunction;
import xxl.core.functions.Function;
import xxl.core.functions.Functional.UnaryFunction;
import xxl.core.indexStructures.*;
import xxl.core.indexStructures.Tree.Query.Candidate;
import xxl.core.io.LRUBuffer;
import xxl.core.io.converters.ConvertableConverter;
import xxl.core.io.converters.Converter;
import xxl.core.io.converters.IntegerConverter;
import xxl.core.io.converters.LongConverter;
import xxl.core.spatial.SpaceFillingCurves;
import xxl.core.spatial.points.DoublePoint;
import xxl.core.spatial.rectangles.DoublePointRectangle;
import xxl.core.spatial.rectangles.Rectangle;
import xxl.core.spatial.rectangles.Rectangles;

import java.io.*;
import java.util.*;

/**
 * Simple Test of the HilbertRTree with random DoublePoint ({@link DoublePoint})
 */
public class WavesMiserySTR {

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
     */
    static public Function getHilbertValue = new AbstractFunction(){
        public Object invoke(Object point){
            if (point instanceof DoublePoint){
                DoublePoint middlePoint = (DoublePoint)point;
                double x = middlePoint.getValue(0);
                double y = middlePoint.getValue(1);
                return new Long(SpaceFillingCurves.hilbert2d((int)
                        (x*WavesMisery.FILLING_CURVE_PRECISION),(int) (y*WavesMisery.FILLING_CURVE_PRECISION)));
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

    /**
     * Starts simpleHilbertRTree Test
     * @param args path To HilbertRTree
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

        if (args.length != 1) {
            System.out.println("Need an input file");
            return;
        }

        PropertyReader reader = new PropertyReader(args[0]);
        reader.print_properties();


        PropertyReader.wsm_print("Init new tree");
        String treePath = reader.m_outputDir + "/" + reader.m_treeName + reader.m_primaryIndexSize;
        Container fileContainer = new BlockFileContainer(treePath, reader.m_blockSize);
        Container ioContainer = new CounterContainer(fileContainer);

        // size of a data entry = size of DoublePoint (1 double per dimension)
        int dataSize = reader.m_dimension * 8;
        // size of a descriptor = size of DoublePointRectangle (2 doubles per dimension)
        int descriptorSize = reader.m_dimension * 2 * 8;
        Tree tree = null;
        CounterContainer lowerCounterContainer;
        Container treeContainer;
        FilteredCounterContainer container;
        FilteredCounterContainer lowerLeafCounter, upperLeafCounter;

        if (reader.m_treeName.equalsIgnoreCase("rstar"))
            tree = new RTree(reader.m_reinsertion);
        else if (reader.m_treeName.equalsIgnoreCase("linear"))
            tree = new LinearRTree();
        else if (reader.m_treeName.equalsIgnoreCase("quadratic"))
            tree = new QuadraticRTree();
        else if (reader.m_treeName.equalsIgnoreCase("rrstar"))
            tree = new RevisedRTree(WavesMisery.dimension);
        else {
            PropertyReader.wsm_print("Tree name not found. Accept only: rstar, linear, quadratic, rrstar.");
            return;
        }

        Converter converter = ((RTree)tree).nodeConverter(new ConvertableConverter(
                                                    DoublePoint.factoryFunctionDoublePoint(reader.m_dimension)),
                                                    reader.m_dimension);
        ConverterContainer converterContainer = new ConverterContainer(fileContainer, converter);
        treeContainer = converterContainer;
        if (reader.m_bufferSize != 0) {
            // buffered container that counts the access to the buffered RTree
            treeContainer = new BufferedContainer(converterContainer,
                                                    new LRUBuffer(reader.m_bufferSize), true);
        }

        // container that counts the access to the RTree
        container = new FilteredCounterContainer(treeContainer,
                                                    FilteredCounterContainer.TREE_LEAF_NODE_COUNTER_FUNCTION);

        ((RTree)tree).initialize(WavesMisery.getEntryMBR, container, reader.m_blockSize, dataSize,
                descriptorSize, reader.m_minMaxFactor);

        tree.outDir = reader.m_outputDir;
        tree.unequalRandomSplit = reader.m_unequalRandSplit;
        tree.regularElectiveSplit = reader.m_regularSplit;
        if (tree.regularElectiveSplit) {
            tree.splitFrequency = reader.m_splitFreq;
//            tree.minElectiveSplit = reader.m_minCap;
            tree.fastElectiveSplit = reader.m_fastVersion;
            tree.overflowPageNum = reader.m_overflowPageNum;
        }


        /** Bulk-loading */
        PropertyReader.wsm_print("Construct the tree by bulk loading");
        long begin = System.nanoTime();

        Converter<DoublePoint> objectConverter =
                new ConvertableConverter<DoublePoint>(DoublePoint.factoryFunctionDoublePoint(reader.m_dimension));

        int[] sortingFunction = {0, 1};
        STRBulkLoader<DoublePoint> strBulkloader = new STRBulkLoader<DoublePoint>((RTree)tree,
                                                                    treePath+"str",
                                                                        reader.m_dimension,
                                                                        reader.m_blockSize,
                                                                        reader.m_minMaxFactor,
                                                                        reader.m_bulkLoadingInitialCap,
                                                                        sortingFunction);

        int objectNumber = Cursors.count(new FileInputCursor<DoublePoint>(objectConverter, new File(reader.m_bulkLoadingFileName)));

        // with this value we provide an available memory for initial run generation;
        int memorySizeForRuns = 1024 * 1024 * 10;

        strBulkloader.init(objectNumber, memorySizeForRuns, dataSize, objectConverter,
                new UnaryFunction<DoublePoint, DoublePointRectangle>() {
                            @Override
                            public DoublePointRectangle invoke(DoublePoint arg) {
                                return new DoublePointRectangle(arg, arg);
                            }
        });
        // conduct bulk-loading
        strBulkloader.buildRTree(new FileInputCursor<>(objectConverter, new File(reader.m_bulkLoadingFileName)));

        tree = strBulkloader.getRTree();

        long done = System.nanoTime();

        PropertyReader.wsm_print("Bulk-loading finishes in ", (done - begin) / (1000.0*1000*1000*60), " min.");
        PropertyReader.wsm_print("Total data in leaf nodes is ", tree.numOfDataInLeafNodes);
        PropertyReader.wsm_print("Total leaf node capacity is ", tree.totalLeafNodeCapacity);
        PropertyReader.wsm_print("Initial leaf node utilization is ", tree.numOfDataInLeafNodes/tree.totalLeafNodeCapacity);
        PropertyReader.wsm_print("Initial leaf node capacity is ", tree.singleLeafNodeCapacity);
        PropertyReader.wsm_print("Total leaf node number is ", tree.leafNodeNum);


        container.reset();

        if (reader.m_bufferSize > 0) {
//            System.out.println("Flushing buffers");
            container.flush();
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

        Random random = new Random(reader.m_insertionSeed);

        ArrayList<Integer> leafSplitsList = new ArrayList<>();
        ArrayList<Integer> nonLeafSplitsList = new ArrayList<>();
        ArrayList<Integer> insertionLeafAccessList = new ArrayList<>();
        ArrayList<Integer> insertionNonLeafAccessList = new ArrayList<>();
        ArrayList<Integer> leafSplitsList2 = new ArrayList<>();
        ArrayList<Integer> nonLeafSplitsList2 = new ArrayList<>();
        ArrayList<Integer> insertionLeafAccessList2 = new ArrayList<>();
        ArrayList<Integer> insertionNonLeafAccessList2 = new ArrayList<>();

        ArrayList<Integer> treeHeightList = new ArrayList<>();
        ArrayList<Integer> leafNumberList = new ArrayList<>();
        ArrayList<Integer> queryResultSizeList = new ArrayList<>();

        ArrayList<ArrayList<Integer>> queryLeafAccessList = new ArrayList<>();
        ArrayList<ArrayList<Integer>> queryNonLeafAccessList = new ArrayList<>();
        ArrayList<ArrayList<Double>> queryTimeList = new ArrayList<>();

        ArrayList<Double> insertTimeList = new ArrayList<>();

        ArrayList<Double> avgLeafUtilList = new ArrayList<>();
        ArrayList<Double> avgLeafDataOverLeafNumList = new ArrayList<>();
        ArrayList<Double> avgTotalDataOverLeafNumList = new ArrayList<>();
        ArrayList<Double> avgTotalDataOverAllPageNumList = new ArrayList<>();

        ArrayList<Double> overlapAreaList = new ArrayList<>();
        ArrayList<Double> areaDifferenceList = new ArrayList<>();

        ArrayList<Double> avgIOPerQuery = new ArrayList<>();
        ArrayList<Double> avgLeafPerQuery = new ArrayList<>();

        ArrayList<Double> splitOverlapList = new ArrayList<>();
        ArrayList<Double> areaDiffList = new ArrayList<>();
        ArrayList<Double> splitOverlapRatioList = new ArrayList<>();
        ArrayList<Double> splitOverlapRatio2List = new ArrayList<>();
        ArrayList<Double> splitRatioList = new ArrayList<>();

        double nanoToSeconds = 1000 * 1000 * 1000;

        double[] point;
        long start;
        long end;

        int leafSplitCount, nonLeafSplitCount;
        int insertionLeafAccess, insertionNonLeafAccess;
        double insertionTime = 0;
        double avgLeafUtilRatio, avgLeafDataOverLeafNum, avgTotalDataOverLeafNum, avgTotalDataOverAllPageNum;
        int rootLevel;
        int leafNodeNum;
        double overlapValue, areaDiff;
        double splitOverlap, splitOverlapRatio, splitOverlapRatio2, splitRatio;

        int numBatches = reader.m_numBatches;
        int batchSize = reader.m_batchSize;
        int dimension = reader.m_dimension;
        boolean isUniform = reader.m_isUniform;
        String line;
        boolean end_early = false, first_query = true;

        for (int j = 0; j < numBatches; j++) {
            if (end_early)
                break;

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
            leafNodeNum = 0;
            overlapValue = 0.0;
            areaDiff = 0.0;
            splitOverlap  = 0.0;
            areaDiff = 0.0;
            splitOverlapRatio = 0.0;
            splitOverlapRatio2 = 0.0;
            splitRatio = 0.0;
//            upperLeafCounter.reset();
            container.reset();
            container.flush();
            for (int k = 0; k < batchSize; k++) {
                // Generate a new point
                point = new double[dimension];
                if (reader.m_realData) {
                    if ((line = br.readLine()) != null) {
                        String[] values = line.split(",");
                        point[0] = Double.parseDouble(values[0]);
                        point[1] = Double.parseDouble(values[1]);
//                        System.out.println(point[0] + ", " + point[1]);
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
                overlapValue += tree.overlapArea;
//                areaDiff += tree.areaDifference;

                insertionLeafAccess += tree.insertionLeafAccess;
                insertionNonLeafAccess += tree.insertionNonLeafAccess;

                rootLevel += tree.rootEntry().level();
                leafNodeNum += tree.leafNodeNum;

                avgLeafUtilRatio += tree.numOfDataInLeafNodes / tree.totalLeafNodeCapacity;
                avgLeafDataOverLeafNum += (tree.numOfDataInLeafNodes - tree.totalOverflowEntryNum)
                        / (1.0 * tree.singleLeafNodeCapacity * tree.leafNodeNum);
                avgTotalDataOverLeafNum += tree.numOfDataInLeafNodes / (1.0 * tree.singleLeafNodeCapacity * tree.leafNodeNum);
                avgTotalDataOverAllPageNum += 1.0 * tree.numOfDataInLeafNodes /
                        (tree.singleLeafNodeCapacity *(tree.leafNodeNum + tree.totalOverflowPageNum));

                splitOverlap += tree.splitOverlap;
                areaDiff += tree.areaDiff;
                splitOverlapRatio += tree.splitOverlapRatio;
                splitOverlapRatio2 += tree.splitOverlapRatio2;
                splitRatio += tree.splitRatio;
            }

            rootLevel /= batchSize;
            treeHeightList.add(rootLevel);

            leafNodeNum /= batchSize;
            leafNumberList.add(leafNodeNum);

            avgLeafUtilRatio /= batchSize;
            avgLeafUtilList.add(avgLeafUtilRatio);

            avgLeafDataOverLeafNum /= batchSize;
            avgLeafDataOverLeafNumList.add(avgLeafDataOverLeafNum);

            avgTotalDataOverLeafNum /= batchSize;
            avgTotalDataOverLeafNumList.add(avgTotalDataOverLeafNum);

            avgTotalDataOverAllPageNum /=batchSize;
            avgTotalDataOverAllPageNumList.add(avgTotalDataOverAllPageNum);

            leafSplitsList.add(leafSplitCount);
            nonLeafSplitsList.add(nonLeafSplitCount);
            insertionLeafAccessList.add(insertionLeafAccess);
            insertionNonLeafAccessList.add(insertionNonLeafAccess);

            leafSplitsList2.add(container.insertPredicates);
            nonLeafSplitsList2.add(container.inserts - container.insertPredicates);
            insertionLeafAccessList2.add(container.getsPredicates);
            insertionNonLeafAccessList2.add(container.gets - container.getsPredicates);
            insertTimeList.add(insertionTime / nanoToSeconds);

            overlapAreaList.add(overlapValue);
            areaDifferenceList.add(areaDiff);

            splitOverlapList.add(splitOverlap);
            areaDiffList.add(areaDiff);
            splitOverlapRatioList.add(splitOverlapRatio);
            splitOverlapRatio2List.add(splitOverlapRatio2);
            splitRatioList.add(splitRatio);

            // We need to perform range query after each bath
            if (!reader.m_performRangeQuery)
                continue;

        }

        WavesMisery.writeIntToFile(reader.m_outputDir + reader.m_treeName + "_leafSplits.txt", leafSplitsList);
        WavesMisery.writeIntToFile(reader.m_outputDir + reader.m_treeName + "_nonLeafSplits.txt", nonLeafSplitsList);
        WavesMisery.writeIntToFile(reader.m_outputDir + reader.m_treeName + "_leafSplits2.txt", leafSplitsList2);
        WavesMisery.writeIntToFile(reader.m_outputDir + reader.m_treeName + "_nonLeafSplits2.txt", nonLeafSplitsList2);
        WavesMisery.writeDoubleToFile(reader.m_outputDir + reader.m_treeName + "_insertTime.txt", insertTimeList);
        WavesMisery.writeIntToFile(reader.m_outputDir + reader.m_treeName + "_insertionLeafAccess.txt", insertionLeafAccessList);
        WavesMisery.writeIntToFile(reader.m_outputDir + reader.m_treeName + "_insertionNonLeafAccess.txt", insertionLeafAccessList);
        WavesMisery.writeIntToFile(reader.m_outputDir + reader.m_treeName + "_insertionLeafAccess2.txt", insertionLeafAccessList2);
        WavesMisery.writeIntToFile(reader.m_outputDir + reader.m_treeName + "_insertionNonLeafAccess2.txt", insertionLeafAccessList2);

        WavesMisery.writeIntToFile(reader.m_outputDir + reader.m_treeName + "_treeHeight.txt", treeHeightList);
        WavesMisery.writeIntToFile(reader.m_outputDir + reader.m_treeName + "_leafNodeNum.txt", leafNumberList);

        WavesMisery.writeIntToFile(reader.m_outputDir + reader.m_treeName + "_queryResultSize.txt", queryResultSizeList);
        WavesMisery.writeIntToFile(reader.m_outputDir + reader.m_treeName + "_queryLeafAccessList.txt", queryLeafAccessList);
        WavesMisery.writeIntToFile(reader.m_outputDir + reader.m_treeName + "_queryNonLeafAccessList.txt", queryNonLeafAccessList);
        WavesMisery.writeDoubleToFile(reader.m_outputDir + reader.m_treeName + "_queryTime.txt", queryTimeList);
        WavesMisery.writeDoubleToFile(reader.m_outputDir + reader.m_treeName + "_queryAvgIO.txt", avgIOPerQuery);
        WavesMisery.writeDoubleToFile(reader.m_outputDir + reader.m_treeName + "_queryAvgLeaf.txt", avgLeafPerQuery);

        /** Averaged leaf utilization ratio per batch */
        // There are 4 different ways to compute average leaf node utilization
        // They only differ when overflowing pages are enabled
        WavesMisery.writeDoubleToFile(reader.m_outputDir + reader.m_treeName + "_avgLeafRatio.txt", avgLeafUtilList);
        WavesMisery.writeDoubleToFile(reader.m_outputDir + reader.m_treeName + "_avgLeafDataOverLeafNumList.txt", avgLeafDataOverLeafNumList);
        WavesMisery.writeDoubleToFile(reader.m_outputDir + reader.m_treeName + "_avgTotalDataOverLeafNumList.txt", avgTotalDataOverLeafNumList);
        WavesMisery.writeDoubleToFile(reader.m_outputDir + reader.m_treeName + "_avgTotalDataOverAllPageNumList.txt", avgTotalDataOverAllPageNumList);

        WavesMisery.writeDoubleToFile(reader.m_outputDir + reader.m_treeName + "_1.txt", splitOverlapList);
        WavesMisery.writeDoubleToFile(reader.m_outputDir + reader.m_treeName + "_11.txt", areaDiffList);
        WavesMisery.writeDoubleToFile(reader.m_outputDir + reader.m_treeName + "_111.txt", splitOverlapRatioList);
        WavesMisery.writeDoubleToFile(reader.m_outputDir + reader.m_treeName + "_1111.txt", splitOverlapRatio2List);
        WavesMisery.writeDoubleToFile(reader.m_outputDir + reader.m_treeName + "_11111.txt", splitRatioList);

        done = System.nanoTime();
        PropertyReader.wsm_print("Done in ", (done - begin) / (1000.0*1000*1000*60), " min.");
    }


}
