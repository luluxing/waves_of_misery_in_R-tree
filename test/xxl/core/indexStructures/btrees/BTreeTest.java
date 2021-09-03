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

package xxl.core.indexStructures.btrees;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;

import xxl.core.collections.containers.Container;
import xxl.core.collections.containers.CounterContainer;
import xxl.core.collections.containers.io.BlockFileContainer;
import xxl.core.collections.containers.io.BufferedContainer;
import xxl.core.collections.containers.io.ConverterContainer;
import xxl.core.comparators.ComparableComparator;
import xxl.core.cursors.Cursors;
import xxl.core.cursors.sorters.MergeSorter;
import xxl.core.cursors.sources.Permutator;
import xxl.core.cursors.sources.PermutatorDouble;
import xxl.core.functions.AbstractFunction;
import xxl.core.functions.Constant;
import xxl.core.functions.Function;
import xxl.core.indexStructures.*;
import xxl.core.indexStructures.rtrees.PropertyReader;
import xxl.core.io.FilesystemOperations;
import xxl.core.io.LRUBuffer;
import xxl.core.io.LogFilesystemOperations;
import xxl.core.io.XXLFilesystem;
import xxl.core.io.converters.DoubleConverter;
import xxl.core.io.converters.IntegerConverter;
import xxl.core.util.Interval1D;

/**
 * Creates and tests a BTree with minimum (maximum) capacity given
 * by the command line argument 1 (2) by inserting the numbers of
 * a permutation of the numbers ranging from 0 (inclusive) to the value
 * of command line argument 3 (exclusive).
 * This applications thereafter performs a range-query [command line 
 * argument 6, command line argument 7].
 * <p>
 * This applications has to be called with 2 to 7 parameters:
 * <ol>
 * <li>minimum capacity of nodes</li>
 * <li>maximum capacity of nodes</li>
 * <li>number of elements to be inserted. Default: 1000</li>
 * <li>type of insertion: tuple or bulk. Default: tuple</li>
 * <li>buffersize (number of node-objects). Default: 100</li>
 * <li>minimum key for range query. Default: 0</li>
 * <li>maximum key for range query. Default: 0</li>
 * </ol>
 */
public class BTreeTest {

	
	/**
	 * A function that creates a descriptor for a given object.
	 */
	public static Function GET_DESCRIPTOR = new AbstractFunction() {
		public Object invoke (Object object) {
			return new Interval1D(object);
		}
	};

	/** 
	 * The main method performing the tests.
	 * 
	 * @param args command line parameters
	 * @throws Exception exceptions are thrown
	 */
	public static void main (String [] args) throws Exception {
//		System.out.println("BTreeTest: an example using xxl.core.indexStructures.BTree\n");
		
//		if (args.length<2) {
//			System.out.println("This applications has to be called with 2 to 7 parameters:");
//			System.out.println("1. minimum capacity of nodes");
//			System.out.println("2. maximum capacity of nodes");
//			System.out.println("3. number of elements to be inserted. Default: 1000");
//			System.out.println("4. type of insertion: tuple or bulk. Default: tuple");
//			System.out.println("5. buffersize (number of node-objects). Default: 100");
//			System.out.println("6. minimum key for range query. Default: 0");
//			System.out.println("7. maximum key for range query. Default: 0");
//			return;
//		}

		PropertyReader reader = new PropertyReader(args[0]);
		reader.print_properties();
		PropertyReader.wsm_print("Init new tree");
		String treePath = reader.m_outputDir + "/" + reader.m_treeName + reader.m_primaryIndexSize;

		BTree btree = new BTree();
		XXLFilesystem fs = new XXLFilesystem(true);
		FilesystemOperations fso = fs.getFilesystemOperations();
		fso = new LogFilesystemOperations(fso, System.out, true);
		BlockFileContainer rawContainer = new BlockFileContainer(treePath, reader.m_blockSize, fso);
		Container ioContainer = new CounterContainer(rawContainer);

		Container converterContainer = new ConverterContainer(ioContainer, btree.nodeConverter
				(DoubleConverter.DEFAULT_INSTANCE, DoubleConverter.DEFAULT_INSTANCE, new ComparableComparator()));
		Container treeContainer = new BufferedContainer(converterContainer, new LRUBuffer(reader.m_bufferSize), true);
		FilteredCounterContainer container = new FilteredCounterContainer(treeContainer, FilteredCounterContainer.TREE_LEAF_NODE_COUNTER_FUNCTION);

		int dataSize = DoubleConverter.SIZE;
		int descriptorSize = DoubleConverter.SIZE * 2; // TODO: ?
		btree.initialize(null, null, GET_DESCRIPTOR, reader.m_blockSize, container, dataSize,
				descriptorSize, reader.m_minMaxFactor);

		// Pass parameters to the tree
		btree.outDir = reader.m_outputDir;
		btree.unequalRandomSplit = reader.m_unequalRandSplit;
		btree.regularElectiveSplit = reader.m_regularSplit;
		if (btree.regularElectiveSplit) {
			btree.splitFrequency = reader.m_splitFreq;
			btree.fastElectiveSplit = reader.m_fastVersion;
			btree.overflowPageNum = reader.m_overflowPageNum;
		}
//		btree.initialize(GET_DESCRIPTOR, container, mincap, maxcap);

//		int mincap = Integer.parseInt(args[0]);
//		int maxcap = Integer.parseInt(args[1]);
//		int elements = (args.length>=3)? Integer.parseInt(args[2]) : 1000;
//		boolean bulk = (args.length>=4) && (args[3].equalsIgnoreCase("bulk"));
//		int bufferSize = (args.length>=5)? Integer.parseInt(args[4]) : 100;
//		int queryMin = (args.length>=6)? Integer.parseInt(args[5]) : 0;
//		int queryMax = (args.length>=7)? Integer.parseInt(args[6]) : 0;


		// an unbuffered container that counts the access to the BTree
//		CounterContainer lowerCounterContainer = new CounterContainer(
//			new ConverterContainer(
//				new BlockFileContainer("BTree", 4+2+16*maxcap),
//				btree.nodeConverter(IntegerConverter.DEFAULT_INSTANCE, IntegerConverter.DEFAULT_INSTANCE, new ComparableComparator())
//			)
//		);
		// a buffered container that count the access to the buffered BTree
//		BufferedContainer bufferedContainer = new BufferedContainer(lowerCounterContainer, new LRUBuffer(bufferSize), true);
//		CounterContainer upperCounterContainer = new CounterContainer(bufferedContainer);

		// the container that stores the content of the BTree
//		Container container = upperCounterContainer;
		// Container container = new MapContainer();
		
		// initialize the BTree with the descriptor-factory method, a
		// container for storing the nodes and the minimum and maximum
		// capacity of them

		Iterator it = new PermutatorDouble(reader.m_primaryIndexSize, reader.m_remedySeed);

		long t1, t2;
		t1 = System.currentTimeMillis();

		// insert an iterator of objects by inserting every single object
		// or by bulk-insertion
//		if (bulk) {
//			it = new MergeSorter(it, new ComparableComparator(), 12, 4*4096, 4*4096);
		int memorySizeForRuns = 1024 * 1024 * 10;
		it = new MergeSorter(it, new ComparableComparator(), dataSize, memorySizeForRuns, memorySizeForRuns);
		new SortBasedBulkLoading(btree, it, new Constant(container), reader.m_bulkLoadingInitialCap);
//		}
//		else {
//			while (it.hasNext())
//				btree.insert(it.next());
//		}

		t2 = System.currentTimeMillis();

		System.out.println("Time for bulk loading: "+(t2-t1));
		System.out.println("Insertion complete, height: "+btree.height()+", universe: \n"+btree.rootDescriptor());
		PropertyReader.wsm_print("Total data in leaf nodes is ", btree.numOfDataInLeafNodes);
		PropertyReader.wsm_print("Total leaf node capacity is ", btree.totalLeafNodeCapacity);
		PropertyReader.wsm_print("Initial leaf node utilization is ", btree.numOfDataInLeafNodes/btree.totalLeafNodeCapacity);
		PropertyReader.wsm_print("Single leaf node capacity is ", btree.singleLeafNodeCapacity);
		PropertyReader.wsm_print("Single index node capacity is ", btree.singleIndexNodeCapacity);
		PropertyReader.wsm_print("Total leaf node number is ", btree.leafNodeNum);
//		System.out.println("\nAccessing the BufferedContainer\n"+upperCounterContainer+"\n");
//		System.out.println("Accessing the ConverterContainer and the BlockFileContainer\n"+lowerCounterContainer+"\n");

		System.out.println("Reset counters");
//		upperCounterContainer.reset();
//		lowerCounterContainer.reset();
		container.reset();
		container.flush();
		
//		System.out.println("Flushing buffers");
//		bufferedContainer.flush();

//		System.out.println("\nAccessing the BufferedContainer\n"+upperCounterContainer+"\n");
//		System.out.println("Accessing the ConverterContainer and the BlockFileContainer\n"+lowerCounterContainer+"\n");
//
//		System.out.print("Checking descriptors... ");
//		btree.checkDescriptors();
//
//		System.out.println("done.\n");
//
//		System.out.println("Reset counters");
//		upperCounterContainer.reset();
//		lowerCounterContainer.reset();
//
//		System.out.println("Performing Query");
//		t1 = System.currentTimeMillis();
//		// perform a range-query
//		int hits = Cursors.count(
//			btree.query(new Interval1D(new Integer(queryMin), new Integer(queryMax)))
//		);
//		t2 = System.currentTimeMillis();
//
//		System.out.println("Time for queries: "+(t2-t1));
//		System.out.println("Number of hits: "+hits);
//
		PropertyReader.wsm_print("Begin inserting data and counting splits");

		long begin = System.nanoTime();
		Random random = new Random(reader.m_insertionSeed);
		double nanoToMilliSeconds = 1000 * 1000;
		int leafSplitCount, nonLeafSplitCount;
		int numBatches = reader.m_numBatches;
		int batchSize = reader.m_batchSize;
		String[] outputFileNames = {
				"leafSplits", "nonLeafSplits", "leafSplits2", "nonLeafSplits2"};
		ArrayList<PrintWriter> pwArray = new ArrayList<>();
		for (int i = 0; i < outputFileNames.length; i++) {
			PrintWriter tempPW = new PrintWriter(reader.m_outputDir + reader.m_treeName + "_" +
					outputFileNames[i] + ".txt");
			pwArray.add(tempPW);
		}
		((BufferedContainer) treeContainer).flush();
		((BufferedContainer) treeContainer).clearBuffers();

		for (int j = 0; j < numBatches; j++) {
			if (j % 50 == 0)
				PropertyReader.wsm_print("Currently running " + j + "th batch.");
			leafSplitCount = 0;
			nonLeafSplitCount = 0;
			container.reset();
			for (int k = 0; k < batchSize; k++) {
				btree.insert(random.nextDouble());
				leafSplitCount += btree.leafSplitCount;
				nonLeafSplitCount += btree.nonLeafSplitCount;
			}
			pwArray.get(0).printf("%d\n", leafSplitCount);
			pwArray.get(1).printf("%d\n", nonLeafSplitCount);
			pwArray.get(2).printf("%d\n", container.insertPredicates);
			pwArray.get(3).printf("%d\n", container.inserts - container.insertPredicates);
			for (int i = 0; i < outputFileNames.length; i++) {
				pwArray.get(i).flush();
			}
		}
		long done = System.nanoTime();
		PropertyReader.wsm_print("Done in ", (done - begin) / (1000.0*1000*1000*60), " min.");

		System.out.println("Closing application");
		container.close();
	}
}
