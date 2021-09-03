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

package xxl.core.indexStructures;

import java.io.*;
import java.lang.management.ManagementFactory;
import java.util.*;

import xxl.core.collections.MapEntry;
import xxl.core.collections.ReversedList;
import xxl.core.collections.containers.Container;
import xxl.core.cursors.AbstractCursor;
import xxl.core.cursors.Cursor;
import xxl.core.cursors.Cursors;
import xxl.core.cursors.filters.Filter;
import xxl.core.cursors.sources.EmptyCursor;
import xxl.core.cursors.sources.SingleObjectCursor;
import xxl.core.cursors.unions.Sequentializer;
import xxl.core.functions.AbstractFunction;
import xxl.core.functions.Constant;
import xxl.core.functions.Function;
import xxl.core.io.converters.ConvertableConverter;
import xxl.core.io.converters.Converter;
import xxl.core.io.converters.MeasuredConverter;
import xxl.core.predicates.AbstractPredicate;
import xxl.core.spatial.points.DoublePoint;
import xxl.core.spatial.points.Point;
import xxl.core.spatial.rectangles.DoublePointRectangle;
import xxl.core.spatial.rectangles.Rectangle;

/**
 * 
 * HilbertRTree with R*tree split(1-2 split policy) and BPlusTree merge
 * strategies. This Implementation is based on BPlusTree. The indexEntries of
 * BPlusTree were extended to manage MBBs. To use spatial semantic of the tree
 * queries run {@link HilbertRTree#queryOR(Rectangle, int)} For a detailed
 * discussion see Ibrahim Kamel, Christos Faloutsos :
 * "Hilbert R-tree: An Improved R-tree Using Fractals" Proceedings of the 20th
 * VLDB Conference Sanitago, Chile, 1994
 * 
 * @see Tree
 * @see BPlusTree
 * @see ORTree
 * @see RTree
 */
public class HilbertRTree extends BPlusTree {
	private int count = 1;
	private int interval = 500;
	private Random rand = new Random(21);
	public String outDir;


	/** Dimension of the data */
	protected int dimension;
	/**
	 * minMaxFactor the quotient between minimum and maximum number of entries
	 * in an node, e.g. 0.5 (BPlusTRee) standard Split or e.g. 1d/3d R*Split
	 */
	protected double minMaxFactor;
	/** Function to compute MBB from data */
	protected Function getRectangle;
	/** Converter to serialize MBBs */
	protected Converter descriptorConverter;
	/** Function to compute middle point of the MBB */
	protected Function computeMiddlePoint;
	/** Normalize function to hypercube [0..1) */
	protected Function normalize;
	/** Universe-hypercube of the data */
	protected Rectangle universe;


	/**
	 * Creates a new <tt>HilbertRTree</tt>. With a default setting of duplicates
	 * = true, and minMaxFactor 0.5 (standard BPlusTree split)
	 * 
	 * @param blockSize
	 * @param universe
	 *            Universe hypercube of the data
	 */
	public HilbertRTree(int blockSize, Rectangle universe) {
		this(blockSize, universe, 0.5d, true, 1);
	}

	/**
	 * Creates a new <tt>HilbertRTree</tt>. With a default setting minMaxFactor
	 * 0.5 (standard BPlusTree split)
	 * 
	 * @param blockSize
	 * @param universe
	 *            Universe hypercube of the data
	 * @param duplicate
	 */
	public HilbertRTree(int blockSize, Rectangle universe, boolean duplicate) {
		this(blockSize, universe, 0.5d, duplicate, 1);
	}

	/**
	 * Creates a new <tt>HilbertRTree</tt>. With a default setting of duplicates
	 * = true.
	 * 
	 * @param blockSize
	 * @param universe
	 *            Universe hypercube of the data
	 * @param minMaxFactor
	 *            minMaxFactor the quotient between minimum and maximum number
	 *            of entries in an node, e.g. 0.5 (BPlusTRee) standard Split or
	 *            e.g. 1d/3d R*Split
	 */
	public HilbertRTree(int blockSize, Rectangle universe, Double minMaxFactor) {
		this(blockSize, universe, minMaxFactor, true, 1);
	}

	/**
	 * Created by LX. The first three param is the same as above
	 * @param blockSize
	 * @param universe
	 * @param minMaxFactor
	 * @param siblingNum: number 2 indicate a deferred split of 2-to-3
	 */
	public HilbertRTree(int blockSize, Rectangle universe, Double minMaxFactor, Integer siblingNum) {
		this(blockSize, universe, minMaxFactor, true, siblingNum);
	}

	/**
	 * Creates a new <tt>HilbertRTree</tt>.
	 * 
	 * @param blockSize
	 * @param universe
	 *            Universe hypercube of the data
	 * @param minMaxFactor
	 *            minMaxFactor the quotient between minimum and maximum number
	 *            of entries in an node, e.g. 0.5 (BPlusTRee) standard Split or
	 *            e.g. 1d/3d R*Split
	 * @param duplicate
	 * @param siblingNum: added by LX. To distinguish 1-to-2 from 2-to-3
	 */
	public HilbertRTree(int blockSize, Rectangle universe, double minMaxFactor,
			boolean duplicate, int siblingNum) {
		super(blockSize, minMaxFactor, siblingNum);
		this.duplicate = duplicate;
		this.minMaxFactor = minMaxFactor;
		this.universe = universe;
		this.dimension = universe.dimensions();
//		this.splitCount = 0;
	}

	/**
	 * Initializes the <tt>HilbertRTree</tt>.
	 * 
	 * @param computeSFCurveValue
	 *            Function which computes SpaceFillingCurve value of the middle
	 *            point of the normalized dataMBR
	 * @param getMBR
	 *            function which computes the MBR of the dataObject
	 * @param container
	 *            container of the tree
	 * @param keyConverter
	 *            converter for keys of the dataObjects
	 * @param dataConverter
	 * @param createORSeparator
	 *            a factory Function to create (Overlaping)Separators
	 * @param createKeyRange
	 *            a factory Function to create (Overlaping)KeyRanges
	 * @return the initialized <tt>HilbertRTree</tt> itself
	 */
	public HilbertRTree initialize(final Function computeSFCurveValue,
			Function getMBR, Container container,
			MeasuredConverter keyConverter, MeasuredConverter dataConverter,
			Function createORSeparator, Function createKeyRange) {
		return initialize(null, null, computeSFCurveValue, getMBR, container,
				keyConverter, dataConverter, createORSeparator, createKeyRange);
	}

	/**
	 * Initializes the <tt>HilbertRTree</tt>.
	 * 
	 * @param computeSFCurveValue
	 *            computeSFCurveValue Function which computes SpaceFillingCurve
	 *            value of the middle Point of the normalized dataMBR
	 * @param getMBR
	 *            function which computes the MBR of the dataObject
	 * @return the initialized <tt>HilbertRTree</tt> itself
	 */
	protected HilbertRTree initialize(final Function computeSFCurveValue,
			Function getMBR) {
		getRectangle = getMBR;
		normalize = new AbstractFunction() {
			public Object invoke(Object entry) {
				if (entry instanceof DoublePointRectangle) {
					DoublePointRectangle normMBR = (DoublePointRectangle) ((DoublePointRectangle) entry)
							.clone();
					return normMBR.normalize(universe);
				} else {
					System.out.println("entry " + entry);
					throw new IllegalArgumentException();
				}
			}
		};
		descriptorConverter = new ConvertableConverter(new AbstractFunction() {
			public Object invoke() {
				return new DoublePointRectangle(dimension);
			}
		});
		// Normalize and compute middle Point
		computeMiddlePoint = new AbstractFunction() {
			public Object invoke(Object entry) {
				Rectangle rectangle = (Rectangle) normalize.invoke(entry);
				double[] pointArray = new double[rectangle.dimensions()];
				for (int i = 0; i < pointArray.length; i++) {
					pointArray[i] = (rectangle.getCorner(true).getValue(i) + rectangle
							.getCorner(false).getValue(i)) / 2d;
				}
				return new DoublePoint(pointArray);
			}
		}; // Edit
		this.getKey = new AbstractFunction() {
			public Object invoke(Object entry) {
				Rectangle rectangle = (Rectangle) getRectangle.invoke(entry);
				Point middlePoint = (Point) computeMiddlePoint
						.invoke(rectangle);
//				System.out.println("(" + middlePoint + ") => " + computeSFCurveValue.invoke(middlePoint));
				return computeSFCurveValue.invoke(middlePoint);
			}
		};
		this.getDescriptor = new AbstractFunction() {
			public Object invoke(Object entry) {
//				System.out.println(entry);
				if (entry instanceof ORSeparator)
					return (ORSeparator) entry;
				if (entry instanceof BPlusTree.IndexEntry) // LX
					return ((BPlusTree.IndexEntry) entry).separator; // LX
				return createORSeparator(key(entry),
						(Descriptor) getRectangle.invoke(entry));
			}
		};
		return this;
	}

	/**
	 * This initialization method mainly is used to restore the persistent tree.
	 * 
	 * @param rootEntry
	 *            rootEntry of the tree
	 * @param rootDescriptor
	 *            the ORKeyRange of the tree {@link HilbertRTree.ORKeyRange}
	 * @param computeSFCurveValue
	 *            Function which computes SpaceFillingCurve value of the middle
	 *            point of the dataMBR
	 * @param getMBR
	 *            function which computes the MBR of the dataObject
	 * @param container
	 *            container of the tree
	 * @param keyConverter
	 *            converter for keys of the dataObjects
	 * @param dataConverter
	 * @param createORSeparator
	 *            a factory Function to create (overlapping)Separators
	 * @param createKeyRange
	 *            a factory Function to create (overlapping)KeyRanges
	 * @return the initialized <tt>HilbertRTree</tt> itself
	 */
	public HilbertRTree initialize(IndexEntry rootEntry,
			Descriptor rootDescriptor, final Function computeSFCurveValue,
			Function getMBR, Container container,
			MeasuredConverter keyConverter, MeasuredConverter dataConverter,
			Function createORSeparator, Function createKeyRange) {
		super.initialize(rootEntry, rootDescriptor, getKey, container,
				keyConverter, dataConverter, createORSeparator, createKeyRange,
				new Constant(minMaxFactor), new Constant(
						((1.0 - minMaxFactor) <= minMaxFactor) ? minMaxFactor
								: 1.0 - minMaxFactor));
		return initialize(computeSFCurveValue, getMBR);

	}

	/**
	 * Creates a new node on a given level.
	 * 
	 * @param level
	 *            the level of the new Node
	 * 
	 * @see xxl.core.indexStructures.Tree#createNode(int)
	 */
	public Tree.Node createNode(int level) {
		Node node = new Node(level);
		return node;
	}

	/**
	 * returns the space filling curve value of data object ( e.g. Hilbert-Curve
	 * value of the MBRs middle point)
	 * 
	 * @param data
	 * @return
	 */
	public Comparable getSFCValue(Object data) {
		return key(data);
	}

	/**
	 * returns comparator based on the functions of the initialized tree
	 * 
	 * @return comparator for universe of the tree
	 */
	public Comparator getDefaultComparator() {
		return new Comparator() {
			public int compare(Object o1, Object o2) {
				ORSeparator s1 = orSeparator(o1);
				ORSeparator s2 = orSeparator(o2);
				return s1.compareTo(s2);
			}
		};
	}

	/**
	 * Returns ORSeparator of the data entry
	 * 
	 * @param entry
	 * @return ORSeparator
	 */
	protected ORSeparator orSeparator(Object entry) {
		return (ORSeparator) separator(entry);
	}

	/**
	 * 
	 * @return descriptor of the entry
	 */
	public Descriptor descriptor(Object entry) {
		return (Descriptor) getDescriptor.invoke(entry);
	}

	/**
	 * creates overlapping separator
	 * 
	 * @param sepValue
	 * @param mbr
	 * @return
	 */
	protected ORSeparator createORSeparator(Comparable sepValue, Descriptor mbr) {
		return (ORSeparator) this.createSeparator.invoke(sepValue, mbr);
	}

	/**
	 * creates overlapping KeyRange
	 * 
	 * @param sepValue
	 * @param mbr
	 * @return
	 */
	protected ORKeyRange createORKeyRange(Comparable sepValue,
			Comparable maxBound, Rectangle mbr) {
		return (ORKeyRange) this.createKeyRange.invoke(Arrays
				.asList(new Object[] { sepValue, maxBound, mbr }));
	}

	/**
     * 
     */
	protected KeyRange createKeyRange(Comparable min, Comparable max) {
		Rectangle rect;
		if (rootDescriptor == null) {
			rect = (Rectangle) ((Descriptor) universe).clone();
		} else {
			rect = (Rectangle) ((ORKeyRange) rootDescriptor).getIndexEntryMBR()
					.clone();
		}

		return createORKeyRange(min, max, rect);
	}

	/**
	 * computes Descriptor(MBB) of the entry
	 * 
	 * @param entry
	 * @return
	 */
	protected Rectangle rectangle(Object entry) {
		// Modified by LX
		if (entry instanceof HilbertRTree.IndexEntry) {
			return (Rectangle) (((IndexEntry) entry).hilbertEntryMBR).clone();
		} else if (entry instanceof BPlusTree.IndexEntry) {
			ORSeparator separator = (ORSeparator) separator(entry);
			return (Rectangle) ((Descriptor) separator.getIndexEntryMBR()).clone();
		}
		else {
			return (Rectangle) getRectangle.invoke(entry);
		}
//		ORSeparator separator = (ORSeparator) separator(entry);
//		return (Rectangle) ((Descriptor) separator.getIndexEntryMBR()).clone();
	}

	/**
	 * Computes the union of the descriptors(MBRs) of the collection's entries.
	 * 
	 * @param entries
	 *            a collection of objects
	 * @return the union of the MBRs of the collection's entries
	 * 
	 */
	public Rectangle computeMBR(Iterator entries) {
		Rectangle descriptor = null;
		if (entries.hasNext()) {
			descriptor = (Rectangle) ((Descriptor) rectangle(entries.next()))
					.clone();
			while (entries.hasNext())
				descriptor.union(rectangle(entries.next()));
		}
		return descriptor;
	}

	/**
     * 
     */
	protected NodeConverter createNodeConverter() {
		return new NodeConverter();
	}

	/**
     * 
     */
	protected void insert(Object data, Descriptor descriptor, int targetLevel) {
//		this.splitCount = 0;
		if (rootEntry() == null) {
			Comparable key = ((Separator) descriptor).sepValue();
			Rectangle indexEntryMBR = ((ORSeparator) descriptor)
					.getIndexEntryMBR();

			Rectangle mbr = (Rectangle) ((Descriptor) indexEntryMBR).clone();
			rootDescriptor = createORKeyRange(key, key, mbr);
			grow(data);
		} else
			super.insert(data, descriptor, targetLevel);
	}

	/**
	 * Original Algorithm from ORTree The cursor does not support remove
	 * operations. To remove the object from the tree call
	 * {@link HilbertRTree#remove(Object)}
	 * 
	 * @see ORTree#query()
	 * @param queryRectangle
	 * @param targetLevel
	 * @return
	 */
	public Cursor queryOR(final Rectangle queryRectangle, final int targetLevel) {

		final Iterator[] iterators = new Iterator[height() + 1];
		Arrays.fill(iterators, EmptyCursor.DEFAULT_INSTANCE);
		if (height() > 0
				&& queryRectangle.overlaps(((ORKeyRange) rootDescriptor())
						.getIndexEntryMBR()))
			iterators[height()] = new SingleObjectCursor(rootEntry());
		return new AbstractCursor() {
			int queryAllLevel = 0;
			Object toRemove = null;
			Stack path = new Stack();

			public boolean hasNextObject() {
				for (int parentLevel = targetLevel;;) {
					if (iterators[parentLevel].hasNext()) {
						if (parentLevel == targetLevel)
							return true;
						else {
							IndexEntry indexEntry = (IndexEntry) iterators[parentLevel]
									.next();
							if (indexEntry.level() >= targetLevel) {
//								long start = System.nanoTime();
								Tree.Node node = indexEntry.get(true);
//								long end = System.nanoTime();
//								System.out.println("Query get takes time: " + (end-start));

								Iterator queryIterator;
								if (parentLevel <= queryAllLevel
										|| queryRectangle
										.contains(((ORSeparator) indexEntry
												.separator())
												.getIndexEntryMBR())) {
									queryIterator = node.entries(); // Falls
									// alle
									// entries
									if (parentLevel > queryAllLevel
											&& !iterators[node.level].hasNext())
										queryAllLevel = node.level;
								} else
									// edit
									queryIterator = ((Node) node)
											.queryOR(queryRectangle);
								iterators[parentLevel = node.level] = iterators[parentLevel]
										.hasNext() ? new Sequentializer(
										queryIterator, iterators[parentLevel])
										: queryIterator;
								path.push(new MapEntry(indexEntry, node));
							}
						}
					} else if (parentLevel == height()) {
						return false;
					} else {
						if (parentLevel == queryAllLevel)
							queryAllLevel = 0;
						if (level(path) == parentLevel)
							path.pop();
						iterators[parentLevel++] = EmptyCursor.DEFAULT_INSTANCE;
					}
				}
			}

			public Object nextObject() {
				return toRemove = iterators[targetLevel].next();
			}

			public void remove() throws UnsupportedOperationException {
				throw new UnsupportedOperationException();
			}

		};
	}

	/**
	 * Added by LX
	 */
	public int queryNum = 1;
	public Cursor queryOR(final Rectangle queryRectangle, final int targetLevel, final int[] nodeAccess) {
		final Iterator[] iterators = new Iterator[height() + 1];
		Arrays.fill(iterators, EmptyCursor.DEFAULT_INSTANCE);
		if (height() > 0
				&& queryRectangle.overlaps(((ORKeyRange) rootDescriptor())
				.getIndexEntryMBR()))
			iterators[height()] = new SingleObjectCursor(rootEntry());
		queryNum++;

		return new AbstractCursor() {
			int queryAllLevel = 0;
			Object toRemove = null;
			Stack path = new Stack();

			public boolean hasNextObject() {
				for (int parentLevel = targetLevel;;) {
					if (iterators[parentLevel].hasNext()) {
						if (parentLevel == targetLevel) {
							return true;
						} else {
							IndexEntry indexEntry = (IndexEntry) iterators[parentLevel].next();
							if (indexEntry.level() >= targetLevel) {

								Tree.Node node = indexEntry.get(false); // modified by LX

								Iterator queryIterator;

								if (parentLevel <= queryAllLevel
										|| queryRectangle.contains(indexEntry.hilbertEntryMBR) ) { // LX
//												(indexEntry.ORseparator()).getIndexEntryMBR())) {
									queryIterator = node.entries(); // Falls
									// alle
									// entries
									if (parentLevel > queryAllLevel && !iterators[node.level].hasNext())
										queryAllLevel = node.level;
								} else {
									// edit
									queryIterator = ((Node) node).queryOR(queryRectangle);
								}

								iterators[parentLevel = node.level] = iterators[parentLevel].hasNext() ?
										new Sequentializer(queryIterator, iterators[parentLevel]) : queryIterator;
								path.push(new MapEntry(indexEntry, node));
							}
						}
					} else if (parentLevel == height()) {
						return false;
					} else {
						if (parentLevel == queryAllLevel)
							queryAllLevel = 0;
						if (level(path) == parentLevel) {
//							path.pop();
							up(path); // LX
						}
						iterators[parentLevel++] = EmptyCursor.DEFAULT_INSTANCE;
					}
				}
			}

			public Object nextObject() {
				return toRemove = iterators[targetLevel].next();
			}

			public void remove() throws UnsupportedOperationException {
				throw new UnsupportedOperationException();
			}
		};
	}

	/**
	 * This method is an implementation of an efficient querying algorithm. The
	 * result is a lazy cursor pointing to all leaf entries whose descriptors
	 * overlap with the given <tt>queryDescriptor</tt>.
	 * 
	 * @param queryDescriptor
	 *            describes the query in terms of a descriptor
	 * @return a lazy <tt>Cursor</tt> pointing to all response objects
	 * @see HilbertRTree#query(Descriptor, int)
	 */
	public Cursor queryOR(Rectangle queryDescriptor) {
		return queryOR(queryDescriptor, 0);
	}

	/**
	 * This method executes a query using the rootDescriptor on a given level.
	 * That means, that the response consists of all entries of the given level.
	 * 
	 * @param level
	 *            the target level of the query
	 * @return a lazy <tt>Cursor</tt> pointing to all response objects
	 */
	public Cursor queryOR(int level) {
		return queryOR(((ORKeyRange) rootDescriptor()).entryMBR, level);
	}

	/**
	 * This method executes a query using the rootDescriptor on the leaf level.
	 * That means, that the response consists of all leaf entries (i.e. data
	 * objects) stored in the tree.
	 * 
	 * @return a lazy Cursor pointing to all response objects
	 */
	public Cursor queryOR() {
		int[] nodeAccess = {0, 0};
		return queryOR(((ORKeyRange) rootDescriptor()).entryMBR, 0, nodeAccess);
	}

	// LX
	public Tree.IndexEntry createIndexEntry(int parentLevel) {
		return new IndexEntry(parentLevel);
	}

	public class IndexEntry extends BPlusTree.IndexEntry {

		public ORSeparator ORDescriptor;

		public Rectangle hilbertEntryMBR;

		public IndexEntry(int level) { super(level); }



	}

	/**
	 * @see BPlusTree.Node
	 * 
	 */
	public class Node extends BPlusTree.Node {
		/**
		 * 
		 * @param level
		 */
		public Node(int level) {
			super(level);
		}

		/**
		 * @see BPlusTree.Node#initialize(Object)
		 */
		public Tree.Node.SplitInfo initialize(Object entry) {
			Stack path = new Stack();
			Rectangle entryMBR = ((ORKeyRange) rootDescriptor).entryMBR;

			if (height() > 0) {
				IndexEntry indexEntry = (IndexEntry) entry;
				if (indexEntry == rootEntry) {
					indexEntry.separator = createORSeparator(
							((ORKeyRange) rootDescriptor).maxBound,
							(Rectangle) ((Descriptor) entryMBR).clone());
				}
			}
			grow(entry, path);
			// separator of new root entry in tree.grow() method
			Comparable bound = ((ORKeyRange) rootDescriptor).maxBound;
			return new SplitInfo(path).initialize(createORSeparator(bound,
					(Rectangle) ((Descriptor) entryMBR).clone()));
		}

		/**
		 * Returns an iterator pointing to entries whose MBRs overlap
		 * <tt>queryMBR</tt>.
		 * 
		 * @param queryRectangle
		 *            the descriptor describing the query
		 * @return an iterator pointing to entries whose MBRs overlap with
		 *         <tt>queryMBR</tt>
		 */
		public Iterator queryOR(final Rectangle queryRectangle) {
			return new Filter(entries(), new AbstractPredicate() {
				public boolean invoke(Object object) {
					return rectangle(object).overlaps(queryRectangle);
				}
			});
		}

		/**
		 * @see BPlusTree.Node#chooseSubtree(Descriptor, Stack)
		 */
		protected Tree.IndexEntry chooseSubtree(Descriptor descriptor,
				Stack path) {
			IndexEntry indexEntry = (IndexEntry) super.chooseSubtree(
					descriptor, path);
			ORSeparator entryDescriptor = (ORSeparator) descriptor;
			ORSeparator indexEntryDescriptor = (ORSeparator) indexEntry
					.separator();
			if (!indexEntryDescriptor.containsMBR(entryDescriptor.getIndexEntryMBR())) {
				indexEntryDescriptor.unionMBR(entryDescriptor);
				update(path);
			}
			return indexEntry;
		}

		/**
		 * Added by LX. Used in deferred split of Hilbert R-tree only.
		 */
		protected Tree.IndexEntry chooseSubtree(Descriptor descriptor,
												Stack path, Stack siblingStack) {
			IndexEntry indexEntry = (IndexEntry) super.chooseSubtree(
					descriptor, path, siblingStack);
			ORSeparator entryDescriptor = (ORSeparator) descriptor;
			ORSeparator indexEntryDescriptor = (ORSeparator) indexEntry
					.separator();
			if (!indexEntryDescriptor.containsMBR(entryDescriptor.getIndexEntryMBR())) {
				indexEntryDescriptor.unionMBR(entryDescriptor);
				update(path);
			}
			return indexEntry;
		}

		/**
		 * This method always returns false.
		 * 
		 * @see BPlusTree.Node#redistributeLeaf(Stack)
		 * 
		 */
		protected boolean redistributeLeaf(Stack path) {
			return false;
		}

		/**
		 * Added by LX
		 * Implements 2-to-3 split
		 * (*this) refers to the newly created node
		 * @param path
		 * @return
		 */
		public Tree.Node.SplitInfo deferSplit(Stack path, Stack siblingStack) {
			Node node = (Node) node(path);
			Tree.IndexEntry indexEntry = indexEntry(path);

			// Find the cooperative sibling
			Tree.IndexEntry coopSibling = (Tree.IndexEntry) ((Map.Entry) siblingStack.peek()).getValue();
			boolean rightSibling = (boolean)((Map.Entry)siblingStack.peek()).getKey();
			Node coopNode = (HilbertRTree.Node) coopSibling.get();

			int desiredEntries = (node.number() + coopNode.number()) / 3;

			List curEntries;
			if (rightSibling) {
				curEntries = node.entries;
				curEntries.addAll(coopNode.entries);
			} else {
				curEntries = coopNode.entries;
				curEntries.addAll(node.entries);
			}
//			System.out.printf("Total entry num %d, desired %d\n", curEntries.size(), desiredEntries);
			Rectangle fstDescriptor, sndDescriptor, thirdDescriptor;
			int count = 0;
			Iterator iter = curEntries.iterator();
			fstDescriptor = new DoublePointRectangle(rectangle(iter.next()));
			count++;
			while (iter.hasNext()) {
				if (count >= desiredEntries)
					break;
				fstDescriptor.union(rectangle(iter.next()));
				count++;
			}

			sndDescriptor = new DoublePointRectangle(rectangle(iter.next()));
			count++;
			while (iter.hasNext()) {
				if (count >= 2 * desiredEntries)
					break;
				sndDescriptor.union(rectangle(iter.next()));
				count++;
			}

			thirdDescriptor = new DoublePointRectangle(rectangle(iter.next()));
			while (iter.hasNext()) {
				thirdDescriptor.union(rectangle(iter.next()));
			}

			Object[] entryArr = curEntries.toArray();
			node.entries.clear();
			coopNode.entries.clear();
			curEntries = Arrays.asList(entryArr.clone());
//			System.out.printf("total entry num %d, desired %d\n", curEntries.size(), desiredEntries);
			if (rightSibling) {
				node.entries.addAll(curEntries.subList(0, desiredEntries));
				coopNode.entries.addAll(curEntries.subList(desiredEntries, 2 * desiredEntries));
				((ORSeparator) ((BPlusTree.IndexEntry) indexEntry).separator).entryMBR = (Rectangle) fstDescriptor.clone();
				((ORSeparator) ((BPlusTree.IndexEntry) coopSibling).separator).entryMBR = (Rectangle) sndDescriptor.clone();
				adjustSeparatorValue(node, (IndexEntry) indexEntry, path);
				adjustSeparatorValue(coopNode, (IndexEntry) coopSibling, path);
			} else {
				coopNode.entries.addAll(curEntries.subList(0, desiredEntries));
				node.entries.addAll(curEntries.subList(desiredEntries, 2 * desiredEntries));
				((ORSeparator) ((BPlusTree.IndexEntry) coopSibling).separator).entryMBR = (Rectangle) fstDescriptor.clone();
				((ORSeparator) ((BPlusTree.IndexEntry) indexEntry).separator).entryMBR = (Rectangle) sndDescriptor.clone();
			}
			adjustSeparatorValue(node, (IndexEntry) indexEntry, path);
			adjustSeparatorValue(coopNode, (IndexEntry) coopSibling, path);
			entries.addAll(curEntries.subList(2 * desiredEntries, curEntries.size()));
			ORSeparator sepNewNode = (ORSeparator) separator(this.getLast())
					.clone();
			sepNewNode.entryMBR = (Rectangle) thirdDescriptor.clone();
			Tree.Node.SplitInfo spi = (new SplitInfo(path)).initialize(sepNewNode);
			coopSibling.update(coopNode, true);
//			System.out.println("After split: " + node.entries.size() +  "," + coopNode.entries.size() + "," + entries.size());
			return spi;
		}

		/**
		 * 
		 * This method is used to split the first <tt>Node</tt> on the path in
		 * two <tt>Nodes</tt>. The current <tt>Node</tt> should be the empty new
		 * <tt>Node</tt>. The method distributes the entries of overflowing
		 * <tt>Node</tt> to the both <tt>Nodes</tt>. The Split-Algorithm (This
		 * is equivalent to the R*Tree SplitAlgorithm after the suitable axis
		 * was choosen) tests all distributions (M-2m+2) according the
		 * SpaceFillingCurve order and chooses distribution with minimal
		 * overlap-, margin- and area Value.
		 * 
		 * @param path
		 *            the <tt>Node</tt> already visited during this insert
		 * @return a <tt>SplitInfo</tt> containig all needed information about
		 *         the split
		 */

		protected Tree.Node.SplitInfo split(Stack path) {
			Node node = (Node) node(path);

			Distribution distribution;
			List<Distribution> distributionsList = new ArrayList<Distribution>();
			int minEntries = node.splitMinNumber();
			int maxEntries = node.splitMaxNumber();
			// LX: spread split
			if (node.number() < minEntries) {
//				cannot be in here
				minEntries = node.number() / 2;
				maxEntries = node.number() - minEntries;
			}

			maxEntries = (maxEntries < node.number() && maxEntries > minEntries) ? maxEntries
					: minEntries;

			// copied from rtree split algorithmus
			Rectangle[][] rectangles = new Rectangle[2][];

			// LX: if this is RES remedy,
			// need to add the entries from the overflow nodes
//			boolean hasOverflowPage = false;
			if (HilbertRTree.this.regularElectiveSplit) {
				if (level() == 0 && HilbertRTree.this.pseudoOverflowNodes.containsKey(node.nodeId)) {
					List overflowPages = (List) HilbertRTree.this.pseudoOverflowNodes.get(node.nodeId);
					for (Object list: overflowPages) {
						for (Object entry: (List) list) {
							node.grow(entry);
						}
						HilbertRTree.this.totalOverflowEntryNum -= ((List<?>) list).size();
					}
					HilbertRTree.this.totalOverflowPageNum -= overflowPages.size();
					HilbertRTree.this.pseudoOverflowNodes.remove(node.nodeId);
//					hasOverflowPage = true;
				}
			}

			for (int k = 0; k < 2; k++) {
				Iterator entries = (k == 0 ? node.entries : new ReversedList(node.entries)).iterator();

				Rectangle rectangle = new DoublePointRectangle(
						rectangle(entries.next()));
				for (int l = (k == 0 ? minEntries : node.number() - maxEntries); --l > 0;) {
					rectangle.union(rectangle(entries.next()));
				}
				(rectangles[k] = new Rectangle[maxEntries - minEntries + 1])[0] = rectangle;
				for (int j = 1; j <= maxEntries - minEntries; rectangles[k][j++] = rectangle) {

					rectangle = Descriptors.union(rectangle,
							rectangle(entries.next()));
				}
			}


			// Creation of the distributions
			for (int j = minEntries; j <= maxEntries; j++)
				distributionsList.add(new Distribution(node.entries.toArray(),
						j, rectangles[0][j - minEntries],
						rectangles[1][maxEntries - j]));

			if (unequalRandomSplit) {
				int idx = rand.nextInt(distributionsList.size());
				distribution = (Distribution) distributionsList.get(idx);
			}
			else {
				// Choose the distributions of the chosen dimension with minimal
				// overlap
				distributionsList = Cursors.minima(distributionsList.iterator(),
						new AbstractFunction() {
							public Object invoke(Object object) {
								return new Double(((Distribution) object)
										.overlapValue());
							}
						});
				// If still more than one distribution has to be considered, choose
				// one
				// with minimal perimeter
				distributionsList = Cursors.minima(distributionsList.iterator(),
						new AbstractFunction() {
							public Object invoke(Object object) {
								// return new
								// Double(((Distribution)object).areaValue());
								return new Double(((Distribution) object)
										.marginValue());
							}
						});
				distributionsList = Cursors.minima(distributionsList.iterator(),
						new AbstractFunction() {
							public Object invoke(Object object) {
								return new Double(((Distribution) object)
										.areaValue());
								// return new
								// Double(((Distribution)object).marginValue());
							}
						});
				distribution = (Distribution) distributionsList.get(0);
			}

			node.entries.clear();

			node.entries.addAll(distribution.entries(false));
			entries.addAll(distribution.entries(true));

			ORSeparator sepNewNode = (ORSeparator) separator(this.getLast())
					.clone();
			sepNewNode.entryMBR = (Rectangle) distribution.descriptor(true)
					.clone();
			// update the descriptor of the old index entry
			((ORSeparator) (((BPlusTree.IndexEntry) (indexEntry(path))).separator)).entryMBR = (Rectangle) distribution
					.descriptor(false).clone();
			((IndexEntry) (indexEntry(path))).hilbertEntryMBR = (Rectangle) distribution.descriptor(false).clone();
//			((IndexEntry) (indexEntry(path))).ORDescriptor = (ORSeparator)(((BPlusTree.IndexEntry) (indexEntry(path))).separator);
//			return (new SplitInfo(path)).initialize(sepNewNode);
			Tree.Node.SplitInfo spi = (new SplitInfo(path)).initialize(sepNewNode);
//			spi.splitOverlap = distribution.overlapValue();
//			spi.areaDiff = distribution.areaDiffValue();
//			spi.splitOverlapRatio = spi.splitOverlap / (Math.min(distribution.firstDescriptor.area(), distribution.secondDescriptor.area()));
//			spi.splitOverlapRatio2 = spi.splitOverlap / (Math.max(distribution.firstDescriptor.area(), distribution.secondDescriptor.area()));
//			spi.splitRatio = (Math.min(distribution.firstDescriptor.area(), distribution.secondDescriptor.area())) /
//					(Math.max(distribution.firstDescriptor.area(), distribution.secondDescriptor.area()));
			// LX: RES remedy
			// need to assign a node id to the newly node
			if (HilbertRTree.this.regularElectiveSplit && node.level() == 0) {
				long nid = node.nodeId;
				long nnid = HilbertRTree.this.leafNodeNum + 1;
				nodeId = nnid;
				Iterator it = node.entries();
				if (it.hasNext()) {
					Object obj = it.next();
					if (obj != null)
						HilbertRTree.this.leafNodeToRange.put(nid, HilbertRTree.this.getDescriptor.invoke(obj));
				}
				it = entries();
				if (it.hasNext()) {
					Object obj = it.next();
					if (obj != null)
						HilbertRTree.this.leafNodeToRange.put(nodeId, HilbertRTree.this.getDescriptor.invoke(obj));
				}
			}
			return spi;
		}

		/**
		 * Added by LX
		 * For deferred split only:
		 * if the node overflows but its cooperative sibling is not full,
		 * redistribute entries evenly among these two nodes (assume 2-to-3 split)
		 */
		boolean rearrangeEntriesWithSibling(Stack path, Stack siblingStack) {
			Tree.IndexEntry parent = indexEntry(path);
			// Find the cooperative sibling
			Tree.IndexEntry coopSibling = (Tree.IndexEntry) ((Map.Entry) siblingStack.peek()).getValue();
			boolean rightSibling = (boolean)((Map.Entry)siblingStack.peek()).getKey();
			if (rightSibling && coopSibling == null)
				return false; // no sibling
			Node coopNode = (HilbertRTree.Node) coopSibling.get();
//			System.out.println("Before: " + ((BPlusTree.IndexEntry) parent).separator + " and " + ((BPlusTree.IndexEntry) coopSibling).separator);
//			System.out.println("Before: " + number() + " and " + coopNode.number());
			if ((level == 0 && coopNode.number() + number() > 2*singleLeafNodeCapacity) ||
					(level > 0 && coopNode.number() + number() > 2*singleIndexNodeCapacity))
				return false;
			int desiredEntries = (number() + coopNode.number()) / 2; // for the left node
//			System.out.println("Before " + ((ORSeparator) ((BPlusTree.IndexEntry) parent).separator).entryMBR + " and " +
//					((ORSeparator) ((BPlusTree.IndexEntry) coopSibling).separator).entryMBR);

			List curEntries;
			if (rightSibling) {
				curEntries = this.entries;
				curEntries.addAll(coopNode.entries);
			} else {
				curEntries = coopNode.entries;
				curEntries.addAll(this.entries);
			}

			Rectangle fstDescriptor, sndDescriptor;
			int count = 0;

			Iterator iter = curEntries.iterator();
			fstDescriptor = new DoublePointRectangle(rectangle(iter.next()));
			while (iter.hasNext()) {
				fstDescriptor.union(rectangle(iter.next()));
				count++;
				if (count >= desiredEntries)
					break;
			}

			sndDescriptor = new DoublePointRectangle(rectangle(iter.next()));
			while (iter.hasNext()) {
				sndDescriptor.union(rectangle(iter.next()));
			}
//			System.out.println("Distribution " + fstDescriptor + " and " + sndDescriptor);
			Distribution distribution = new Distribution(curEntries.toArray(), desiredEntries, fstDescriptor, sndDescriptor);
			this.entries.clear();
			coopNode.entries.clear();
			if (rightSibling) {
				this.entries.addAll(distribution.entries(false));
				coopNode.entries.addAll(distribution.entries(true));
				((ORSeparator) ((BPlusTree.IndexEntry) parent).separator).entryMBR = (Rectangle) distribution
						.descriptor(false).clone();
//				adjustSeparatorValue(this, (IndexEntry) parent, path);
				Comparable sepLastValue = ((Separator)separator(this.getLast()).clone()).sepValue();
				((BPlusTree.IndexEntry) parent).separator.updateSepValue(sepLastValue); // update
				((ORSeparator) ((BPlusTree.IndexEntry) coopSibling).separator).entryMBR = (Rectangle) distribution
						.descriptor(true).clone();
			} else {
				this.entries.addAll(distribution.entries(true));
				coopNode.entries.addAll(distribution.entries(false));
				((ORSeparator) ((BPlusTree.IndexEntry) parent).separator).entryMBR = (Rectangle) distribution
						.descriptor(true).clone();
				((ORSeparator) ((BPlusTree.IndexEntry) coopSibling).separator).entryMBR = (Rectangle) distribution
						.descriptor(false).clone();
				Comparable sepLastValue = ((Separator)separator(coopNode.getLast()).clone()).sepValue();
				((BPlusTree.IndexEntry) coopSibling).separator.updateSepValue(sepLastValue); // update
//				adjustSeparatorValue(coopNode, (IndexEntry) coopSibling, path);
			}
//			System.out.println("After " + ((ORSeparator) ((BPlusTree.IndexEntry) parent).separator).entryMBR + " and " +
//					((ORSeparator) ((BPlusTree.IndexEntry) coopSibling).separator).entryMBR);
			coopSibling.update(coopNode, true);
//			System.out.println("After: " + ((BPlusTree.IndexEntry) parent).separator + " and " + ((BPlusTree.IndexEntry) coopSibling).separator);
//			System.out.println("After rearrange: " + number() + " and " + coopNode.number());
			return true;
		}

		/**
		 * Updates MBR after the entry was deleted but the node was not overflow
		 */
		protected boolean adjustSeparatorValue(BPlusTree.Node node,
				IndexEntry indexEntry, Stack path) {
			super.adjustSeparatorValue(node, indexEntry, path);
			if (!path.isEmpty()) {// indexEntry.separator() != null ){
				Rectangle mbr = (Rectangle) ((Descriptor) computeMBR(node
						.entries())).clone();
				((ORSeparator) indexEntry.separator).updateMBR(mbr);
				indexEntry.hilbertEntryMBR = mbr; // LX
				return true;
			}
			return false;
		}

		/**
		 * @see BPlusTree.Node#redressUnderflow(xxl.core.indexStructures.BPlusTree.IndexEntry,
		 *      xxl.core.indexStructures.BPlusTree.Node, Stack, boolean)
		 */
		protected boolean redressUnderflow(IndexEntry indexEntry,
				BPlusTree.Node node, Stack path, boolean up) {
			Node parentNode = (Node) node(path);
			// The merge...
			MergeInfo mergeInfo = node.merge(indexEntry, path);
			// Post the changes to the parent Node...
			if (mergeInfo.isMerge()) {
				if (mergeInfo.isRightSide()) {
					ORSeparator separatorNode = (ORSeparator) mergeInfo
							.newSeparator();
					Rectangle mbrNode = (Rectangle) ((Descriptor) computeMBR(mergeInfo
							.node().entries())).clone();
					separatorNode.entryMBR = mbrNode;
					((IndexEntry) parentNode.getEntry(mergeInfo.index())).separator = separatorNode;
					mergeInfo.indexEntry().update(mergeInfo.node(), true);
					((IndexEntry) parentNode.remove(mergeInfo.index() + 1))
							.remove();
				} else {
					ORSeparator separatorSibling = (ORSeparator) mergeInfo
							.siblingNewSeparator();
					Rectangle mbrSiblingNode = (Rectangle) ((Descriptor) computeMBR(mergeInfo
							.siblingNode().entries())).clone();
					separatorSibling.entryMBR = mbrSiblingNode;
					((IndexEntry) parentNode.getEntry(mergeInfo.index() - 1)).separator = separatorSibling;
					mergeInfo.siblingIndexEntry().update(
							mergeInfo.siblingNode(), true);
					((IndexEntry) parentNode.remove(mergeInfo.index()))
							.remove();
				}
			} else {
				ORSeparator separatorNode = (ORSeparator) mergeInfo
						.newSeparator();
				Rectangle mbrNode = (Rectangle) ((Descriptor) computeMBR(mergeInfo
						.node().entries())).clone();
				separatorNode.entryMBR = mbrNode;
				ORSeparator separatorSibling = (ORSeparator) mergeInfo
						.siblingNewSeparator();
				Rectangle mbrSiblingNode = (Rectangle) ((Descriptor) computeMBR(mergeInfo
						.siblingNode().entries())).clone();
				separatorSibling.entryMBR = mbrSiblingNode;
				((IndexEntry) parentNode.getEntry(mergeInfo.index())).separator = separatorNode;
				int sibIndex = mergeInfo.isRightSide() ? mergeInfo.index() + 1
						: mergeInfo.index() - 1;
				((IndexEntry) parentNode.getEntry(sibIndex)).separator = separatorSibling;
				if (up) {
					mergeInfo.indexEntry().update(mergeInfo.node(), true);
					mergeInfo.siblingIndexEntry().update(
							mergeInfo.siblingNode(), true);
				}
			}
			return true;
		}

		/**
		 * <tt>Distribution</tt> is the class used to represent the distribution
		 * of entries of a node of the <tt>HilbertRTree</tt> into two partitions
		 * used for a split.
		 * 
		 * @see RTree.Node.Distribution
		 */
		protected class Distribution {
			/**
			 * Entries stored in this distribution.
			 */
			protected Object[] entries;

			/**
			 * Start index of the second part of the distribution.
			 */
			protected int secondStart;

			/**
			 * Bounding Rectangle of the first part of the distribution.
			 */
			protected Rectangle firstDescriptor;

			/**
			 * Bounding Rectangle of the first part of the distribution.
			 */
			protected Rectangle secondDescriptor;

			/**
			 * @param entries
			 *            an array containing all entries to be distributed
			 * @param secondStart
			 *            the start index of the second partition
			 * @param firstDescriptor
			 *            the descriptor for the first partition
			 * @param secondDescriptor
			 *            the descriptor for the second partition
			 */
			protected Distribution(Object[] entries, int secondStart,
					Rectangle firstDescriptor, Rectangle secondDescriptor) {
				this.entries = entries;
				this.secondStart = secondStart;
				this.firstDescriptor = firstDescriptor;
				this.secondDescriptor = secondDescriptor;
			}

			/**
			 * Returns one of the partitions of this distribution.
			 * 
			 * @param second
			 *            a <tt>boolean</tt> value determining if the second
			 *            partition should be returned
			 * @return the entries of the first (if <tt>second == false</tt>) or
			 *         of the second partition (if <tt>second == true</tt>)
			 */
			protected List entries(boolean second) {
				return Arrays.asList(entries).subList(second ? secondStart : 0,
						second ? entries.length : secondStart);
			}

			/**
			 * Returns a descriptor of one of the partitions of this
			 * distribution.
			 * 
			 * @param second
			 *            a <tt>boolean</tt> value determining if the descriptor
			 *            of second partition should be returned
			 * @return the descriptor of the first (if <tt>second == false</tt>)
			 *         or of the second partition (if <tt>second == true</tt>)
			 */
			protected Descriptor descriptor(boolean second) {
				return second ? secondDescriptor : firstDescriptor;
			}

			/**
			 * Returns the sum of the margins of the two partitions.
			 * 
			 * @return the sum of the margins of the two partitions
			 */
			protected double marginValue() {
				return firstDescriptor.margin() + secondDescriptor.margin();
			}

			/**
			 * Returns the overlap of the two partitions.
			 * 
			 * @return the overlap of the two partitions
			 */
			protected double overlapValue() {
				return firstDescriptor.overlap(secondDescriptor);
			}

			/**
			 * Returns the sum of the areas of the two partitions.
			 * 
			 * @return the sum of the areas of the two partitions
			 */
			protected double areaValue() {
				return firstDescriptor.area() + secondDescriptor.area();
			}

			/**
			 * Added by LX
			 * Return the area difference
			 */
			protected double areaDiffValue() {
				return Math.abs(firstDescriptor.area() - secondDescriptor.area());
			}
		}
	}

	/**
	 * Extends Separator with Rectangle field which represents the MBR.
	 * 
	 * @see Separator
	 */
	public static abstract class ORSeparator extends Separator {
		/** */
		protected Rectangle entryMBR;

		/**
		 * 
		 * @param separatorValue
		 * @param entryMBR
		 */
		public ORSeparator(Comparable separatorValue, Rectangle entryMBR) {
			super(separatorValue);
			this.entryMBR = entryMBR;
		}

		/**
		 * 
		 */
		public Rectangle getIndexEntryMBR() {
			return entryMBR;
		}

		/**
		 * 
		 * @param descriptor
		 * @return
		 */
		public boolean overlapsMBR(Descriptor descriptor) {
			if (descriptor instanceof ORSeparator) {
				return this.entryMBR.overlaps(((ORSeparator) descriptor)
						.getIndexEntryMBR());
			} else if (descriptor instanceof Rectangle) {
				return this.entryMBR.overlaps(descriptor);
			} else {
				throw new IllegalArgumentException();
			}
		}

		/**
		 * 
		 * @param descriptor
		 * @return
		 */
		public boolean containsMBR(Descriptor descriptor) {
			if (descriptor instanceof ORSeparator) {
				return this.entryMBR.contains(((ORSeparator) descriptor)
						.getIndexEntryMBR());
			} else if (descriptor instanceof Rectangle) {
				return this.entryMBR.contains(descriptor);
			} else {
				throw new IllegalArgumentException();
			}

		}

		/**
		 * 
		 * @param descriptor
		 * @return
		 */
		public void unionMBR(Descriptor descriptor) {
			if (descriptor instanceof ORSeparator) {
				this.entryMBR.union(((ORSeparator) descriptor)
						.getIndexEntryMBR());
			} else if (descriptor instanceof Rectangle) {
				this.entryMBR.union(descriptor);
			} else {
				throw new IllegalArgumentException();
			}
		}

		public void updateMBR(Rectangle mbr) {
			this.entryMBR = mbr;
		}

	}

	/**
	 * @see BPlusTree.KeyRange
	 */
	public static abstract class ORKeyRange extends KeyRange {

		protected Rectangle entryMBR;

		public ORKeyRange(Comparable min, Comparable max, Rectangle entryMBR) {
			super(min, max);
			this.entryMBR = entryMBR;
		}

		/**
		 * 
		 */
		public Descriptor getIndexEntryMBR() {
			return this.entryMBR;
		}

		/**
		 * 
		 */
		public void union(Descriptor descriptor) {
			if (!(descriptor instanceof ORSeparator))
				return;
			super.union(descriptor);
			if (descriptor instanceof ORKeyRange) {
				this.entryMBR.union(((ORKeyRange) descriptor)
						.getIndexEntryMBR());
			} else if (descriptor instanceof ORSeparator) {
				this.entryMBR.union(((ORSeparator) descriptor)
						.getIndexEntryMBR());
			}
		}

		/**
		 * 
		 */
		public void updateMBR(Rectangle mbr) {
			this.entryMBR = mbr;
		}

	}

	/**
	 * @see BPlusTree.NodeConverter
	 */
	public class NodeConverter extends BPlusTree.NodeConverter {
		@Override
		public Object read(DataInput dataInput, Object object)
				throws IOException {
			int level = dataInput.readInt();
			Node node = (Node) createNode(level);
			int number = dataInput.readInt();
			node.nodeId = dataInput.readLong(); // LX
			boolean readNext = dataInput.readBoolean();
			if (readNext) {
				node.nextNeighbor = (IndexEntry) createIndexEntry(level + 1);
				Container container = node.nextNeighbor.container();
				node.nextNeighbor.id = container.objectIdConverter().read(
						dataInput, null);
			} else {
				node.nextNeighbor = null;
			}

			for (int i = 0; i < number; i++) {
				Object entry;
				if (node.level == 0)
					entry = dataConverter.read(dataInput, null);
				else
					entry = readIndexEntry(dataInput, createIndexEntry(level));
				node.entries.add(i, entry);
			}
			return node;
		}

		public void write(DataOutput dataOutput, Object object)
				throws IOException {
			Node node = (Node) object;
			// 2x Integer
			dataOutput.writeInt(node.level);
			dataOutput.writeInt(node.number());
			dataOutput.writeLong(node.nodeId); // LX
			// Boolean
			dataOutput.writeBoolean(node.nextNeighbor != null);
			// ID
			if (node.nextNeighbor != null) {
				Container container = node.nextNeighbor.container();
				container.objectIdConverter().write(dataOutput,
						node.nextNeighbor.id());
			}
//
			// entries
			for (int i = 0; i < node.number(); i++) {
				if (node.level == 0)
					dataConverter.write(dataOutput, node.getEntry(i));
				else
					writeIndexEntry(dataOutput, node.getEntry(i));
			}
		}

		protected void writeIndexEntry(DataOutput dataOutput, Object entry)
				throws IOException {
			BPlusTree.IndexEntry indexEntry = (BPlusTree.IndexEntry) entry;
			Container container = indexEntry.container();
			container.objectIdConverter().write(dataOutput, indexEntry.id());
			keyConverter.write(dataOutput, indexEntry.separator().sepValue);
			ORSeparator orSeparator = (ORSeparator) indexEntry.separator();
			descriptorConverter.write(dataOutput, orSeparator.entryMBR);
		}

		protected Object readIndexEntry(DataInput dataInput, Object entry)
				throws IOException {
			BPlusTree.IndexEntry indexEntry = (BPlusTree.IndexEntry) entry; // LX
			Container container = indexEntry.container();
			indexEntry.id = container.objectIdConverter().read(dataInput, null);
			Comparable sepValue = (Comparable) keyConverter.read(dataInput,
					null);
			Rectangle mbr = (Rectangle) descriptorConverter.read(dataInput,
					null);
			indexEntry.separator = createORSeparator(sepValue, mbr);
			((IndexEntry) indexEntry).hilbertEntryMBR = mbr; // LX
			((IndexEntry) indexEntry).ORDescriptor = (ORSeparator) indexEntry.separator; // LX
			return indexEntry;
		}

		/**
		 * Computes the maximal size (in bytes) of an <tt>IndexEntry</tt>. It
		 * calls the method getIdSize() of the tree container. If the tree
		 * container has not been initialized a NullPointerException is thrown.
		 * 
		 * @return the maximal size of an <tt>IndexEntry</tt>
		 */
		protected int indexEntrySize() {
			return super.indexEntrySize() + dimension * 2 * 8;
		}

		protected int leafEntrySize() {
			return dataConverter.getMaxObjectSize();
		}
	}
}
