# Experimental Evaluation and Investigation of Waves of Misery in R-trees (Experiments, Analysis & Benchmark)
 
## Description
Waves of misery is a phenomenon observed mainly in insertion-heavy disk-based 
tree-like indexes, include B-tree and R-tree

## Prerequisites
- java

## Usage
### Data generation
Generate point files for bulk loading, uniform or non-uniform query files with 
different range sizes, transform a csv file into a point file etc. 
A sample file is provided as `SampleFileWriter.properties`. Note that fixed 
selectivity queries for real data are generated separately.

Assuming in the same directory as `README.md` 
```java
java -classpath target/test-classes:target/classes xxl.core.indexStructures.rtrees.WritePointOrQueryFile SampleFileWriter.properties
```

### Experiment 1
Most of the experiments are run this way. An example parameter file is provided as
`SampleExperimentParams.properties`. 

Assuming in the same directory as `README.md`.
```java
java -classpath target/test-classes:target/classes xxl.core.indexStructures.rtrees.WavesMisery SampleExperimentParams.properties
```

### Experiment 2
Experiment with bulk loading method STR. An example parameter file is provided as
`SampleExperimentParams.properties`, which is the same as the above.

Assuming in the same directory as `README.md`.
```java
java -classpath target/test-classes:target/classes xxl.core.indexStructures.rtrees.WavesMiserySTR SampleExperimentParams.properties
```

### Experiment 3
Experiment to measure tail latency. An example parameter file is provided as
`SampleExperimentParams.properties`. Note that `Query_file` is only the directory 
for all the 6 range sizes. The default query files' names in the directory are: 
- `query_0.05_num_1k`
- `query_0.01_num_1k`
- `query_0.001_num_1k`
- `query_0.0001_num_1k`
- `query_0.00001_num_1k`
- `query_0.000001_num_1k`

Assuming in the same directory as `README.md`.
```java
java -classpath target/test-classes:target/classes xxl.core.indexStructures.rtrees.WavesMiseryTailLatency SampleExperimentParams.properties
```

### Experiment 4
Generate fixed range queries as well as issue these queries. An example parameter 
file is provided as `SampleExperimentParams.properties`. There will be one query file
per batch. Note that `kNN_query` should
be true to generate these queries and when issuing `Query_file` is the directory
that contains the queries.

Assuming in the same directory as `README.md`.
```java
java -classpath target/test-classes:target/classes xxl.core.indexStructures.rtrees.WavesMisery SampleExperimentParams.properties
```
