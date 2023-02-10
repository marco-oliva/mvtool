# MVTool

Tool to extract lifting and sequences index from a set of VCF files to be used by `moni-align`


### Docker ###
PFP is available on docker:

```bash
docker pull moliva3/mvtool:latest
docker run moliva3/mvtool:latest mvtool --help
```

If using singularity:
```bash
singularity pull mvtool_sif docker://moliva3/mvtool:latest
./mvtool_sif mvtool --help
```

### Build ###

#### Dependencies ####

* Htslib
* OpenMP

#### Build Instructions ####

```
git clone https://github.com/marco-oliva/moni_vcf_tools.git
cd pfp
mkdir build && cd build
cmake ..
make
```

### Usage ###

```shell
MVTool
Usage: ./mvtool [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -v,--vcf TEXT ... REQUIRED  List of comma ',' separated vcf files. Assuming in genome order!
  -r,--ref TEXT ... REQUIRED  List of comma ',' separated reference files. Assuming in genome order!
  -w,--window-size UINT:INT in [3 - 200] REQUIRED
                              Sliding window size.
  -m,--max UINT               Max number of samples to analyze
  -o,--out-prefix TEXT REQUIRED
                              Output prefix
  -S,--samples TEXT           File containing the list of samples to parse
  -H,--haplotype TEXT         Haplotype: [1,2,12].
  --tmp-dir TEXT:DIR          Temporary files directory.

```
 
