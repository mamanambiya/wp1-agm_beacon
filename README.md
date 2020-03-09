# Inserting CHILD WP3 data into Beacon-2.x

This script will insert the CHILD pseudodata and some sample variant data into the postgres database
backing the reference Beacon-2.x implementation.

* Start and run the Beacon-2.x reference deployment [here](https://github.com/EGA-archive/beacon-2.x)
* Insert the corresponding CHILD data with:

```
gunzip child.vcf
python3 ./insert_child.py --vcffile child.vcf fake_child.json
```

The fake_child data should be replaced with the WP3-mapped child pseudodata.
