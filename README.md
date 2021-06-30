# Fix_RNACentral_Taxonomy
Uses curated taxonomy and various software to clean up and deduplicate RNAcentral's database

## First Get The Curated Taxonomy
Follow these instructions to get the curated Taxonomy Database:
https://github.com/TealFurnholm/Universal-Taxonomy-Database

## Run the script
 - download the script in this repository
 - run: perl FIX_RNACENTAL_HEADERS.pl

## Script Details
 - The downloaded sequences are 
   - fixed: IUPAC -> "N"
   - sorted by length decending
   - deduplicated @ 100% identity, including containments
   - The headers of duplicates/containments are concatenated 
 - There are sometimes hundreds of identical sequences, but occasionally have a different RNA type annotation
   - The count of each RNA type is tracked
   - Only annotations within 75% of the max count are kept
   - eg. if 3 types: 7002 rRNA, 50 miRNA, 1 snoRNA - only the rRNA annotation is kept
 - There are often mixed "unclassified organism" and know species annotations (again hundreds of identical seqs)
   - If any know organisms, the LCA is derived from those
   - The LCA is made from only the top 25 taxon ids
     - The taxon ids are in order from longest to shortest
     - 25 is more than plenty to get to the lowest common ancestor, mostly lots and lots of strains of the same species
   - The header is the RNAC identifier of the longest-archaetype sequence.


