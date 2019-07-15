# vcf-parser
Parse VCF files to [r/qtl2](https://kbroman.org/qtl2/assets/vignettes/input_files.html) file format

# build
The makefile only works for the macos Dlang compiler (ldc2) currently.
```sh 
git clone https://github.com/rroutsong/vcf-parser.git
cd vcf-parser
make
```

# Run
The --output flag defaults to *vcf_csv_out.txt*

## Imputed
```sh
./build/vcf_parser -i --file=data/gatk_wgs.vcf --output=output.txt
```

## No imputation
```sh
./build/vcf_parser --file=data/gatk_wgs.vcf --output=output.txt
```
