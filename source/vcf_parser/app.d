module vcf_parser.main;

import std.getopt;
import std.stdio;
import std.array : uninitializedArray;
import std.random;
import std.conv;
import std.string;
import std.process;

import vcf_parser.header;

// Main

void main(string[] args){
  string file_name;
  size_t BATCH_SIZE = 1000;

  getopt(args,
    // TODO: simple BIMBAM vs mean genotype BIMBAM
    "file"    , &file_name,
    "batch"   , &BATCH_SIZE
  );

  writeln("Processing file => ", file_name);
  auto vcfstats = getvcfcounts(file_name);
  writeln("Variants in VCF file: " ~ to!string(vcfstats["variants"]));
  writeln("Samples in VCF file: " ~ to!string(vcfstats["samples"]));

  File file;

  if(file_name[$-2..$] == "gz"){
    auto pipe = pipeShell("gunzip -c " ~ file_name);
    // TODO:
    // executeShell vs. pipeShell?
    // file below is of type pipe(?), if used want to strip white space
    // need to be string
    file = pipe.stdout;
  }
  else{
    file = File(file_name);
  }

  size_t index;
  int size = 0;

  size_t write_index = 0;
  string geno_file_name = "out.geno." ~ to!string(write_index) ~ ".bimbam";
  string pos_file_name = "out.pos.bimbam";

  // BIMBAM basic genotype format
  File genofile = File(geno_file_name, "w");
  genofile.write(to!string(vcfstats["samples"]) ~ "\n");
  genofile.write(to!string(vcfstats["variants"]) ~ "\n");

  // BIMBAM SNP location file
  File snp_pos = File(pos_file_name, "w");

  size_t batch_index;
  size_t line_index = 1;

  auto raw_genotype_matrix = IntMatrix2D(to!int(vcfstats["variants"]), to!int(vcfstats["samples"]));

  foreach(line; file.byLine){
    auto tokens = line.split("\t");
    if(tokens.length == 0 || tokens.length < 6){
      if(tokens[0].startsWith("##")) {
        continue;
      } else {
        writeln("File line: " ~ to!string(line_index) ~ "\n");
        throw new MisformatVCFLine("Header line not commented out properly or misformatted");
      }
    }

    
    if(tokens[0] == "#CHROM"){
      // Catch VCF header line
      char[][] samples = tokens[9..$];
      raw_genotype_matrix.setci(samples);

      genofile.write("IND");
      foreach(sample; tokens[9..$]){
        genofile.write("," ~ sample);
      }
      genofile.write("\n");
      continue;
    }

    // Chromosome/Contig ID
    auto chrom = to!string(tokens[0]);
    
    // Base pair position validation
    int pos;
    try {
      pos = to!int(tokens[1]);
    } catch (ConvException e) {
      throw new BPPOSNotNumeric("POS column must only contain integers.");
    }

    // SNP ID validation
    auto ids = tokens[2].split(";");
    string snp;
    if (ids[0] == ".") {
      snp = chrom ~ "_" ~ to!string(pos).idup;
    } else {
      snp = to!string(ids[0]);
    }
    
    if (ids.length > 1) {
      writeln("Only the first SNP id will be used in SNP position file. Please make sure this is the ID you want.");
    }

    // Reference and alterantive alleles
    auto reference = to!string(tokens[3]);
    /* 
      TODO: address multiple alternate alleles here:
        string.indexOf(',', alternate) for multiple allele testing , seperated alleles
    */
    auto alternate = to!string(tokens[4]); 

    // write snp position information
    snp_pos.write(to!string(snp) ~ "," ~ to!string(pos) ~ "," ~ to!string(chrom) ~ "\n");

    // write snp id to start line of genotype file
    genofile.write(to!string(snp));

    raw_genotype_matrix.setri(to!string(snp));

    // find GT index
    auto format = tokens[8].split(":");
    // TODO: remove this, per VCFv4.2 specification page 5, bottom of the page, 'GT' is first :
    int gt_index = 0;

    foreach(field; format) {
      if(field == "GT") {
        break;
      } else {
        ++gt_index;
      }
    }

    // Process sample GT data
    int[] gtrow;
    foreach(sample; tokens[9..$]) {
      auto sample_data = sample.split(':');
      auto gt_data = to!string(sample_data[gt_index].strip());
      auto score = gt_to_score(gt_data);

      gtrow ~= [to!int(score)];

      genofile.write("," ~ to!string(score));
    }

    raw_genotype_matrix.setrow(gtrow);
  
    // file line index for error reporting
    genofile.write("\n");
    ++line_index;

    // batch processing
    batch_index++;
    if(batch_index % BATCH_SIZE == 0){
      batch_index = 0;
      ++write_index;
      geno_file_name = "out.geno." ~ to!string(write_index) ~ ".bimbam";
      genofile = File(geno_file_name, "w");
    }
  }

  auto rgmfields = __traits(allMembers, typeof(raw_genotype_matrix));
  auto rgmvalues = raw_genotype_matrix.tupleof;

  foreach(rgmindex, value; rgmvalues) {
    writef("\n%-15s %s", rgmfields[rgmindex], value);
  }
}
