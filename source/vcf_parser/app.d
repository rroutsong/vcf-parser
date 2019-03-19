module vcf_parser.main;

import std.getopt;
import std.stdio;
import std.random;
import std.conv;
import std.string;
import std.process;
import std.exception;

static class QualNotNumeric : Expcetion {
  mixin basicExceptionCtors;
}

void main(string[] args){
  string file_name;
  size_t BATCH_SIZE = 1000;

  getopt(args,
    "file"    , &file_name,
    "batch"   , &BATCH_SIZE
  );

  writeln("Processing file => ", file_name);

  File file;

  if(file_name[$-2..$] == "gz"){
    auto pipe = pipeShell("gunzip -c " ~ file_name);
    file = pipe.stdout;
  }
  else{
    file = File(file_name);
  }

  size_t index;
  int size = 0;

  size_t write_index = 0;
  string out_file_name = "out." ~ to!string(write_index) ~ ".bimbam";
  File outfile = File(out_file_name, "w");

  size_t batch_index;
  foreach(line; file.byLine){
    auto tokens = line.split("\t");
    if(tokens.length == 0 || tokens.length < 6){continue;}
    if(tokens[0] == "#CHROM"){continue;}
    auto chrom = to!string(tokens[0]);
    auto pos = to!int(tokens[1]);
    auto ids = tokens[2].split(";");
    auto id = ids[0];
    auto reference = to!string(tokens[3]);
    auto alternate = to!string(tokens[4]);
    double qual;
    try {
      qual = to!double(tokens[5]);
    } catch (ConvException e) {
      throw new QualNotNumeric("QUAL column must only contain numeric data. No symbols or strings.");
    }
    auto filter = to!string(tokens[6]);
    auto info = tokens[7].split(";");

    outfile.write("chrom => ", chrom);
    outfile.write("pos=> " , pos);
    //outfile.write("ids => ", ids);
    outfile.write("id => ", id);
    outfile.write("reference => ", reference);
    outfile.write("alternate => ", alternate);
    outfile.write("qual => ", qual);
    outfile.write("filter => ", filter);
    outfile.write("\n");

    double[string] map_info_to_val;
    foreach(item; info){
      auto dict =  item.split("=");
      double value;
      if(dict.length == 2){
        if(dict[1].isNumeric )
          value = to!double(dict[1]);
        else                                 // TODO: check for set and missing value;
          value = 0;
      }else{
        value = 0;
      }
      map_info_to_val[to!string(dict[0])] = value;
    }
    outfile.write("info =>", map_info_to_val);

    auto format = tokens[8].split(":");


    outfile.write("format => ", format);
    outfile.write("samples => ");
    foreach(sample; tokens[9..$]){
      string[string] sample_store;
      auto sample_vals = sample.split(":");
      size_t index2 = 0;
      foreach(key; format){
        sample_store[to!string(key)] = to!string(sample_vals[index]);
        index2++;
      }
      outfile.write(sample_store);
    }
    outfile.write("=================================");
    batch_index++;
    if(batch_index % BATCH_SIZE == 0){
      batch_index = 0;
      write_index++;
      out_file_name = "out." ~ to!string(write_index) ~ ".bimbam";
      outfile = File(out_file_name, "w");
    }
  }

}
