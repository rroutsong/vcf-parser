module vcf_parser.main;

import std.getopt;
import std.stdio;
import std.random;
import std.conv;
import std.string;
import std.process;

void main(string[] args){
  string file_name;

  getopt(args,
    "file"    , &file_name,
  );

  writeln("Processing file => ", file_name);

  File file = File(file_name);
  //auto pipe = pipeShell("gunzip -c " ~ file_name);
  //File file = pipe.stdout;

  size_t index;
  int size = 0;
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
    auto qual = to!double(tokens[5]);
    auto filter = to!string(tokens[6]);
    auto info = tokens[7].split(";");

    writeln("chrom => ", chrom);
    writeln("pos=> " , pos);
    //writeln("ids => ", ids);
    writeln("id => ", id);
    writeln("reference => ", reference);
    writeln("alternate => ", alternate);
    writeln("qual => ", qual);
    writeln("filter => ", filter);

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
    writeln("info =>", map_info_to_val);

    auto format = tokens[8].split(":");


    writeln("format => ", format);
    writeln("samples => ");
    foreach(sample; tokens[9..$]){
      string[string] sample_store;
      auto sample_vals = sample.split(":");
      size_t index2 = 0;
      foreach(key; format){
        sample_store[to!string(key)] = to!string(sample_vals[index]);
        index2++;
      }
      writeln(sample_store);
    }
    writeln("=================================");

  }

}
