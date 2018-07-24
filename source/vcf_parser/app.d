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
  foreach(line; file.byLine){
    auto tokens = line.split("\t");
    foreach(token; tokens){
      writeln(token);
    }
    writeln("=================================");
  }

}
