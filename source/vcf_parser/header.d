module vcf_parser.header;

import std.exception : basicExceptionCtors;
import std.conv;
import std.string;
import std.process;

// Custom exceptions

static class QualNotNumeric : Exception {
  mixin basicExceptionCtors;
}

static class BPPOSNotNumeric : Exception {
  mixin basicExceptionCtors;
}

static class MisformatVCFLine : Exception {
  mixin basicExceptionCtors;
}

// User defined structures

struct FloatMatrix2D {
  int rows;
  int cols;
  string[] row_ids;   
  string[] col_ids;        
  float[][] values;

  this(int r, int c) {
    this.rows = r;                                
    this.cols = c;
  }

  void setrow(int r, float[] vals) {
  	this.values[r] = vals;
  }

  void setci(char[][] ci) {
    string[] cids;
    string cid;

    foreach(cidcharr; ci[0..$]) {
      string current_cid = "";
      foreach (cich; cidcharr[0..$]) {
        current_cid ~= cich;
      }
      cids[$] = current_cid;
    }

    this.col_ids = cids;
  }

  void setri(string ri) {
  	this.row_ids[$] = ri;
  }
}

struct IntMatrix2D {
  int rows;
  int cols;
  string[] row_ids;   
  string[] col_ids;        
  int[][] values;

  this(int r, int c) {
    this.rows = r;                                
    this.cols = c;
  }

  void setrow(int r, int[] vals) {
  	this.values[r] = vals;
  }

  void setci(char[][] ci) {
    string[] cids;
    string cid;

    foreach(cidcharr; ci[0..$]) {
      string current_cid = "";
      foreach (cich; cidcharr[0..$]) {
        current_cid ~= cich;
      }
      cids ~= [current_cid];
    }

  	this.col_ids = cids;
  }

  void setri(string ri) {
  	this.row_ids[row_ids.length] = ri;
  }
}

// Auxillary functions

ulong[string] getvcfcounts(string vcffile) {
  /*
  // This is a rudementary non-gnu coreutils dependent
  // variant counter for vcf files
  auto f = File(vcffile);
  size_t c = 0; // count of lines
  size_t h = 0; // count of # 
  auto buffer = uninitializedArray!(ubyte[])(1024);
  foreach (chunk; f.byChunk(buffer)) {
    foreach(ch; chunk) {
      if(ch == '\n') {
        ++c;
      }
      if(ch == '#') {
        ++h;
      }
    }
  }

  size_t vars = c-(h/2)-1;
  return vars;
  */

  // This sample and variant counter for vcf files depends on
  // gnu core utils: grep, awk, wc, and cut
  auto grep = executeShell("grep ^[^##] " ~ vcffile ~ " | wc -l");
  auto goutput = strip(grep.output);
  auto vars = to!ulong(goutput);

  auto countsamps = executeShell("grep \"#CHROM\" " ~ vcffile ~ " | cut -f 10- | awk -F '\t' '{ print NF }'");
  auto coutput = strip(countsamps.output);
  auto samps = to!ulong(coutput);

  ulong[string] vcfinfo;
  vcfinfo["samples"] = samps;
  vcfinfo["variants"] = vars;

  return vcfinfo;
}

string gt_to_score(string gt) {
  if(gt.length > 1) {
    auto gts = gt.split("|");
    if (count(gts) < 2)
      gts = gt.split("/");

    if(gts[0] == "." || gts[1] == ".") {
      return "NA";
    }
    int par_one = to!int(gts[0]);
    int par_two = to!int(gts[1]);

    return to!string(par_one*(par_one+1)/2+par_two);
  } else if (gt.length == 1) {
    // haploids
    return "0";
  }
  return "0";
}

FloatMatrix2D compute_grm_from_rgm(IntMatrix2D rgm) {
  FloatMatrix2D grm;
  return grm;
}


