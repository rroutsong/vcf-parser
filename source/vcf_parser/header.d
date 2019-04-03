module vcf_parser.header;

import std.exception : basicExceptionCtors;
import std.conv;
import std.string;
import std.process;
import std.stdio;
import std.math : sqrt;

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
  string[] row_allele;  
  string[] col_ids;        
  float[][] values;

  this(int r, int c) {
    this.rows = r;                                
    this.cols = c;
  }

  void setrow(float[] vals) {
  	this.values ~= vals;
  }

  void setci(string[] ci) {
    this.col_ids = ci;
  }

  void setri(string[] ri) {
  	this.row_ids = ri;
  }

  void setrallele(string[] ral) {
    this.row_allele = ral;
  }
}

struct IntMatrix2D {
  int rows;
  int cols;
  string[] row_ids;
  string[] row_allele;   
  string[] col_ids;        
  int[][] values;

  this(int r, int c) {
    this.rows = r;                                
    this.cols = c;
  }

  void setrow(int[] vals) {
  	this.values ~= vals;
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
  	this.row_ids ~= [ri];
  }

  void setrallele(string allelestring) {
    this.row_allele ~= [allelestring];
  }

  int[][] getvalues() {
    return this.values;
  }

  int[] getrowvaluebyindex(int r) {
    return this.values[r];
  }

  string getallelebyrowindex(int r) {
    return this.row_allele[r];
  }

  int getrowcount() {
    return this.rows;
  }

  int getcolcount() {
    return this.cols;
  }

  string[] getrowids() {
    return this.row_ids;
  }

  string[] getcolids() {
    return this.col_ids;
  }

  float getpj(int col) {
    float colmean = 0.0;
    int notmissing = 0;

    foreach (i;0..this.rows) {
      if(this.values[i][col] != -1) {
        colmean = colmean + to!float(this.values[i][col]);
        ++notmissing;
      }
    }

    return (colmean/to!float(notmissing))/2;
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

int gt_to_score(string gt) {
  if(gt.length > 1) {
    auto gts = gt.split("|");
    if (count(gts) < 2)
      gts = gt.split("/");

    if(gts[0] == "." || gts[1] == ".") {
      return -1;
    }
    int par_one = to!int(gts[0]);
    int par_two = to!int(gts[1]);

    return par_one*(par_one+1)/2+par_two;
  } else if (gt.length == 1) {
    // haploids
    return 0;
  }
  return 0;
}

FloatMatrix2D compute_grm_from_rgm(IntMatrix2D rgm) {
  // calculations based on hail formula:
  // https://hail.is/docs/0.2/methods/genetics.html#hail.methods.genetic_relatedness_matrix
  auto grm = FloatMatrix2D(rgm.getrowcount(), rgm.getcolcount());
  auto rowcount = rgm.getrowcount();
  auto cids = rgm.getcolids();
  auto rids = rgm.getrowids();

  grm.setri(rids);
  grm.setci(cids);

  foreach(i; 0..rowcount) {
    auto allelestring = rgm.getallelebyrowindex(i);
    auto rowdata = rgm.getrowvaluebyindex(i);
    float[] newrowdata;

    foreach(j; 0..rowdata.length) {
      auto alleles = allelestring.split(':');
      auto ref_al = alleles[0];
      float alt_count;

      if(indexOf(alleles[1], ',') == -1) {
        auto alt_al = alleles[1];
        alt_count = 1.00;
      } else {
        auto alt_al = alleles[1].split(',');
        alt_count = to!float(alt_al.length);
      }

      if (rowdata[j] == -1) {
        newrowdata ~= [0.00];
      } else {
        float pj = rgm.getpj(to!int(j));
        float m = to!float(rowcount);
        float val_ij = (alt_count-(2*pj))/(sqrt((2*pj)*(1-pj)*m));

        newrowdata ~= [val_ij];
      }
    }

    grm.setrow(newrowdata);
  }

  return grm;
}


