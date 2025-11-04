#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <zlib.h>
#include <filesystem>
#include <Rcpp.h>
using namespace Rcpp;

// This is a simple function using Rcpp that creates an R list
// containing a character vector and a numeric vector.
//
// Learn more about how to use Rcpp at:
//
//   https://www.rcpp.org/
//   https://adv-r.hadley.nz/rcpp.html
//
// and browse examples of code using Rcpp at:
//
//   https://gallery.rcpp.org/
//

// [[Rcpp::export]]
Rcpp::List getvcfinfo(CharacterVector vcf_name0) {
  std::string vcf_name_str = Rcpp::as<std::string>(vcf_name0[0]);
  const char* vcf_name = vcf_name_str.c_str();
  std::vector<std::string> chrom;
  std::vector<int> pos;
  std::vector<std::string> id;
  std::vector<std::string> ref;
  std::vector<std::string> alt;
  const char* c;
  std::string c_s;
  std::string id_s;
  std::string ref_s;
  std::string alt_s;
  std::string alt_s2;

  std::filesystem::path relativePath = vcf_name_str;
  std::filesystem::path absolutePath = std::filesystem::absolute(relativePath);
  std::string absolutePathStr = absolutePath.string();

  for (char &cp : absolutePathStr) {
    if (cp == '\\') {
      cp = '/';
    }
  }

  htsFile *inf = bcf_open(vcf_name, "r");
  if (inf == NULL) {
    Rcpp::Rcout << "File not opened" << std::endl;
    return NULL;
  }

  bcf_hdr_t *hdr = bcf_hdr_read(inf);

  bcf1_t *rec = bcf_init();

  Rcpp::CharacterVector sid;
  for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i)
    sid.push_back(hdr->samples[i]);

  int n = 0;
  while (bcf_read(inf, hdr, rec) == 0) {
    ++n;
    bcf_unpack(rec, BCF_UN_STR);
    bcf_unpack(rec, BCF_UN_INFO);

    //Process CHROM
    c = bcf_hdr_id2name(hdr, rec->rid);
    c_s = c;
    chrom.push_back(c_s);
    //Process POS
    pos.push_back(rec->pos+1);
    //Process ID
    id_s = rec->d.id;
    id.push_back(id_s);
    //process REF
    ref_s = rec->d.allele[0];
    ref.push_back(ref_s);
    //process ALT
    if(rec->n_allele > 1) {
      alt_s = rec->d.allele[1];
    } else {
      alt_s = "None";
    }
    //printf("%lu\n", (unsigned long)rec->n_allele);
    if(rec->n_allele > 2) {
      for (int i=2; i<rec->n_allele; ++i) {
        alt_s = alt_s + ' ';
        alt_s2 = rec->d.allele[i];
        alt_s = alt_s + alt_s2;
      }
    }
    alt.push_back(alt_s);
  }

  bcf_destroy(rec);
  bcf_hdr_destroy(hdr);

  if (hts_close(inf))
    Rcpp::Rcerr << "Error closing VCF file" << std::endl;

  Rcpp::DataFrame snps = Rcpp::DataFrame::create(
    Rcpp::Named("chromosome") = chrom,
    Rcpp::Named("location") = pos,
    Rcpp::Named("snpid") = id,
    Rcpp::Named("reference") = ref,
    Rcpp::Named("alternate") = alt);

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("filename") = absolutePathStr,
    Rcpp::Named("sid") = sid,
    Rcpp::Named("snps") = snps);

  return result;
}
