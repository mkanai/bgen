
#include "reader.h"

namespace bgen {

CppBgenReader::CppBgenReader(std::string path, std::string sample_path, bool delay_parsing) {
  if (path != "/dev/stdin") {
    handle = new std::ifstream(path, std::ios::in | std::ios::binary);
  } else {
    is_stdin = true;
    handle = &std::cin;
  }
  if (handle->fail()) {
    throw std::invalid_argument("error reading from '" + path + "'");
  }
  header = Header(handle);
  if (header.has_sample_ids) {
    samples = Samples(handle, header.nsamples);
  } else if (sample_path.size() > 0) {
    samples = Samples(sample_path, header.nsamples);
  } else {
    samples = Samples(header.nsamples);
  }
  
  offset = header.offset + 4;
  if (!delay_parsing) {
    parse_all_variants();
  }
}

Variant CppBgenReader::next_var() {
  if (handle->eof()) {
    throw std::out_of_range("reached end of file");
  }
  Variant var(handle, offset, header.layout, header.compression, header.nsamples, is_stdin);
  offset = var.next_variant_offset;
  return var;
}

/// load all variants into memory at once
void CppBgenReader::parse_all_variants() {
  if (variants.size() == header.nvariants) {
    return;
  }
  offset = header.offset + 4;
  variants.clear();
  variants.resize(header.nvariants);
  int idx = 0;
  while (true) {
    try {
      variants[idx] = next_var();
      idx += 1;
    } catch (const std::out_of_range & e) {
      break;
    }
  }
  // finally reset the offset position to the first variant, so we can iterate
  // over variants more easily in python
  offset = header.offset + 4;
}

/// drop a subset of variants passed in by indexes
void CppBgenReader::drop_variants(std::vector<int> indices) {
  // sort indices in descending order, so dropping elemtns doesn't affect later items
  std::sort(indices.rbegin(), indices.rend());
  
  auto it = std::unique(indices.begin(), indices.end());
  if (it != indices.end()) {
    throw std::invalid_argument("can't drop variants with duplicate indices");
  }
  
  for (auto idx : indices) {
    variants[idx] = variants.back();
    variants.pop_back();
  }
  variants.shrink_to_fit();
  
  // and sort the variants again afterward
  std::sort(variants.begin(), variants.end(),
          [] (Variant const& a, Variant const& b) { return a.pos < b.pos; });
}

/// get all the IDs for the variants in the bgen file
std::vector<std::string> CppBgenReader::varids() {
  parse_all_variants();
  std::vector<std::string> varid(variants.size());
  for (std::uint32_t x=0; x<variants.size(); x++) {
    varid[x] = variants[x].varid;
  }
  return varid;
}

/// get all the rsIDs for the variants in the bgen file
std::vector<std::string> CppBgenReader::rsids() {
  parse_all_variants();
  std::vector<std::string> rsid(variants.size());
  for (std::uint32_t x=0; x<variants.size(); x++) {
    rsid[x] = variants[x].rsid;
  }
  return rsid;
}

/// get all the chroms for the variants in the bgen file
std::vector<std::string> CppBgenReader::chroms() {
  parse_all_variants();
  std::vector<std::string> chrom(variants.size());
  for (std::uint32_t x=0; x<variants.size(); x++) {
    chrom[x] = variants[x].chrom;
  }
  return chrom;
}

/// get all the positions for the variants in the bgen file
std::vector<std::uint32_t> CppBgenReader::positions() {
  parse_all_variants();
  std::vector<std::uint32_t> position(variants.size());
  for (std::uint32_t x=0; x<variants.size(); x++) {
    position[x] = variants[x].pos;
  }
  return position;
}

/// read variants at specific file offsets
std::vector<Variant> CppBgenReader::read_variants_at_offsets(const std::vector<std::uint64_t>& offsets) {
  if (is_stdin) {
    throw std::invalid_argument("cannot seek to offsets when reading from stdin");
  }

  std::vector<Variant> result;
  result.reserve(offsets.size());

  // Cast to ifstream to use seekg
  std::ifstream* file_handle = static_cast<std::ifstream*>(handle);

  for (std::uint64_t file_offset : offsets) {
    // Seek to the specific offset
    file_handle->clear(); // Clear any error flags
    file_handle->seekg(file_offset, std::ios::beg);

    if (file_handle->fail()) {
      throw std::runtime_error("failed to seek to offset " + std::to_string(file_offset));
    }

    // Read variant at this position
    try {
      Variant var(handle, file_offset, header.layout, header.compression, header.nsamples, is_stdin);
      result.push_back(std::move(var));
    } catch (const std::exception& e) {
      throw std::runtime_error("failed to read variant at offset " + std::to_string(file_offset) + ": " + e.what());
    }
  }

  return result;
}

} // namespace bgen
