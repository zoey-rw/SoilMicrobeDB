#include "taxonomy.h"
#include "reports.h"

typedef std::vector<std::string> vs;
typedef std::unordered_map<uint64_t, uint64_t> taxon_counts_t;

void split(std::string& line, vs& tokens, const char* deli)
{
  char* c_line = const_cast<char*>(line.data());
  char* pch = strtok(c_line, deli);
  while (pch != NULL)
  {
    tokens.emplace_back(std::string(pch));
    pch = strtok(NULL, deli);
  }
}

int main(int argc, char** argv)
{
  kraken2::Taxonomy taxonomy(argv[1], false);
  std::ifstream kraken2_result(argv[2]);

  std::string buffer;
  uint64_t total_classified = 0;
  uint64_t total_unclassified = 0;
  std::unordered_map<uint64_t, uint> taxon;
  taxon_counts_t call_counts;

  for(uint64_t i = 0; i < taxonomy.node_count(); i++)
    taxon[taxonomy.nodes()[i].external_id] = i;

  while (kraken2_result.peek() != EOF)
  {
    getline(kraken2_result, buffer);

    if(buffer[0] == 'C')
    {
      vs tokens;
      split(buffer, tokens, "\t");
      call_counts[taxon[std::stoul(tokens[2])]]++;
      total_classified++;
    }
    else if(buffer[0] == 'U')
    {
      total_unclassified++;
    }
  }

  uint64_t total_sequences = total_unclassified + total_classified;
  kraken2::ReportKrakenStyle(argv[3], false, taxonomy,
    call_counts, total_sequences, total_unclassified);
  
  return 0;
}