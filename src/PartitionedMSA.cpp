#include "PartitionedMSA.hpp"

using namespace std;

PartitionedMSA::PartitionedMSA(const NameList& taxon_names)
{
  set_taxon_names(taxon_names);
}

PartitionedMSA& PartitionedMSA::operator=(PartitionedMSA&& other)
{
  _part_list = std::move(other._part_list);
  _full_msa = std::move(other._full_msa);
  _taxon_names = std::move(other._taxon_names);
  _taxon_id_map = std::move(other._taxon_id_map);
   return *this;
}

void PartitionedMSA::set_taxon_names(const NameList& taxon_names)
{
  _taxon_names.assign(taxon_names.cbegin(), taxon_names.cend());
  for (size_t i = 0; i < _taxon_names.size(); ++i)
    _taxon_id_map[_taxon_names[i]] = i;

  assert(_taxon_names.size() == taxon_names.size() && _taxon_id_map.size() == taxon_names.size());
}

uintVector PartitionedMSA::get_site_part_assignment() const
{
  const size_t full_len = _full_msa.num_sites();

  uintVector spa(full_len);

  size_t p = 0;
  for (auto& pinfo: _part_list)
  {
    try
    {
      pinfo.mark_partition_sites(p+1, spa);
    }
    catch (MultiplePartitionForSiteException& e)
    {
      e.pinfo2(_part_list.at(spa[e.site()-1]-1));
      throw e;
    }
    p++;
  }

  /* check if all sites were assigned to partitions */
  MissingPartitionForSiteException e_unassinged;
  for (size_t i = 0; i < full_len; ++i)
  {
    if (!spa[i])
      e_unassinged.add_unassigned_site(i+1);
  }

  if (e_unassinged.count() > 0)
    throw e_unassinged;

  return spa;
}

const uintVector& PartitionedMSA::site_part_map() const
{
  if (_site_part_map.empty())
    _site_part_map = get_site_part_assignment();

  return _site_part_map;
}

size_t PartitionedMSA::full_msa_site(size_t index, size_t site) const
{
  if (part_count() == 1)
    return site;
  else
  {
    size_t cur_site = site;
    auto index_map = site_part_map();

    // TODO: this can be optimized
    for (size_t i = 0; i < index_map.size(); ++i)
    {
      if (index_map[i] == index+1)
      {
        if (!cur_site)
          return i;

        cur_site--;
      }
    }

    throw runtime_error("Site " + to_string(site+1) +
                        " not found in partition " +  to_string(index+1));
  }
}


void PartitionedMSA::full_msa(MSA&& msa)
{
  _full_msa = std::move(msa);

  set_taxon_names(_full_msa.labels());
}

void PartitionedMSA::split_msa()
{
  if (part_count() > 1)
  {
    const bool is_probabilistic = _full_msa.probabilistic();
    const size_t num_states = _full_msa.states();
    const auto& spm = site_part_map();
    const size_t num_taxa = _full_msa.size();
    const size_t full_length = _full_msa.length();

    /* split MSA into partitions */
    pll_msa_t ** part_msa_list =
        pllmod_msa_split(_full_msa.pll_msa(), spm.data(), part_count());

    /* if probabilistic, prepare to split probability vectors */
    std::vector<ProbVectorList> part_probs;
    std::vector<size_t> part_site_counts;

    if (is_probabilistic)
    {
      part_probs.resize(part_count());
      part_site_counts.resize(part_count(), 0);

      /* count sites per partition */
      for (size_t site = 0; site < full_length; ++site)
      {
        size_t part_idx = spm[site] - 1;
        part_site_counts[part_idx]++;
      }

      /* initialize probability vectors for each partition */
      for (size_t p = 0; p < part_count(); ++p)
      {
        part_probs[p].resize(num_taxa);
        for (size_t t = 0; t < num_taxa; ++t)
        {
          part_probs[p][t].resize(part_site_counts[p] * num_states);
        }
      }

      /* split probability vectors by site */
      std::vector<size_t> part_offsets(part_count(), 0);
      for (size_t site = 0; site < full_length; ++site)
      {
        size_t part_idx = spm[site] - 1;
        size_t part_site = part_offsets[part_idx];

        for (size_t t = 0; t < num_taxa; ++t)
        {
          auto src_it = _full_msa.probs(t, site);
          auto dst_it = part_probs[part_idx][t].begin() + part_site * num_states;
          std::copy(src_it, src_it + num_states, dst_it);
        }

        part_offsets[part_idx]++;
      }
    }

    for (size_t p = 0; p < part_count(); ++p)
    {
      part_msa(p, part_msa_list[p]);
      pll_msa_destroy(part_msa_list[p]);

      /* set probability data for probabilistic alignments */
      if (is_probabilistic)
      {
        _part_list[p].msa().set_probs(num_states, std::move(part_probs[p]));
      }
    }
    free(part_msa_list);
  }
  else
  {
    if (_part_list[0].range_string().empty())
      _part_list[0].range_string("1-" + to_string(_full_msa.num_sites()));
    part_msa(0, std::move(_full_msa));
  }
}

void PartitionedMSA::compress_patterns()
{
  for (PartitionInfo& pinfo: _part_list)
  {
    pinfo.compress_patterns();
  }
}

size_t PartitionedMSA::total_length() const
{
  size_t sum = 0;

  for (const auto& pinfo: _part_list)
  {
    sum += pinfo.msa().length();
  }

  return sum;
}

size_t PartitionedMSA::total_sites() const
{
  size_t sum = 0;

  for (const auto& pinfo: _part_list)
  {
    sum += pinfo.stats().site_count;
  }

  return sum;
}

size_t PartitionedMSA::total_patterns() const
{
  size_t sum = 0;

  for (const auto& pinfo: _part_list)
  {
    sum += pinfo.stats().pattern_count;
  }

  return sum;
}

size_t PartitionedMSA::total_free_model_params() const
{
  size_t sum = 0;

  for (const auto& pinfo: _part_list)
  {
    sum += pinfo.model().num_free_params();
  }

  return sum;
}

size_t PartitionedMSA::taxon_clv_size() const
{
  size_t clv_size = 0;

  for (const auto& pinfo: _part_list)
  {
    clv_size += pinfo.taxon_clv_size();
  }

  return clv_size;
}

void PartitionedMSA::set_model_empirical_params()
{
  for (PartitionInfo& pinfo: _part_list)
  {
    pinfo.set_model_empirical_params();
  }
}

std::ostream& operator<<(std::ostream& stream, const PartitionedMSA& part_msa)
{
  for (size_t p = 0; p < part_msa.part_count(); ++p)
  {
    const PartitionInfo& pinfo = part_msa.part_info(p);
    const auto pstats = pinfo.stats();
    stream << "Partition " << p << ": " << pinfo.name() << endl;
    stream << "Model: " << pinfo.model().to_string() << endl;
    if (pinfo.msa().num_patterns())
    {
      stream << "Alignment sites / patterns: " << pstats.site_count <<
          " / " << pstats.pattern_count << endl;
    }
    else
      stream << "Alignment sites: " << pinfo.msa().num_sites() << endl;

//    stream << fixed;
    stream << "Gaps: " << setprecision(2) << (pstats.gap_prop * 100) << " %" << endl;
    stream << "Invariant sites: " << setprecision(2) << (pstats.inv_prop * 100) << " %" << endl;
    stream << endl;
  }

  return stream;
}


