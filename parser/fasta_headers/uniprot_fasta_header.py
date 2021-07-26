import re

HEADER_REGEX = '(\w{2})\|(\w*)\|(\w*) (.*) (OS=.*)'

class UniprotFastaHeader:
  def doesMatch(self, headerString):
    return re.compile(HEADER_REGEX).match(headerString)

  def fromString(self, headerString):
    general_structure_parsed = re.search(HEADER_REGEX, headerString)

    key_value_regex = '(\w+)=(\w+)'
    key_pair_parsed = re.findall(key_value_regex, general_structure_parsed[5])

    key_pair_dictionary = {}
    for i in key_pair_parsed:
      iterable = iter(i)
      key = next(iterable)
      value = next(iterable)

      key_pair_dictionary[key] = value

    gene = None
    gene_exists = 'GN' in key_pair_dictionary
    if gene_exists == True:
      gene = key_pair_dictionary['GN']

    return {
      'proteinId': general_structure_parsed[2],
      'proteinName': general_structure_parsed[4],
      'gene': gene,
    }
