import re

HEADER_REGEX = '(\w*)\|(\w*) (.*)'

class UPSFastaHeader:
  def doesMatch(self, headerString):
    return re.compile(HEADER_REGEX).match(headerString)

  def fromString(self, headerString):
    general_structure_parsed = re.search(HEADER_REGEX, headerString)

    return {
      'proteinId': general_structure_parsed[1],
      'proteinName': general_structure_parsed[3],
      'gene': None,
    }
