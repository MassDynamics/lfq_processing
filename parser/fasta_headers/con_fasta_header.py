import re

HEADER_REGEX = '(CON__.*)'

class CONFastaHeader:
  def doesMatch(self, headerString):
    return re.compile(HEADER_REGEX).match(headerString)

  def fromString(self, headerString):
    general_structure_parsed = re.search(HEADER_REGEX, headerString)

    return {
      'proteinId': general_structure_parsed[1],
      'proteinName': None,
      'gene': None,
    }
