class NotFoundFastaHeader:
  def doesMatch(self, headerString):
    return True

  def fromString(self, headerString):
    return {
      'proteinId': headerString,
      'proteinName': None,
      'gene': None,
    }
