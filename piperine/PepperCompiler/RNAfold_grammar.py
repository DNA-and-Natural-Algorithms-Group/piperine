"""Grammar for Vienna RNAfold output format"""

from pyparsing import *

Map = lambda f: (lambda s, t, l: map(f, l))

ParserElement.setDefaultWhitespaceChars(" \t")


float_ = Word(nums+"-.").setParseAction(Map(float))

seq = Word("ATUCG"+"N") # Note that we use N for now as a dummy element
struct = Word(".()")

# (structure, dG) - ignoring input sequences
statement = Suppress(seq) + Suppress(LineEnd()) + \
            struct + Suppress("(") + float_ + Suppress(")") + Optional(Suppress(LineEnd()))

document = StringStart() + statement + StringEnd()

parseFile = document.parseFile
