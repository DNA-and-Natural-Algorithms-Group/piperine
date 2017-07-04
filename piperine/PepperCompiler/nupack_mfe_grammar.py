"""Grammar for NUPACK's mfe used by DNAfold"""

from pyparsing import *

Map = lambda f: (lambda s, t, l: map(f, l))

ParserElement.setDefaultWhitespaceChars(" \t\n")


nucleotides  = int_   =  Word(nums).setParseAction(Map(int))
energy       = float_ =  Word(nums+"-.").setParseAction(Map(float))
mfe_structure = Word(".()+")

comment = "%" + SkipTo( LineEnd() | StringEnd() , include=True)

document = StringStart() + Suppress(nucleotides) + energy + mfe_structure + Suppress(ZeroOrMore(int_ + int_)) + StringEnd()
document.ignore(comment)

parseFile = document.parseFile
