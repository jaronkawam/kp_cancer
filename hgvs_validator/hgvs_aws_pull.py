#starter: https://github.com/biocommons/hgvs/blob/e6dbd1e7933ab96a1fc0a8ace9b46da62d244eae/src/hgvs/extras/babelfish.py#L39
import os

from bioutils.assemblies import make_ac_name_map, make_name_ac_map
from bioutils.sequences import reverse_complement
import hgvs
import hgvs.normalizer
from hgvs.edit import NARefAlt
from hgvs.location import Interval, SimplePosition
from hgvs.normalizer import Normalizer
from hgvs.posedit import PosEdit
from hgvs.sequencevariant import SequenceVariant
import hgvs.dataproviders.uta

import copy
import logging
import re
import bioutils.sequences
import ometa.runtime
import parsley
from hgvs import edit, enums, hgvsposition, location, posedit, sequencevariant
from hgvs.exceptions import HGVSParseError
from src.hgvs.generated.hgvs_grammar import createParserClass

hp = hgvs.parser.Parser()
hdp = hgvs.dataproviders.uta.connect()
# vm = hgvs.variantmapper.VariantMapper(hdp, assembly_name="GRCh38")

def _as_interbase(posedit):
    if posedit.edit.type == "ins":
        # ins coordinates (only) exclude left position
        start_i = posedit.pos.start.base
        end_i = posedit.pos.end.base - 1
    else:
        start_i = posedit.pos.start.base - 1
        end_i = posedit.pos.end.base
    return (start_i, end_i)

class Babelfish:
    def __init__(self, hdp, assembly_name):
        self.assembly_name = assembly_name
        self.hdp = hdp
        self.hn = hgvs.normalizer.Normalizer(
            hdp, cross_boundaries=False, shuffle_direction=5, validate=False
        )
        self.ac_to_name_map = make_ac_name_map(assembly_name)
        self.name_to_ac_map = make_name_ac_map(assembly_name)
        # We need to accept accessions as chromosome names, so add them pointing at themselves
        self.name_to_ac_map.update({ac: ac for ac in self.name_to_ac_map.values()})

    def vcf_to_g_hgvs(self, chrom, position, ref, alt):
        # VCF spec https://samtools.github.io/hts-specs/VCFv4.1.pdf
        # says for REF/ALT "Each base must be one of A,C,G,T,N (case insensitive)"
        ref = ref.upper()
        alt = alt.upper()

        ac = self.name_to_ac_map[chrom]

        if ref != alt:
            # Strip common prefix
            if len(alt) > 1 and len(ref) > 1:
                pfx = os.path.commonprefix([ref, alt])
                lp = len(pfx)
                if lp > 0:
                    ref = ref[lp:]
                    alt = alt[lp:]
                    position += lp
            elif alt == ".":
                alt = ref

        if ref == "":  # Insert
            # Insert uses coordinates around the insert point.
            start = position - 1
            end = position
        else:
            start = position
            end = position + len(ref) - 1

        var_g = SequenceVariant(
            ac=ac,
            type="g",
            posedit=PosEdit(
                Interval(start=SimplePosition(start), end=SimplePosition(end), uncertain=False),
                NARefAlt(ref=ref or None, alt=alt or None, uncertain=False),
            ),
        )
        n = Normalizer(self.hdp)
        return n.normalize(var_g)

    def hgvs_to_vcf(self, var_g):
        """**EXPERIMENTAL**

        converts a single hgvs allele to (chr, pos, ref, alt) using
        the given assembly_name. The chr name uses official chromosome
        name (i.e., without a "chr" prefix).
        """

        if var_g.type != "g":
            raise RuntimeError("Expected g. variant, got {var_g}".format(var_g=var_g))
        # if var_g.type == "c":
        #     var_g = vm.c_to_g(var_g)

        vleft = self.hn.normalize(var_g)

        (start_i, end_i) = _as_interbase(vleft.posedit)

        chrom = self.ac_to_name_map[vleft.ac]

        typ = vleft.posedit.edit.type

        if typ == "dup":
            start_i -= 1
            alt = self.hdp.seqfetcher.fetch_seq(vleft.ac, start_i, end_i)
            ref = alt[0]
        elif typ == "inv":
            ref = vleft.posedit.edit.ref
            alt = reverse_complement(ref)
        else:
            alt = vleft.posedit.edit.alt or ""

            if typ in ("del", "ins"):  # Left anchored
                start_i -= 1
                ref = self.hdp.seqfetcher.fetch_seq(vleft.ac, start_i, end_i)
                alt = ref[0] + alt
            else:
                ref = vleft.posedit.edit.ref
                if ref == alt:
                    alt = "."
        return chrom, start_i + 1, ref, alt, typ


class Parser:
    """Provides comprehensive parsing of HGVS variant strings (*i.e.*,
    variants represented according to the Human Genome Variation
    Society recommendations) into Python representations.  The class
    wraps a Parsing Expression Grammar, exposing rules of that grammar
    as methods (prefixed with `parse_`) that parse an input string
    according to the rule.  The class exposes all rules, so that it's
    possible to parse both full variant representations as well as
    components, like so:

    >>> hp = Parser()
    >>> v = hp.parse_hgvs_variant("NM_01234.5:c.22+1A>T")
    >>> v
    SequenceVariant(ac=NM_01234.5, type=c, posedit=22+1A>T, gene=None)
    >>> v.posedit.pos
    BaseOffsetInterval(start=22+1, end=22+1, uncertain=False)
    >>> i = hp.parse_c_interval("22+1")
    >>> i
    BaseOffsetInterval(start=22+1, end=22+1, uncertain=False)

    The `parse_hgvs_variant` and `parse_c_interval` methods correspond
    to the `hgvs_variant` and `c_interval rules` in the grammar,
    respectively.

    As a convenience, the Parser provides the `parse` method as a
    shorthand for `parse_hgvs_variant`:
    >>> v = hp.parse("NM_01234.5:c.22+1A>T")
    >>> v
    SequenceVariant(ac=NM_01234.5, type=c, posedit=22+1A>T, gene=None)

    Because the methods are generated on-the-fly and depend on the
    grammar that is loaded at runtime, a full list of methods is not
    available in the documentation.  However, the list of
    rules/methods is available via the `rules` instance variable.

    A few notable methods are listed below:

    `parse_hgvs_variant()` parses any valid HGVS string supported by the grammar.

      >>> hp.parse_hgvs_variant("NM_01234.5:c.22+1A>T")
      SequenceVariant(ac=NM_01234.5, type=c, posedit=22+1A>T, gene=None)
      >>> hp.parse_hgvs_variant("NP_012345.6:p.Ala22Trp")
      SequenceVariant(ac=NP_012345.6, type=p, posedit=Ala22Trp, gene=None)

    The `hgvs_variant` rule iteratively attempts parsing using the
    major classes of HGVS variants. For slight improvements in
    efficiency, those rules may be invoked directly:

      >>> hp.parse_p_variant("NP_012345.6:p.Ala22Trp")
      SequenceVariant(ac=NP_012345.6, type=p, posedit=Ala22Trp, gene=None)

    Similarly, components of the underlying structure may be parsed
    directly as well:

      >>> hp.parse_c_posedit("22+1A>T")
      PosEdit(pos=22+1, edit=A>T, uncertain=False)
      >>> hp.parse_c_interval("22+1")
      BaseOffsetInterval(start=22+1, end=22+1, uncertain=False)

    """

    def __init__(self, grammar_fn=None, expose_all_rules=False):
        bindings = {"hgvs": hgvs, "bioutils": bioutils, "copy": copy}
        if grammar_fn is None:
            self._grammar = parsley.wrapGrammar(
                createParserClass(ometa.runtime.OMetaGrammarBase, bindings)
            )
        else:
            # Still allow other grammars if you want
            with open(grammar_fn, "r") as grammar_file:
                self._grammar = parsley.makeGrammar(grammar_file.read(), bindings)
        self._logger = logging.getLogger(__name__)
        self._expose_rule_functions(expose_all_rules)

    def parse(self, v):
        """parse HGVS variant `v`, returning a SequenceVariant

        :param str v: an HGVS-formatted variant as a string
        :rtype: SequenceVariant

        """
        return self.parse_hgvs_variant(v)

    def _expose_rule_functions(self, expose_all_rules=False):
        """add parse functions for public grammar rules

        Defines a function for each public grammar rule, based on
        introspecting the grammar. For example, the `c_interval` rule
        is exposed as a method `parse_c_interval` and used like this::

          Parser.parse_c_interval('26+2_57-3') -> Interval(...)

        """

        def make_parse_rule_function(rule_name):
            "builds a wrapper function that parses a string with the specified rule"

            def rule_fxn(s):
                try:
                    return self._grammar(s).__getattr__(rule_name)()
                except ometa.runtime.ParseError as exc:
                    raise HGVSParseError(
                        "{s}: char {exc.position}: {reason}".format(
                            s=s, exc=exc, reason=exc.formatReason()
                        )
                    )

            rule_fxn.__doc__ = "parse string s using `%s' rule" % rule_name
            return rule_fxn

        exposed_rule_re = re.compile(
            r"hgvs_(variant|position)|(c|g|m|n|p|r)"
            r"_(edit|hgvs_position|interval|pos|posedit|variant)"
        )
        exposed_rules = [
            m.replace("rule_", "")
            for m in dir(self._grammar._grammarClass)
            if m.startswith("rule_")
        ]
        if not expose_all_rules:
            exposed_rules = [
                rule_name for rule_name in exposed_rules if exposed_rule_re.match(rule_name)
            ]
        for rule_name in exposed_rules:
            att_name = "parse_" + rule_name
            rule_fxn = make_parse_rule_function(rule_name)
            self.__setattr__(att_name, rule_fxn)
        self._logger.debug(
            "Exposed {n} rules ({rules})".format(
                n=len(exposed_rules), rules=", ".join(exposed_rules)
            )
        )


if __name__ == "__main__":
    hdp = hgvs.dataproviders.uta.connect()
    hp = Parser()
    test = Babelfish(hdp, 'GRCh38')
    print("vcf to hgvs: " + str(Babelfish.vcf_to_g_hgvs(test, chrom='20', position=10001019, ref='T', alt='G')))
    print("hgvs to vcf: " + str(Babelfish.hgvs_to_vcf(test, hp.parse_hgvs_variant('NC_000020.11:g.10001019T>G'))))
    print(Babelfish.hgvs_to_vcf(test, hp.parse_hgvs_variant('NC_000004.12:c.941A>G')))
