import vcf, hgvs, vcf.utils

__author__ = 'Tamas Gutsohn'

class Assignment3:
    def __init__(self):
        ## Check if pyvcf is installed
        print("PyVCF version: %s" % vcf.VERSION)
        ## Check if hgvs is installed
        print("HGVS version: %s" % hgvs.__version__)

    def get_total_number_of_variants_mother(self, vcf):
        numbers = []
        for record in vcf:
            numbers.append(record.QUAL)
        totalnumber = len(numbers)
        print('\nTotal number of variants of mother: ')
        print(totalnumber)

    def get_total_number_of_variants_father(self, vcf):
        numbers = []
        for record in vcf:
            numbers.append(record.QUAL)
        totalnumber = len(numbers)
        print('\nTotal number of variants of father: ')
        print(totalnumber)

    def get_variants_shared_by_father_and_son(self, vcf1, vcf2):
        father_son_vcf = vcf.utils.walk_together(vcf1,vcf2)
        numbers = 0
        for record in father_son_vcf:
            if not record[0] is None and not record[1] is None:
                numbers += 1
        print('\nTotal number of shared variants of father and son: ')
        print(numbers)


    def get_variants_shared_by_mother_and_son(self, vcf1, vcf2):
        mother_son_vcf = vcf.utils.walk_together(vcf1,vcf2)
        numbers = 0
        for record in mother_son_vcf:
            if not record[0] is None and not record[1] is None:
                numbers += 1
        print('\nTotal number of shared variants of mother and son: ')
        print(numbers)

    def get_variants_shared_by_trio(self, vcf1, vcf2, vcf3):
        mother_father_son_vcf = vcf.utils.walk_together(vcf1,vcf2,vcf3)
        numbers = 0
        for record in mother_father_son_vcf:
            if not record[0] is None and not record[1] is None and not record[2] is None:
                numbers += 1
        print('\nTotal number of shared variants of mother, father and son: ')
        print(numbers)

    def merge_mother_father_son_into_one_vcf(self):
        merge_file = open("mother_father_son.vcf", "w")
        writer = vcf.Writer(merge_file, vcf_mother, "\n")
        for records in vcf.utils.walk_together(vcf_mother, vcf_father, vcf_son):
            for entry in records:
                if entry is not None:
                    writer.write_record(entry)
        print("Merge done!")

    def convert_first_variants_of_son_into_HGVS(self):
        '''
        Convert the first 100 variants identified in the son into the corresponding transcript HGVS.
        Each variant should be mapped to all corresponding transcripts. Pointer:
        - https://hgvs.readthedocs.io/en/master/examples/manuscript-example.html#project-genomic-variant-to-a-new-transcript
        :return:
        '''
        # import hgvs.dataproviders.uta
        # import hgvs.parser
        # from bioutils.assemblies import make_name_ac_map
        #
        # ## Connect to UTA
        # hdp = hgvs.dataproviders.uta.connect()
        #
        # ## Used to get the transcripts
        # assembly_mapper = hgvs.assemblymapper.AssemblyMapper(hdp)  # EasyVariantMapper before
        #
        # ## Used for parsing
        # hgvsparser = hgvs.parser.Parser()  # Parser
        #
        # ## Now for each variant
        #
        # ## Get chromosome mapping
        # refseq_nc_number = make_name_ac_map("GRCh37.p13")[record.CHROM[3:]]
        # ## Format: nc_number :g. position reference > alternative
        # genome_hgvs = "%s:g.%s%s>%s" % (refseq_nc_number, str(record.POS), str(record.REF), str(record.ALT[0]))
        #
        # ## Now parse the variant
        # ## http://hgvs.readthedocs.io/en/master/modules/io.html?highlight=parser_hgvs

    def print_summary(self):
        print("Print all results here")


if __name__ == '__main__':
    vcf_son = vcf.Reader(open('AmpliseqExome.20141120.NA24385.vcf', 'r'))
    vcf_mother = vcf.Reader(open('AmpliseqExome.20141120.NA24143.vcf', 'r'))
    vcf_father = vcf.Reader(open('AmpliseqExome.20141120.NA24149.vcf', 'r'))
    print("Assignment 3")
    assignment1 = Assignment3()
    assignment1.print_summary()
    assignment1.get_total_number_of_variants_mother(vcf_mother)
    assignment1.get_total_number_of_variants_father(vcf_father)
    assignment1.get_variants_shared_by_father_and_son(vcf_father, vcf_son)
    assignment1.get_variants_shared_by_mother_and_son(vcf_mother,vcf_son)
    assignment1.get_variants_shared_by_trio(vcf_mother,vcf_father,vcf_son)
    #assignment1.merge_mother_father_son_into_one_vcf()
