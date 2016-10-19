from unittest import TestCase
import Assemble

class TestAssembler(TestCase):

    def test_read_fasta(self):
        pass #self.fail()

    def test_assemble(self):
        """
        Tests for assembly
        """

        # Small given test case
        dbg_assembler = Assemble.Assembler(fragments_fasta=file("Data/coding_challenge_data_set_example.txt"))
        assembled_sequence = dbg_assembler.assemble()
        self.assertEqual(assembled_sequence, "ATTAGACCTGCCGGAATAC")








